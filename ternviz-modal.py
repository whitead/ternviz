import gzip
import io
import math
import os
import shutil
import tarfile
import uuid
import warnings
from glob import glob
from pathlib import Path
from typing import List

import modal


def _resolve_vmd_archive() -> Path:
    candidates = []

    configured = os.environ.get("TERNVIZ_VMD_ARCHIVE_PATH")
    if configured:
        candidates.append(Path(configured).expanduser())

    candidates.append(Path.home() / "Downloads" / "vmd-2.0.0.bin.LINUXAMD64.tar")

    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()

    raise RuntimeError(
        "Set TERNVIZ_VMD_ARCHIVE_PATH to a local VMD 2 Linux tarball before deploying. "
        "Modal will upload that archive into the image build with add_local_file(copy=True)."
    )

VMD_REMOTE_ARCHIVE = "/tmp/vmd.tar"
VMD_ARCHIVE_PATH = _resolve_vmd_archive() if modal.is_local() else None
GPU_RENDERER = "TachyonLOptiXInternal"
GPU_RENDER_FORMAT = "PPM"
GPU_FRAME_EXT = "ppm"
MAX_RENDER_CONTAINERS = 4
TARGET_FRAMES_PER_RENDER_JOB = 30
MOVIE_RENDER_CPU = 16

app = modal.App("ternviz-pdb")
volume = modal.Volume.from_name("ternviz-scratch", create_if_missing=True)
volumes = {"/var/tmp": volume}

image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install(
        "ffmpeg",
        "fonts-dejavu-core",
        "libgl1",
        "libglu1-mesa",
        "libgomp1",
        "libice6",
        "libatk1.0-0",
        "libgtk-3-0",
        "libsm6",
        "libx11-6",
        "libxcursor1",
        "libxext6",
        "libxft2",
        "libxi6",
        "libxinerama1",
        "libxmu6",
        "libxrender1",
        "libxt6",
        "libxxf86vm1",
        "make",
    )
    .pip_install("fastapi")
    .pip_install_from_pyproject("pyproject.toml")
)

if VMD_ARCHIVE_PATH is not None:
    image = image.add_local_file(VMD_ARCHIVE_PATH, VMD_REMOTE_ARCHIVE, copy=True)

image = image.run_commands(
    "mkdir -p /opt/vmd",
    f"tar -xf {VMD_REMOTE_ARCHIVE} -C /opt/vmd",
    "cd /opt/vmd/vmd-* && ./configure",
    "cd /opt/vmd/vmd-*/src && make install",
    "cp /opt/vmd/vmd-*/lib/redistrib/lib_LINUXAMD64/liboptix*.so* /usr/local/lib/vmd/",
    "test -x /usr/local/bin/vmd",
).add_local_dir("ternviz", "/root/ternviz", copy=True)

with image.imports():
    import MDAnalysis as mda
    from MDAnalysis.analysis import align

    import ternviz


def get_backbone_selection(
    universe: mda.Universe, reference_universe: mda.Universe
) -> str:
    largest_chain = max(universe.segments, key=lambda s: len(s.residues))
    return f"segid {largest_chain.segid} and backbone"


def align_pdbs_to_reference(pdb_files: List[bytes], reference_pdb: bytes) -> List[bytes]:
    ref = mda.Universe(io.StringIO(reference_pdb.decode("utf-8")), format="PDB")
    aligned_pdb_strings = []

    for pdb_data in pdb_files:
        u = mda.Universe(io.StringIO(pdb_data.decode("utf-8")), format="PDB")
        ref_selection = get_backbone_selection(u, ref)

        alignment = align.AlignTraj(u, ref, select=ref_selection, in_memory=True)
        alignment.run()

        buffer = io.StringIO()
        with mda.Writer(buffer, n_atoms=u.atoms.n_atoms, format="PDB") as writer:
            writer.write(u.atoms)
            aligned_pdb_strings.append(buffer.getvalue().encode("utf-8"))

    return aligned_pdb_strings


def expand_render_files(
    render_files: List[bytes],
    target_total_frames: int,
    max_frames_per_job: int = TARGET_FRAMES_PER_RENDER_JOB,
    max_render_jobs: int = MAX_RENDER_CONTAINERS,
) -> List[bytes]:
    desired_render_jobs = max(
        len(render_files),
        min(max_render_jobs, math.ceil(target_total_frames / max_frames_per_job)),
    )
    if len(render_files) >= desired_render_jobs:
        return render_files
    return [
        render_files[i % len(render_files)] for i in range(desired_render_jobs)
    ]


async def extract_gz_content(file_content: bytes) -> bytes:
    with gzip.open(io.BytesIO(file_content), "rb") as gz:
        return gz.read()


async def extract_tar_content(file_content: bytes) -> List[tuple]:
    pdb_files = []
    with tarfile.open(fileobj=io.BytesIO(file_content), mode="r:*") as tar:
        for member in tar.getmembers():
            if member.isfile() and member.name.endswith(".pdb"):
                extracted_file = tar.extractfile(member)
                if extracted_file:
                    pdb_files.append((member.name, extracted_file.read()))
    return pdb_files


async def process_upload_file(upload_file) -> List[bytes]:
    pdb_files = {}
    file_content = await upload_file.read()

    file_name = upload_file.filename
    if file_name.endswith(".pdb"):
        pdb_files[file_name] = file_content
    elif file_name.endswith(".gz") and not file_name.endswith(".tar.gz"):
        extracted_content = await extract_gz_content(file_content)
        if file_name[:-3].endswith(".pdb"):
            pdb_files[file_name[:-3]] = extracted_content
    elif file_name.endswith(".tar") or file_name.endswith(".tar.gz"):
        extracted_pdb_files = await extract_tar_content(file_content)
        for fname, content in extracted_pdb_files:
            pdb_files[fname] = content

    return [pdb_files[fname] for fname in sorted(pdb_files)]


def build_web_app():
    from fastapi import Depends, FastAPI, HTTPException, UploadFile, status
    from fastapi.responses import Response
    from fastapi.security import HTTPAuthorizationCredentials, HTTPBearer

    web_app = FastAPI()
    auth_scheme = HTTPBearer()

    @web_app.post("/render/")
    async def render(
        pdb_file: UploadFile,
        message: str = "",
        render_first_struct: bool = True,
        target_total_frames: int = 360,
        token: HTTPAuthorizationCredentials = Depends(auth_scheme),
    ):
        if token.credentials != os.environ["AUTH_TOKEN"]:
            raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail="Incorrect bearer token",
                headers={"WWW-Authenticate": "Bearer"},
            )

        if not pdb_file.filename.endswith((".pdb", ".gz", ".tar", ".tar.gz")):
            raise HTTPException(status_code=400, detail="Unsupported file type")

        pdb_files = await process_upload_file(pdb_file)
        if not pdb_files:
            raise HTTPException(
                status_code=400, detail="No PDB files found in the uploaded file"
            )

        render_files = pdb_files
        if len(pdb_files) > 1:
            with warnings.catch_warnings():
                render_files = align_pdbs_to_reference(pdb_files, pdb_files[0])
            if not render_first_struct and len(render_files) > 1:
                render_files = render_files[1:]

        coloring = "Structure" if len(render_files) == 1 else "Chain"
        target_total_frames = max(target_total_frames, len(render_files))
        job_id = str(uuid.uuid4()).split("-")[0]
        render_files = expand_render_files(render_files, target_total_frames)

        n_frames = len(render_files) * math.ceil(
            target_total_frames / len(render_files)
        )
        rotation = 360.0 / n_frames
        batch_size = n_frames // len(render_files)

        print(
            "Will render",
            n_frames,
            "frames in",
            len(render_files),
            "batches of",
            batch_size,
            "frames each.",
        )

        def render_jobs():
            index = 0
            for pdb_bytes in render_files:
                yield job_id, pdb_bytes, index, index + batch_size, rotation, coloring
                index += batch_size

        list(pdb_render.starmap(render_jobs()))

        blob = movie_render.remote(message, job_id)
        return Response(
            content=blob,
            media_type="media/mp4",
            headers={"Content-Disposition": "attachment; filename=ternviz.mp4"},
        )

    return web_app


@app.function(
    cpu=4,
    gpu="l4",
    image=image,
    volumes=volumes,
    max_containers=MAX_RENDER_CONTAINERS,
)
def pdb_render(
    pdb_id: str,
    pdb_file: bytes,
    start: int,
    stop: int,
    rotation: float,
    color_method: str,
):
    with open("/tmp/pdb.pdb", "wb") as f:
        f.write(pdb_file)
    ternviz.render(
        "/tmp/pdb.pdb",
        800,
        id=pdb_id,
        vmd="/usr/local/bin/vmd",
        script_name="render-pdb.vmd",
        color="black",
        args=[str(start), str(stop), str(rotation), color_method],
        renderer=GPU_RENDERER,
        render_format=GPU_RENDER_FORMAT,
        frame_ext=GPU_FRAME_EXT,
        frame_dir="/tmp",
    )
    for rendered_frame in glob(f"/tmp/{pdb_id}*"):
        shutil.move(rendered_frame, f"/var/tmp/{Path(rendered_frame).name}")
    print(sorted(glob(f"/var/tmp/{pdb_id}*"))[:10])
    volume.commit()
    return True


@app.function(
    cpu=MOVIE_RENDER_CPU,
    image=image,
    volumes=volumes,
    timeout=1200,
)
def movie_render(name: str, pdb_id: str):
    volume.reload()
    movie_path = ternviz.movie(
        pdb_id,
        ffmpeg="ffmpeg",
        short_name=name,
        color="white",
    )
    with open(movie_path, "rb") as f:
        result = f.read()
    for entry in Path("/var/tmp").glob(f"{pdb_id}*"):
        if entry.is_dir():
            shutil.rmtree(entry)
        else:
            entry.unlink()
    volume.commit()
    return result


@app.function(
    secrets=[modal.Secret.from_name("web-auth-token")],
    image=image,
    timeout=1000,
)
@modal.asgi_app(custom_domains=["viz.proteincrow.ai"])
def modal_app():
    return build_web_app()
