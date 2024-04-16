from modal import Image, Volume, Stub, asgi_app, Secret
from fastapi import FastAPI, status, UploadFile, HTTPException, Depends
from fastapi.responses import Response
from fastapi.security import HTTPBearer, HTTPAuthorizationCredentials
import os
import uuid
import math
import gzip
import tarfile
from typing import List
import io
import warnings

VMD_URL = "https://storage.googleapis.com/fh-modal-artifacts/vmd-1.9.4a55.bin.LINUXAMD64-CUDA102-OptiX650-OSPRay185-RTXRTRT.opengl.tar.gz"
stub = Stub("ternviz-pdb")
vol = Volume.from_name("ternviz-scratch", create_if_missing=True)
volume_map = {"/var/tmp": vol}
web_app = FastAPI()
auth_scheme = HTTPBearer()

image = (
    Image.debian_slim(python_version="3.11")
    .apt_install("git", "wget", "libgomp1", "ffmpeg")
    .pip_install("git+https://github.com/whitead/ternviz.git", "mdanalysis")
    .run_commands(f"wget {VMD_URL} -O /tmp/vmd.tar.gz && tar -xvf /tmp/vmd.tar.gz -C /tmp")
    .run_commands("cd /tmp/vmd-* && ./configure && cd src && make install && hash")
)

with image.imports():
    import ternviz
    import MDAnalysis as mda
    from MDAnalysis.analysis import align

def get_backbone_selection(universe: mda.Universe, reference_universe: mda.Universe) -> str:
    largest_chain = max(universe.segments, key=lambda s: len(s.residues))
    return f"segid {largest_chain.segid} and backbone"

def align_pdbs_to_reference(pdb_files: List[bytes], reference_pdb: bytes) -> List[bytes]:
    ref = mda.Universe(io.StringIO(reference_pdb.decode('utf-8')), format='PDB')
    aligned_pdb_strings = []

    for pdb_data in pdb_files:
        u = mda.Universe(io.StringIO(pdb_data.decode('utf-8')), format='PDB')
        ref_selection = get_backbone_selection(u, ref)

        alignment = align.AlignTraj(u, ref, select=ref_selection, in_memory=True)
        alignment.run()

        buffer = io.StringIO()
        with mda.Writer(buffer, n_atoms=u.atoms.n_atoms, format='PDB') as writer:
            writer.write(u.atoms)
            aligned_pdb_strings.append(buffer.getvalue().encode('utf-8'))

    return aligned_pdb_strings


async def extract_gz_content(file_content: bytes) -> bytes:
    """Extract and return the content of a .gz file from bytes."""
    with gzip.open(io.BytesIO(file_content), 'rb') as gz:
        return gz.read()

async def extract_tar_content(file_content: bytes) -> List[tuple]:
    """Extract and return a list of tuples (filename, content) of .pdb files from a .tar or .tar.gz file in bytes."""
    pdb_files = []
    with tarfile.open(fileobj=io.BytesIO(file_content), mode='r:*') as tar:
        for member in tar.getmembers():
            if member.isfile() and member.name.endswith('.pdb'):
                extracted_file = tar.extractfile(member)
                if extracted_file:
                    pdb_files.append((member.name, extracted_file.read()))
    return pdb_files

async def process_upload_file(upload_file: UploadFile) -> List[tuple]:
    """Process an uploaded file and return a sorted list of tuples (filename, binary content) of .pdb files.
    
    Args:
        upload_file (UploadFile): The uploaded file object from FastAPI.
    
    Returns:
        List[tuple]: Sorted list of (filename, binary content) of .pdb files.
    """
    pdb_files = {}
    file_content = await upload_file.read()
    
    file_name = upload_file.filename
    if file_name.endswith('.pdb'):
        pdb_files[file_name] = file_content
    elif file_name.endswith('.gz') and not file_name.endswith('.tar.gz'):
        extracted_content = await extract_gz_content(file_content)
        if file_name[:-3].endswith('.pdb'):
            pdb_files[file_name[:-3]] = extracted_content
    elif file_name.endswith('.tar') or file_name.endswith('.tar.gz'):
        extracted_pdb_files = await extract_tar_content(file_content)
        for fname, content in extracted_pdb_files:
            pdb_files[fname] = content
    
    # Sort files by filename and return the tuples of filename and binary content
    sorted_pdb_files = [pdb_files[fname] for fname in sorted(pdb_files)]
    return sorted_pdb_files

@stub.function(cpu=4, image=image, volumes=volume_map, concurrency_limit=25)
def pdb_render(pdb_id: str, pdb_file: str, start: int, stop: int, rotation: float, color_method: str):
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
    )
    vol.commit()
    return True

@stub.function(cpu=8, image=image, volumes=volume_map, timeout=1200)
def movie(name: str, pdb_id: str):
    vol.reload()
    m = ternviz.movie(
            pdb_id,
            ffmpeg="ffmpeg",
            short_name=name,
            color="white",
    )
    # return the binary file
    with open(m, "rb") as f:
        result = f.read()
    for f in os.listdir("/var/tmp"):
        os.remove(f"/var/tmp/{f}")
    vol.commit()
    return result

@web_app.post("/render/")
async def render(pdb_file: UploadFile, 
                message: str = "",
                render_first_struct: bool = True,  
                target_total_frames: int = 360,
                token: HTTPAuthorizationCredentials = Depends(auth_scheme)):
    """Render a PDB file or a tarball of PDB files to a movie."""
    if token.credentials != os.environ["AUTH_TOKEN"]:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Incorrect bearer token",
            headers={"WWW-Authenticate": "Bearer"},
        )
    
    if not pdb_file.filename.endswith(('.pdb', '.gz', '.tar', '.tar.gz')):
        raise HTTPException(status_code=400, detail="Unsupported file type")

    pdb_files = await process_upload_file(pdb_file)
    if not pdb_files:
        raise HTTPException(status_code=400, detail="No PDB files found in the uploaded file")
    
    aligned_pdb_files = pdb_files
    if len(pdb_files) > 1:
        with warnings.catch_warnings():
            aligned_pdb_files = align_pdbs_to_reference(pdb_files, pdb_files[0])
        # if we are not rendering first
        if not render_first_struct and len(aligned_pdb_files) > 1:
            aligned_pdb_files = aligned_pdb_files[1:]

    coloring = "Structure" if len(pdb_files) == 1 else "Chain"

    # logic to schedule the render job

    target_total_frames = max(target_total_frames, len(pdb_files))
    max_frames_per_job = 30
    job_id = str(uuid.uuid4()).split("-")[0]

    # want to have enough structures so that when we split
    # we do not have too many structures on one node
    while target_total_frames // len(pdb_files) > max_frames_per_job:
        # increase by repeating each structure
        pdb_files = [pdb_files[i // 2] for i in range(len(pdb_files) * 2)]

    # make sure n_frames is exact multiple of len(pdb_files)
    # example:
    #   len(pdb_files) = 48
    #   target_total_frames = 360
    #   n_frames = 336 (closest multiple of 48 to 360)
    # example 2:
    #   len(pdb_files) = 450
    #   target_total_frames = 450
    #   n_frames = 450
    # example 3:
    #   len(pdb_files) = 1
    #   target_total_frames = 360
    #   while condition above will move up until
    #   len(pdb_files) = 2, 4, 8, 16, 32
    #   n_frames = 352
    n_frames = len(pdb_files) * math.ceil(target_total_frames / len(pdb_files))
    rotation = 360.0 / n_frames
    batch_size = n_frames // len(pdb_files)

    print("Will render", n_frames, "frames in", len(pdb_files), "batches of", batch_size, "frames each.")

    def _render_gen():
        index = 0
        for i, s in enumerate(pdb_files):
            yield job_id, s, index, index + batch_size, rotation, coloring
            index += batch_size
    
    # render map takes a generator
    _ = list(pdb_render.starmap(_render_gen()))

    blob = movie.remote(message, job_id)
    return Response(content=blob, media_type="media/mp4", headers={"Content-Disposition": f"attachment; filename=ternviz.mp4"})
    
@stub.function(secrets=[Secret.from_name("web-auth-token")], image=image, timeout=1000)
@asgi_app(custom_domains=["viz.proteincrow.ai"])
def app():
    return web_app