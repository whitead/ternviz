import os
import re
import shutil
import subprocess
import sys
import tempfile
import urllib.request
from importlib.resources import files
from io import BytesIO, StringIO
from pathlib import Path

import requests
from PIL import Image
from rdkit import Chem
from rdkit.Chem import AllChem

FRAME_EXTENSIONS = ("bmp", "png", "ppm", "tga")


def canonicalize(smiles):
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    smiles = Chem.MolToSmiles(mol, isomericSmiles=False, canonical=True)
    return smiles


def sdf2pdb(sdf_path):
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb")
    m = Chem.MolFromMolFile(sdf_path)
    Chem.MolToPDBFile(m, tmp.name)
    return tmp


def check_smiles(s):
    sio = sys.stderr = StringIO()
    Chem.MolFromSmiles(s)
    s = sio.getvalue()
    result = []
    if len(s) > 0:
        for si in s.split("\n"):
            result.append(si.split("SMILES Parse Error:")[-1])
    sys.stderr = sys.__stderr__
    return "\n".join(result)


def get_name(s):
    try:
        url = "https://api.leruli.com/v21_4/graph-to-name"
        reply = requests.post(url, json={"graph": s})
        data = reply.json()
        if data["reference"] == "wikidata":
            return data["name"].replace(" ", "-")
    except Exception:
        pass
    return Chem.rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(s))


def vmd_script(
    width,
    id,
    script_name="render.vmd",
    color="black",
    strip_display_resize=False,
    renderer="TachyonInternal",
    render_format="BMP",
    frame_ext="bmp",
    frame_dir="/var/tmp",
):
    fp = files("ternviz.vmd").joinpath(script_name)
    result = []
    with fp.open("r") as f:
        for line in f.readlines():
            if strip_display_resize and line.lstrip().startswith("display resize "):
                continue
            result.append(
                line.replace("WIDTH", str(width))
                .replace("__ID__", id)
                .replace("BACKGROUND", color)
                .replace("RENDERER", renderer)
                .replace("RENDER_FORMAT", render_format)
                .replace("FRAME_EXT", frame_ext)
                .replace("FRAME_DIR", frame_dir)
            )
    return result


def _find_vmd_bundle_binary(app_path):
    contents_dir = app_path / "Contents"
    binaries = sorted(contents_dir.glob("vmd*/vmd_*"))
    return binaries[0] if binaries else None


def _macos_vmd_env(vmd_path):
    vmd_path = Path(vmd_path)
    parts = vmd_path.parts
    if ".app" not in vmd_path.as_posix() or "Contents" not in parts:
        return {}
    if not vmd_path.name.startswith("vmd_"):
        return {}

    vmddir = vmd_path.parent
    arch = vmd_path.name.split("vmd_", 1)[1]
    python_dir = vmddir / "scripts" / "python"
    existing_pythonpath = os.environ.get("PYTHONPATH")
    env = {
        "VMDDIR": str(vmddir),
        "TCL_LIBRARY": str(vmddir / "scripts" / "tcl"),
        "PYTHONPATH": (
            f"{python_dir}:{existing_pythonpath}"
            if existing_pythonpath
            else str(python_dir)
        ),
    }
    for env_name, base_name in (
        ("STRIDE_BIN", "stride"),
        ("SURF_BIN", "surf"),
        ("TACHYON_BIN", "tachyon"),
    ):
        candidate = vmddir / f"{base_name}_{arch}"
        if candidate.exists():
            env[env_name] = str(candidate)
    return env


def resolve_vmd(vmd):
    candidate = Path(vmd).expanduser()
    if candidate.is_dir() and candidate.suffix == ".app":
        bundle_binary = _find_vmd_bundle_binary(candidate)
        if bundle_binary is not None:
            return str(bundle_binary), _macos_vmd_env(bundle_binary)
    if candidate.exists():
        return str(candidate), _macos_vmd_env(candidate)
    if shutil.which(vmd):
        return vmd, {}

    if sys.platform == "darwin" and vmd == "vmd":
        for pattern in ("VMD2*.app", "VMD*.app"):
            for app_path in sorted(Path("/Applications").glob(pattern), reverse=True):
                bundle_binary = _find_vmd_bundle_binary(app_path)
                if bundle_binary is not None:
                    return str(bundle_binary), _macos_vmd_env(bundle_binary)

    return vmd, {}


def find_template(smiles, count=10):
    # clean-up smiles
    smiles = requests.utils.quote(smiles)
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/smiles/{smiles}/cids"
    try:
        reply = requests.get(
            url,
            params={"Threshold": 95, "MaxRecords": count},
            headers={"accept": "application/json"},
            timeout=10,
        )
    except requests.exceptions.Timeout:
        print("Pubchem seems to be down right now")
        return None
    try:
        data = reply.json()
    except Exception:
        print("Could not find a match")
        return None
    if "IdentifierList" not in data:
        print("Could not find a match")
        return None
    cids = data["IdentifierList"]["CID"]

    # download and get tanimoto (3D first)
    def get_record(cid):
        def try_get(cid, require_3d):
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}"
            try:
                reply = requests.get(
                    url,
                    params={"record_type": "3d" if require_3d else "2d"},
                    headers={"accept": "chemical/x-mdl-sdfile"},
                    timeout=10,
                )
            except requests.exceptions.Timeout:
                print("Pubchem seems to be down right now")
                return None
            if reply.status_code != 200:
                return None
            data = reply.text
            with Chem.ForwardSDMolSupplier(BytesIO(data.encode())) as fsuppl:
                for mol in fsuppl:
                    if mol:
                        return mol
            return None

        m = try_get(cid, True)
        if m:
            return m, True
        return try_get(cid, False), False

    for cid in cids:
        m, is_3d = get_record(cid)
        if is_3d:
            yield m
        else:
            m3 = Chem.AddHs(m)
            if AllChem.EmbedMolecule(m3) == 0:
                yield m3
            else:
                yield m


def gen_coords(s, name=None, template=None):
    s = canonicalize(s)
    m = Chem.MolFromSmiles(s)
    m = Chem.AddHs(m)
    status = -1
    if not template:
        # try pubchem templates (for a while)
        for t in find_template(s, count=10):
            try:
                if t:
                    m = AllChem.ConstrainedEmbed(m, t)
                    status = 0
                    break
            except Exception as e:
                print("Failed on template")
                print(e)
                continue
    else:
        try:
            m = AllChem.ConstrainedEmbed(m, template)
            status = 0
        except Exception:
            print("Trying ignore template")
    if status == -1:
        print("Failed to match template, embedding alone")
        status = AllChem.EmbedMolecule(m)
    if status == -1:
        # try random coords
        print("Trying random coordinates and fancy coord gen")
        params = AllChem.ETKDGv3()
        params.useSmallRingTorsions = True
        params.useRandomCoords = True
        status = AllChem.EmbedMolecule(m, params)
    try:
        AllChem.UFFOptimizeMolecule(m)
    except Exception:
        pass
    try:
        AllChem.MMFFOptimizeMolecule(m)
    except Exception:
        pass

    if name is None:
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb")
    else:
        tmp = open(os.path.join("/var/tmp", name + ".pdb"), "w")
    Chem.MolToPDBFile(m, tmp.name)
    return tmp


def render(
    pdb_path,
    width,
    id="movie",
    script_name="render.vmd",
    color="black",
    vmd="vmd",
    args=[],
    renderer="TachyonInternal",
    render_format="BMP",
    frame_ext="bmp",
    height=None,
    frame_dir="/var/tmp",
):
    vmd_executable, vmd_env = resolve_vmd(vmd)
    # Batch renders always use `-dispdev text -size ...`, so `display resize`
    # is redundant and can terminate VMD 2 scripts before rendering.
    strip_display_resize = True
    script = tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".tcl")
    try:
        script.writelines(
            vmd_script(
                width,
                id,
                script_name,
                color=color,
                strip_display_resize=strip_display_resize,
                renderer=renderer,
                render_format=render_format,
                frame_ext=frame_ext,
                frame_dir=frame_dir,
            )
        )
        script.close()
        if height is None:
            height = max(1, int(width * 0.75))
        command = [
            vmd_executable,
            "-dispdev",
            "text",
            "-size",
            str(width),
            str(height),
            "-e",
            script.name,
            "-eofexit",
            pdb_path,
        ]
        if args:
            command.extend(["-args", *args])
        env = os.environ.copy()
        env.update(vmd_env)
        subprocess.run(command, check=True, env=env)
    finally:
        os.unlink(script.name)


def _frame_files(name, frame_dir="/var/tmp"):
    frame_dir = Path(frame_dir)
    pattern = re.compile(rf"^{re.escape(name)}\.\d{{4}}\.(?P<ext>[A-Za-z0-9]+)$")
    frames = []
    for candidate in sorted(frame_dir.iterdir()):
        match = pattern.match(candidate.name)
        if match and match.group("ext").lower() in FRAME_EXTENSIONS:
            frames.append(candidate)
    return frames


def first_frame(name, frame_dir="/var/tmp"):
    frames = _frame_files(name, frame_dir=frame_dir)
    if not frames:
        raise FileNotFoundError(f"No rendered frames found for {name!r} in {frame_dir}")
    return frames[0]


def frame_pattern(name, frame_dir="/var/tmp"):
    first = first_frame(name, frame_dir=frame_dir)
    suffix = first.suffix
    return str(Path(frame_dir) / f"{name}.%04d{suffix}")


def _resolve_movie_font_path(frame_dir="/var/tmp"):
    for candidate in (
        Path("/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf"),
        Path("/System/Library/Fonts/Menlo.ttc"),
        Path("/System/Library/Fonts/Supplemental/Courier New.ttf"),
    ):
        if candidate.exists():
            return str(candidate)

    font_path = Path(frame_dir) / "CourierPrime-Regular.ttf"
    if not font_path.exists():
        urllib.request.urlretrieve(
            "https://github.com/google/fonts/raw/main/ofl/courierprime/CourierPrime-Regular.ttf",
            font_path,
        )
    return str(font_path)


def get_pdb(query_string, name=None):
    url = "https://search.rcsb.org/rcsbsearch/v2/query?json={search-request}"
    query = {
        "query": {
            "type": "terminal",
            "service": "full_text",
            "parameters": {"value": query_string},
        },
        "return_type": "entry",
    }
    r = requests.post(url, json=query)
    if "result_set" in r.json() and len(r.json()["result_set"]) > 0:
        pdbid = r.json()["result_set"][0]["identifier"]
        url = f"https://files.rcsb.org/download/{pdbid}.cif"
        pdb = requests.get(url)
        if name is None:
            tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".cif")
        else:
            tmp = open(os.path.join("/var/tmp", name + ".cif"), "wb")
        tmp.write(pdb.text.encode())
        return pdbid, tmp
    return None, None


def movie(name, short_name="molecule", color="white", ffmpeg="ffmpeg"):
    out = os.path.join("/var/tmp", f"{name}.mp4")
    frame_input = frame_pattern(name)
    font_path = _resolve_movie_font_path()
    filter_complex = (
        f"[0:v]drawtext=text='{short_name}':fontsize=48:x=(w-text_w)/2:"
        f"y=(2*text_h):fontcolor={color}:fontfile={font_path}[c];"
        "[c]eq=saturation=0.8[d];"
        "[d]format=yuv420p[out]"
    )
    subprocess.run(
        [
            ffmpeg,
            "-framerate",
            "30",
            "-y",
            "-f",
            "image2",
            "-i",
            frame_input,
            "-c:v",
            "libx264",
            "-preset",
            "veryfast",
            "-crf",
            "20",
            "-movflags",
            "+faststart",
            "-filter_complex",
            filter_complex,
            "-map",
            "[out]",
            out,
        ],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return out


def multiplex(videos, name, ffmpeg="ffmpeg"):
    assert len(videos) == 2
    out = os.path.join("/var/tmp", f"{name}.mp4")
    os.system(
        f"{ffmpeg} -i {videos[0]} -i {videos[1]} -filter_complex hstack=inputs=2 {out} > /dev/null"
    )
    return out


def concat(videos, name, ffmpeg="ffmpeg"):
    assert len(videos) == 2
    out = os.path.join("/var/tmp", f"{name}.mp4")
    os.system(
        f'{ffmpeg} -i {videos[0]} -i {videos[1]} -y -filter_complex "[0:v] [1:v]  \
concat" {out} > /dev/null'
    )
    return out


def bmp2png(bmp_path, png_path):
    im = Image.open(bmp_path)
    im.save(png_path, "png")


def align(ref, sel, *args):
    import MDAnalysis
    from MDAnalysis.analysis import align

    ref = MDAnalysis.Universe(ref)

    for f in args:
        if os.path.exists(f):
            out_f = f.split(".pdb")[0] + "-align.pdb"
            p = MDAnalysis.Universe(f)
            alignment = align.AlignTraj(
                p,
                ref,
                select=sel,
                # weights=p.select_atoms("backbone").bfactors / 100,
                filename=out_f,
            )
            alignment.run()
        else:
            raise FileNotFoundError("Could not find structure", f)
