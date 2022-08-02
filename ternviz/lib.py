import sys
from rdkit import Chem
from rdkit.Chem import AllChem
import tempfile
import os
import sys
from io import StringIO, BytesIO
import requests
import subprocess
from rdkit.Chem import AllChem
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity
import urllib.request
from PIL import Image

Chem.WrapLogs()


def canonicalize(smiles):
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    smiles = Chem.MolToSmiles(mol, isomericSmiles=False, canonical=True)
    return smiles


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
    except:
        pass
    return Chem.rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(s))


def vmd_script(width, id, script_name="render.vmd", color="black"):
    from importlib_resources import files
    import ternviz.vmd

    fp = files(ternviz.vmd).joinpath(script_name)
    result = []
    with fp.open("r") as f:
        for line in f.readlines():
            result.append(
                line.replace("WIDTH", str(width))
                .replace("__ID__", id)
                .replace("BACKGROUND", color)
            )
    return result


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
    except:
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
        except Exception as e:
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
    except:
        pass
    try:
        AllChem.MMFFOptimizeMolecule(m)
    except:
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
):
    with tempfile.NamedTemporaryFile(delete=False) as script:
        with open(script.name, "w") as f:
            f.writelines(vmd_script(width, id, script_name, color=color))
        if len(args) > 0:
            arg_str = f'-args {" ".join(args)}'
        else:
            arg_str = ""
        subprocess.run(
            f"'{vmd}' -dispdev text  -eofexit '{pdb_path}' {arg_str} < '{script.name}'",
            shell=True,
        )


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
        print(pdbid)
        url = f"https://files.rcsb.org/download/{pdbid}.cif"
        pdb = requests.get(url)
        if name is None:
            tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".cif")
        else:
            tmp = open(os.path.join("/var/tmp", name + ".cif"), "wb")
        tmp.write(pdb.text.encode())
        print(tmp.name)
        return pdbid, tmp
    return None, None


def movie(name, short_name="molecule", color="white", ffmpeg="ffmpeg"):
    out = os.path.join("/var/tmp", f"{name}.mp4")
    font_path = os.path.join("/var/tmp", "CourierPrime-Regular.ttf")
    if not os.path.exists(font_path):
        urllib.request.urlretrieve(
            "https://github.com/google/fonts/raw/main/ofl/courierprime/CourierPrime-Regular.ttf",
            font_path,
        )
    subprocess.run(
        f"{ffmpeg} -framerate 30 -y -f image2 -i /var/tmp/{name}.%04d.bmp -c:v h264 -crf 9 "
        "-c:v libx264 -movflags +faststart -filter_complex "
        f"\"[0:v]drawtext=text='{short_name}':fontsize=48:x=(w-text_w)/2:y=(2*text_h):fontcolor={color}:fontfile={font_path}[c];"
        "[c]eq=saturation=0.8[d];"
        '[d]format=yuv420p[out]" '
        f'-map "[out]" {out} > /dev/null',
        shell=True,
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
