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


def vmd_script(width, id, high_quality=False):
    from importlib_resources import files
    import ternviz.vmd

    if high_quality:
        fp = files(ternviz.vmd).joinpath("render.vmd")
    else:
        fp = files(ternviz.vmd).joinpath("render-lq.vmd")
    result = []
    with fp.open("r") as f:
        for line in f.readlines():
            result.append(line.replace("WIDTH", str(width)).replace("__ID__", id))
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
            print("Trying 3D molecule", cid)
            yield m
        else:
            print("Trying 2D molecule", cid)
            m3 = Chem.AddHs(m)
            if AllChem.EmbedMolecule(m3) == 0:
                print("Embedded it")
                yield m3
            else:
                print("Passed as 2D", cid)
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


def render(pdb_path, width, id="movie", vmd="vmd", high_quality=True):
    with tempfile.NamedTemporaryFile() as script:
        with open(script.name, "w") as f:
            f.writelines(vmd_script(width, id, high_quality))
        subprocess.run(
            f"{vmd} -dispdev text -eofexit {pdb_path} < {script.name} > /dev/null",
            shell=True,
        )


def movie(name, short_name="molecule", ffmpeg="ffmpeg"):
    out = os.path.join("/var/tmp", f"{name}.mp4")
    font_path = os.path.join(
        os.getenv("CONDA_PREFIX"), "fonts", "open-fonts", "IBMPlexMono-Light.ttf"
    )
    # os.system(
    #     f'{ffmpeg} -framerate 60 -f image2 -i /var/tmp/{name}.%04d.bmp -c:v h264 -crf 9 '
    #     '-c:v libx264 -movflags +faststart -filter_complex '
    #     '"[0:v]tpad=stop_mode=clone:stop_duration=2[b];'
    #     f'[b]drawtext=text=\'{short_name}\':fontsize=36:x=(w-text_w)/2:y=(2*text_h):fontcolor=white:fontfile={font_path}[c];'
    #     '[c]format=yuv420p[out]" '
    #     f'-map "[out]" {out} > /dev/null')
    subprocess.run(
        f"{ffmpeg} -framerate 60 -f image2 -i /var/tmp/{name}.%04d.bmp -c:v h264 -crf 9 "
        "-c:v libx264 -movflags +faststart -filter_complex "
        f"\"[0:v]drawtext=text='{short_name}':fontsize=36:x=(w-text_w)/2:y=(2*text_h):fontcolor=white:fontfile={font_path}[c];"
        '[c]format=yuv420p[out]" '
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
