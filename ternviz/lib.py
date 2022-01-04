import sys
from rdkit import Chem
from rdkit.Chem import AllChem
import tempfile
import os
import sys
from io import StringIO, BytesIO
import requests
import subprocess

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


def find_template(smiles):
    # clean-up smiles
    smiles = requests.utils.quote(smiles)
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/smiles/{smiles}/cids"
    try:
        reply = requests.get(
            url,
            params={"MaxRecords": 1},
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
    if len(data["IdentifierList"]["CID"]) == 0:
        print("Could not find a match")
        return None
    cid = data["IdentifierList"]["CID"][0]
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}"
    try:
        reply = requests.get(
            url,
            headers={"accept": "chemical/x-mdl-sdfile"},
            timeout=10,
        )
    except requests.exceptions.Timeout:
        print("Pubchem seems to be down right now")
        return None
    data = reply.text
    with Chem.ForwardSDMolSupplier(BytesIO(data.encode())) as fsuppl:
        for mol in fsuppl:
            if mol:
                return mol
    return None


def gen_coords(s, name=None, template=None):
    s = canonicalize(s)
    m = Chem.MolFromSmiles(s)
    m = Chem.AddHs(m)
    if template:
        m = AllChem.ConstrainedEmbed(m, template)
        status = 0
    else:
        status = AllChem.EmbedMolecule(m)
    if status == -1:
        # try random coords
        print("Trying random coordinates")
        ps = AllChem.ETKDGv2()
        ps.useRandomCoords = True
        status = AllChem.EmbedMolecule(m, ps)
    if status == -1:
        # try to find template molecule
        print("Trying to find template molecule instead")
        return gen_coords(s, name=name, template=find_template(s))
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
