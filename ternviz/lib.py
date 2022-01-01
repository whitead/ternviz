import sys
from rdkit import Chem
from rdkit.Chem import AllChem
import tempfile
import os


def vmd_script(width):
    from importlib_resources import files
    import ternviz.vmd
    fp = files(ternviz.vmd).joinpath(
        'render.vmd')
    result = []
    with fp.open('r') as f:
        for line in f.readlines():
            result.append(line.replace('WIDTH', str(width)))
    return result


def gen_coords(s):
    m = Chem.MolFromSmiles(s)
    m = Chem.AddHs(m)
    AllChem.EmbedMolecule(m, randomSeed=0xf00d)
    AllChem.MMFFOptimizeMolecule(m)

    tmp = tempfile.NamedTemporaryFile(delete=False, suffix='.pdb')
    Chem.MolToPDBFile(m, tmp.name)
    return tmp


def render(pdb_path, width, vmd='vmd'):
    with tempfile.NamedTemporaryFile() as script:
        with open(script.name, 'w') as f:
            f.writelines(vmd_script(width))
        os.system(
            f'{vmd} -dispdev text -e {script.name} {pdb_path} > /dev/null')


def movie(name, ffmpeg='ffmpeg'):
    os.system(f'{ffmpeg} -framerate 60 -f image2 -i /var/tmp/mov.%04d.bmp -c:v h264 -crf 9 -c:v libx264 -movflags +faststart -vf format=yuv420p {name}.mp4')
    return f'{name}.mp4'


def multiplex(videos, names):
    assert len(videos) == len(names)
    # $1 - i $2 - i $3 - filter_complex "[0]drawtext=text='${N1}':fontsize=36:x=(w-text_w)/2:y=(h-2*text_h):fontcolor=white[v0];[1]drawtext=text='${N2}':fontsize=36:x=(w-text_w)/2:y=(h-2*text_h):fontcolor=white[v1];[v0][v1]hstack=inputs=2[v]" - map "[v]" multi.mp4
