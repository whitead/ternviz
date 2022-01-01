import sys
from rdkit import Chem
from rdkit.Chem import AllChem
import tempfile
import os
import sys
from io import StringIO
Chem.WrapLogs()


def check_smiles(s):
    sio = sys.stderr = StringIO()
    Chem.MolFromSmiles(s)
    s = sio.getvalue()
    result = []
    if len(s) > 0:
        for si in s.split('\n'):
            result.append(si.split('SMILES Parse Error:')[-1])
    sys.stderr = sys.__stderr__
    return '\n'.join(result)


def vmd_script(width, id, high_quality=False):
    from importlib_resources import files
    import ternviz.vmd
    if high_quality:
        fp = files(ternviz.vmd).joinpath(
            'render.vmd')
    else:
        fp = files(ternviz.vmd).joinpath(
            'render-lq.vmd')
    result = []
    with fp.open('r') as f:
        for line in f.readlines():
            result.append(line.replace(
                'WIDTH', str(width)).replace('__ID__', id))
    return result


def gen_coords(s, name=None):
    m = Chem.MolFromSmiles(s)
    m = Chem.AddHs(m)
    AllChem.EmbedMolecule(m, randomSeed=0xf00d)
    AllChem.MMFFOptimizeMolecule(m)

    if name is None:
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix='.pdb')
    else:
        tmp = open(os.path.join('/var/tmp', name + '.pdb'), 'w')
    Chem.MolToPDBFile(m, tmp.name)
    return tmp


def render(pdb_path, width, id='movie', vmd='vmd', high_quality=True):
    with tempfile.NamedTemporaryFile() as script:
        with open(script.name, 'w') as f:
            f.writelines(vmd_script(width, id, high_quality))
        os.system(
            f'{vmd} -dispdev text -e {script.name} {pdb_path} > /dev/null')


def movie(name, ffmpeg='ffmpeg'):
    out = os.path.join('/var/tmp', f'{name}.mp4')
    os.system(
        f'{ffmpeg} -framerate 60 -f image2 -i /var/tmp/{name}.%04d.bmp -c:v h264 -crf 9 -c:v libx264 -movflags +faststart -filter_complex "[0:v]tpad=stop_mode=clone:stop_duration=2[a];[a]format=yuv420p[out]" -map "[out]" {out} > /dev/null')
    return out


def multiplex(videos, names):
    assert len(videos) == len(names)
    # $1 - i $2 - i $3 - filter_complex "[0]drawtext=text='${N1}':fontsize=36:x=(w-text_w)/2:y=(h-2*text_h):fontcolor=white[v0];[1]drawtext=text='${N2}':fontsize=36:x=(w-text_w)/2:y=(h-2*text_h):fontcolor=white[v1];[v0][v1]hstack=inputs=2[v]" - map "[v]" multi.mp4
