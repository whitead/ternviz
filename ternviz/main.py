import click
import os
from .lib import gen_coords, movie, render


@click.command()
@click.argument('smiles', nargs=-1)
@click.option('--names', default=None, help='comma separated names')
@click.option('--vmd', default='vmd')
@click.option('--ffmpeg', default='ffmpeg')
def main(smiles, names, vmd, ffmpeg):
    if len(smiles) == 0:
        return
    if type(smiles) == str:
        smiles = [smiles]
    if names is None:
        names = [f'Molecule {i+1}' for i in range(len(smiles))]
    if type(names) == str:
        names = names.split(',')
    assert len(names) == len(smiles)
    width = 800 // len(smiles)
    for i, (n, s) in enumerate(zip(names, smiles)):
        print('Processing SMILES', i+1, 'of', len(smiles))
        p = gen_coords(s)
        print('Rendering', s)
        render(p.name, width, vmd=vmd)
        print('Making Movie for', s)
        m = movie(n, ffmpeg=ffmpeg)
        print(m)
        p.close()
        os.unlink(p.name)
