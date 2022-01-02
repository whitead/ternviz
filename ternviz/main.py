import click
import os
from .lib import gen_coords, movie, render, check_smiles, get_name, multiplex


@click.command()
@click.argument("smiles", nargs=-1)
@click.option("--names", default=None, help="comma separated names")
@click.option("--vmd", default="vmd")
@click.option("--ffmpeg", default="ffmpeg")
def main(smiles, names, vmd, ffmpeg):
    if len(smiles) == 0:
        return
    if type(smiles) == str:
        smiles = [smiles]
    if names is None:
        names = [get_name(s) for s in smiles]
    if type(names) == str:
        names = names.split(",")
    assert len(names) == len(smiles)
    width = 800 // len(smiles)
    movies = []
    for i, (n, s) in enumerate(zip(names, smiles)):
        print("Processing SMILES", i + 1, "of", len(smiles))
        status = ""
        sm = check_smiles(s)
        if sm:
            raise ValueError(sm)
        p = gen_coords(s)
        if not p:
            raise ValueError("Failed to generate coordinates")
        print("Rendering", s)
        render(p.name, width, id=n, vmd=vmd)
        print("Making Movie for", s)
        m = movie(n, ffmpeg=ffmpeg, short_name=n)
        movies.append(m)
        p.close()
        os.unlink(p.name)
    if len(movies) == 2:
        return multiplex(movies, "out")

    return movies[0]
