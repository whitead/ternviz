import click
import os
from .lib import (
    gen_coords,
    movie,
    render,
    check_smiles,
    get_name,
    multiplex,
    get_pdb,
    align,
)


@click.command()
@click.argument("smiles", nargs=-1)
@click.option("--names", default=None, help="comma separated names")
@click.option("--vmd", default="vmd")
@click.option("--ffmpeg", default="ffmpeg")
@click.option("--low-quality", is_flag=True, default=False)
def main(smiles, names, vmd, ffmpeg, low_quality):
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
        render(
            p.name,
            width,
            id=n,
            vmd=vmd,
            script_name="render-lq.vmd" if low_quality else "render.vmd",
        )
        print("Making Movie for", s)
        m = movie(n, ffmpeg=ffmpeg, short_name=n)
        movies.append(m)
        p.close()
        os.unlink(p.name)
    if len(movies) == 2:
        return multiplex(movies, "out")

    return movies[0]


@click.command()
@click.argument("pdb-query")
@click.option("--vmd", default="vmd")
@click.option("--color", default="black")
@click.option("--ffmpeg", default="ffmpeg")
@click.option("--width", default=800)
@click.option("--name", default=None)
def pdb_main(pdb_query, vmd, color, ffmpeg, width, name):
    if ".pdb" in pdb_query and os.path.exists(pdb_query):
        print("Assuming you meant me to use path")
        pdb_id = pdb_query.split(".pdb")[0]
        p = open(pdb_query, "r")
    else:
        pdb_id, p = get_pdb(pdb_query)
    if not pdb_id:
        raise ValueError("Failed to find pdb")
    print("Rendering", pdb_id, "in file", p.name)
    render(p.name, width, id=pdb_id, vmd=vmd, script_name="render.vmd", color=color)
    print("Making Movie for", pdb_id)
    p.close()
    # os.unlink(p.name)
    m = movie(
        pdb_id,
        ffmpeg=ffmpeg,
        short_name=pdb_id if name is None else name,
        color="black" if color == "white" else "white",
    )
    return m


@click.command()
@click.argument("ref")
@click.argument("aligns", nargs=-1)
def align_cmd(ref, aligns):
    align(ref, *aligns)
