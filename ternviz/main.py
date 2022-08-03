from pydoc import cli
import click
import os
from .lib import (
    bmp2png,
    gen_coords,
    movie,
    render,
    check_smiles,
    get_name,
    multiplex,
    get_pdb,
    align,
    bmp2png,
)


@click.command()
@click.argument("smiles", nargs=-1)
@click.option("--names", default=None, help="comma separated names")
@click.option("--vmd", default="vmd")
@click.option("--ffmpeg", default="ffmpeg")
@click.option("--color", default="black")
@click.option("--low-quality", is_flag=True, default=False)
def main(smiles, names, vmd, ffmpeg, low_quality, color):
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
            color=color,
            script_name="render-lq.vmd" if low_quality else "render.vmd",
        )
        print("Making Movie for", s)
        m = movie(
            n,
            ffmpeg=ffmpeg,
            short_name=n,
            color="black" if color == "white" else "white",
        )
        movies.append(m)
        p.close()
        os.unlink(p.name)
    if len(movies) == 2:
        return multiplex(movies, "out")

    os.rename(movies[0], f"{n}.mp4")
    return movies[0]


@click.command()
@click.argument("pdb-query", nargs=-1)
@click.option("--vmd", default="vmd")
@click.option("--color", default="black")
@click.option("--scolor", default="Structure")
@click.option("--name", default=None)
@click.option("--ffmpeg", default="ffmpeg")
@click.option("--name", default=None)
@click.option("--frames", default=None)
def pdb_main(pdb_query, vmd, color, ffmpeg, name, scolor, frames=None):
    multi = False
    rot_i = 1.0
    if len(pdb_query) == 1:
        pdb_query = pdb_query[0]
        # check if it's an actual file
        if not os.path.isfile(pdb_query):
            pdb_id, p = get_pdb(pdb_query)
        else:
            pdb_id = "pdb"
            p = open(pdb_query, "r")
        path = p.name
        if frames is None:
            frames = 1
        else:
            rot_i = 360.0 / int(frames)
    else:
        pdb_id = str(1)
        multi = True
        path = pdb_query
        if frames is None:
            frames = 600
    if not pdb_id:
        raise ValueError("Failed to find pdb")

    print("Rendering", pdb_id, "in file", path, "multi = ", multi)
    if multi:
        N = len(pdb_query)
        for i in range(N):
            render(
                os.path.abspath(path[i]),
                800,
                id=pdb_id,
                vmd=vmd,
                script_name="render-pdb.vmd",
                color=color,
                args=[
                    str(i * frames // N),
                    str((i + 1) * frames // N),
                    str(rot_i),
                    scolor,
                ],
            )
    else:
        render(
            path,
            800,
            id=pdb_id,
            vmd=vmd,
            script_name="render-pdb.vmd",
            color=color,
            args=["0", str(frames), str(rot_i), scolor],
        )
    if not multi:
        p.close()
    if frames == 1:
        m = bmp2png(f"/var/tmp/{pdb_id}.0000.bmp", f"{pdb_id}.png")
    else:
        print("Making Movie for", pdb_id)
        m = movie(
            pdb_id,
            ffmpeg=ffmpeg,
            short_name=name if name else pdb_id,
            color="black" if color == "white" else "white",
        )
        m = os.rename(m, f"{name if name else pdb_id}.mp4")
    return m


@click.command()
@click.argument("ref")
@click.option("--sel", default="protein", help="Atom selection")
@click.argument("aligns", nargs=-1)
def align_cmd(ref, sel, aligns):
    align(ref, sel, *aligns)
