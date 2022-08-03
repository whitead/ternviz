import tweepy
import click
from dotenv import load_dotenv
from rdkit import Chem
import os
import ternviz
import multiprocessing
import time

TIMEOUT_COORDS = 500
TIMEOUT_RENDER = 5000


def coord_error(api, id, s):
    return api.update_status(
        status=f"Failed to generate coordinates for your query 😢\n Your query was '{s[:80]}'",
        in_reply_to_status_id=id,
        auto_populate_reply_metadata=True,
    )


def render_error(api, id):
    api.update_status(
        status=f"Was able to make coordinates, but rendering took too long 😢 It might be that @andrewwhite01 is broke and cannot afford enough server time🤏🏻",
        in_reply_to_status_id=id,
        auto_populate_reply_metadata=True,
    )


def full_text(api, status_id):
    status = api.get_status(status_id, tweet_mode="extended")
    try:
        text = status.retweeted_status.full_text
    except AttributeError:  # Not a Retweet
        text = status.full_text
    return text


def get_pdb_name(pdb_file):
    with open(pdb_file, "r") as f:
        return f.readline().split()[-1].split("_")[-1]


def make_pdb(api, id, query_text, template=None, protein=False):
    extension = "pdb"
    if protein:
        extension = "cif"
        print("Finding protein query", query_text)
        p = multiprocessing.Process(target=ternviz.get_pdb, args=(query_text, str(id)))
    else:
        smiles = query_text
        smiles_status = ternviz.check_smiles(smiles)
        if smiles_status:
            se = coord_error(api, id, smiles)
            api.update_status(
                status=smiles_status[:170],
                in_reply_to_status_id=se.id,
                auto_populate_reply_metadata=True,
            )
            return None, None
        if template:
            print("Loading template")
            template_mol = Chem.MolFromPDBFile(template)
        else:
            template_mol = None
        print("Determined SMILES are", smiles)
        p = multiprocessing.Process(
            target=ternviz.gen_coords, args=(smiles, str(id), template_mol)
        )
    p.start()
    p.join(TIMEOUT_COORDS)
    if p.is_alive():
        p.terminate()
        p.join()
        coord_error(api, id, query_text)
        print("Failed on timeout")
        return None, None, None
    out = os.path.join("/var/tmp", f"{id}.{extension}")
    if not os.path.exists(out):
        print("Could not find file - ", out)
        coord_error(api, id, query_text)
        return None, None
    return out, get_pdb_name(out) if protein else ternviz.get_name(smiles)


def do_render(api, pdb_file, id, short_name, width, protein, color):
    p = multiprocessing.Process(
        target=ternviz.render,
        args=(
            pdb_file,
            width,
            str(id),
            "render-pdb.vmd" if protein else "render.vmd",
            color,
        ),
    )
    p.start()
    p.join(TIMEOUT_RENDER)
    if p.is_alive():
        p.terminate()
        p.join()
        render_error(api, id)
        return None
    p2 = multiprocessing.Process(
        target=ternviz.movie,
        args=(str(id), short_name, "white" if color == "black" else "black"),
    )
    p2.start()
    p2.join(TIMEOUT_RENDER)
    out = os.path.join("/var/tmp", f"{id}.mp4")
    if not os.path.exists(out):
        render_error(api, id)
        return None
    return out


def render(
    api, status, smiles, width=800, template_pdb=None, protein=False, color="black"
):
    pdb_file, short_name = make_pdb(
        api, status.id, smiles, template=template_pdb, protein=protein
    )
    if not pdb_file:
        return None, None, None
    try:
        if protein:
            msg = f"{short_name} is sure a pretty protein 😍 Working on rendering 🎨... This can take a few minutes. I will delete this message when I'm done."
        else:
            msg = f"These SMILES look great 😁 Working on rendering {short_name} 🎨... This can take a few minutes. I will delete this message when I'm done."
        init_status = api.update_status(
            status=msg,
            in_reply_to_status_id=status.id,
            auto_populate_reply_metadata=True,
        )
    except tweepy.os.error.HTTPException as e:
        return None, None
    movie = do_render(
        api,
        pdb_file,
        str(status.id) + short_name,
        short_name,
        width=width,
        protein=protein,
        color=color,
    )
    api.destroy_status(init_status.id)
    return movie, short_name, pdb_file


@click.command()
@click.argument("user", type=str, default="ternviz")
def bot(user):
    load_dotenv(".secrets")
    consumer_key = os.environ.get("CONSUMER_KEY")
    consumer_secret = os.environ.get("CONSUMER_SECRET")
    access_token = os.environ.get("ACCESS_TOKEN")
    access_token_secret = os.environ.get("ACCESS_TOKEN_SECRET")
    auth = tweepy.OAuthHandler(consumer_key, consumer_secret)
    auth.set_access_token(access_token, access_token_secret)

    api = tweepy.API(auth)
    my_user_id = api.get_user(screen_name="ternviz").id

    # Subclass Stream to print IDs of Tweets received
    class IDPrinter(tweepy.Stream):
        def on_status(self, status):
            print("Processing", status.id)
            call_type = ""
            protein = False
            color = "white"
            text = full_text(api, status.id).replace("please", "")
            if "render" in text:
                call_type = "render"
            elif "compare" in text:
                call_type = "compare"
            if "protein" in text:
                protein = True
                text = text.replace("protein", "")
            if "white" in text:
                color = "white"
                text = text.replace("white", "")
            if "black" in text:
                color = "black"
                text = text.replace("black", "")
            if not call_type:
                print("Determined to be reply with no render or compare in body")
                return
            if status.user.id == my_user_id:
                print("Determined to be self-reply")
                return
            if text.startswith("rt @") == True:
                print("Determined to be RT")
                return
            start = time.time_ns()
            if call_type == "render":
                query_text = text.split("render")[-1].strip()
                movie, short_name, _ = render(
                    api, status, query_text, protein=protein, color=color
                )
            elif call_type == "compare":
                query_text = list(text.split("compare")[-1].split())
                m1, s1, pdb_file = render(
                    api,
                    status,
                    query_text[0],
                    width=800 // 2,
                    protein=protein,
                    color=color,
                )
                m2, s2, _ = render(
                    api,
                    status,
                    query_text[1],
                    width=800 // 2,
                    template_pdb=pdb_file,
                    protein=protein,
                    color=color,
                )
                movie = ternviz.multiplex([m1, m2], status.id)
                short_name = s1 + " vs. " + s2
            if not movie:
                return
            print("Completed Rendering")
            media = api.media_upload(movie, media_category="tweet_video")
            duration = (time.time_ns() - start) / (10**9)
            try:
                api.update_status(
                    status=f"Rendered {short_name} in {duration:.2f} seconds.",
                    media_ids=[media.media_id],
                    in_reply_to_status_id=status.id,
                    auto_populate_reply_metadata=True,
                )
            except tweepy.os.error.HTTPException as e:
                return

    # Initialize instance of the subclass
    printer = IDPrinter(
        consumer_key, consumer_secret, access_token, access_token_secret
    )

    # Filter realtime Tweets by keyword
    print("Set-up done, starting bot")

    while True:
        try:
            p = printer.filter(track=["@ternviz"], threaded=True)
            p.join()
        except Exception as e:
            time.sleep(1.0)
            print("Died from", e)
            print("restarting")
            continue
