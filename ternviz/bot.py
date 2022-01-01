import tweepy
import click
from dotenv import load_dotenv
from rdkit import Chem
import os
import ternviz
import multiprocessing
import time
from retrying import retry

TIMEOUT_COORDS = 10
TIMEOUT_RENDER = 1000


def smiles_error(api, id, s):
    return api.update_status(status=f"Failed to generate coordinates for your SMILES 😢\n Your SMILES were '{s}'",
                             in_reply_to_status_id=id, auto_populate_reply_metadata=True)


def render_error(api, id):
    api.update_status(status=f"Was able to make coordinates, but rendering took too long 😢 It might be that @andrewwhite01 is broke and cannot afford enough server time🤏🏻",
                      in_reply_to_status_id=id, auto_populate_reply_metadata=True)


def make_pdb(api, id, text):
    smiles = text.split('render')[-1].strip()
    print('Tweet text', text)
    smiles_status = ternviz.check_smiles(smiles)
    if smiles_status:
        se = smiles_error(api, id, smiles)
        api.update_status(status=smiles_status,
                          in_reply_to_status_id=se.id, auto_populate_reply_metadata=True)
        return None

    print('Determined SMILES are', smiles)
    p = multiprocessing.Process(
        target=ternviz.gen_coords, args=(smiles, str(id)))
    p.start()
    p.join(TIMEOUT_COORDS)
    if p.is_alive():
        p.terminate()
        p.join()
        smiles_error(api, id, smiles)
        return None
    out = os.path.join('/var/tmp', f'{id}.pdb')
    if not os.path.exists(out):
        smiles_error(api, id, smiles)
        return None
    return out


def render(api, pdb_file, id, width=800):
    p = multiprocessing.Process(
        target=ternviz.render, args=(pdb_file, width, str(id)))
    p.start()
    p.join(TIMEOUT_RENDER)
    if p.is_alive():
        p.terminate()
        p.join()
        render_error(api, id)
        return None
    p2 = multiprocessing.Process(
        target=ternviz.movie, args=(str(id),))
    p2.start()
    p2.join(TIMEOUT_RENDER)
    out = os.path.join('/var/tmp', f'{id}.mp4')
    if not os.path.exists(out):
        render_error(api, id)
        return None
    return out


@click.command()
@click.argument('user', type=str, default='ternviz')
def bot(user):
    load_dotenv('.secrets')
    consumer_key = os.environ.get('CONSUMER_KEY')
    consumer_secret = os.environ.get('CONSUMER_SECRET')
    access_token = os.environ.get('ACCESS_TOKEN')
    access_token_secret = os.environ.get('ACCESS_TOKEN_SECRET')
    auth = tweepy.OAuthHandler(consumer_key, consumer_secret)
    auth.set_access_token(access_token, access_token_secret)

    api = tweepy.API(auth)
    my_user_id = api.get_user(screen_name='ternviz').id

    # Subclass Stream to print IDs of Tweets received
    class IDPrinter(tweepy.Stream):

        def on_status(self, status):
            print('Processing', status.id)
            if 'render' not in status.text:
                print('Determined to be reply with no render in body')
                return
            if status.user.id == my_user_id:
                print('Determined to be self-reply')
                return
            start = time.time_ns()
            pdb_file = make_pdb(api, status.id, status.text)
            if not pdb_file:
                return
            init_status = api.update_status(status="These SMILES look great 😁 Working on rendering 🎨... This can take a few minutes. I will delete this message when I'm done.",
                                            in_reply_to_status_id=status.id, auto_populate_reply_metadata=True)
            movie = render(api, pdb_file, status.id)
            api.destroy_status(init_status.id)
            if not movie:
                return
            print('Completed Rendering')
            media = api.media_upload(movie)
            duration = (time.time_ns() - start) / (10**9)
            init_status = api.update_status(status=f"Rendered in {duration:.2f} seconds", media_ids=[media.media_id],
                                            in_reply_to_status_id=status.id, auto_populate_reply_metadata=True)

    # Initialize instance of the subclass
    printer = IDPrinter(
        consumer_key, consumer_secret,
        access_token, access_token_secret
    )

    # Filter realtime Tweets by keyword
    print('Set-up done, starting bot')

    @retry
    def start():
        printer.filter(track=["@ternviz"],  threaded=True)
    start()
