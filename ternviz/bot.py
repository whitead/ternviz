import tweepy
import click
from dotenv import load_dotenv
from rdkit import Chem
import os
import ternviz
import multiprocessing
import time

TIMEOUT_COORDS = 10
TIMEOUT_RENDER = 1000


def smiles_error(api, id, s):
    return api.update_status(status=f"Failed to generate coordinates for your SMILES 😢\n Your SMILES were '{s}'",
                             in_reply_to_status_id=id, auto_populate_reply_metadata=True)


def render_error(api, id):
    api.update_status(status=f"Was able to make coordinates, but rendering took too long 😢 It might be that @andrewwhite01 is broke and cannot afford enough server time🤏🏻",
                      in_reply_to_status_id=id, auto_populate_reply_metadata=True)


def full_text(api, status_id):
    status = api.get_status(status_id, tweet_mode="extended")
    try:
        text = status.retweeted_status.full_text
    except AttributeError:  # Not a Retweet
        text = status.full_text
    return text


def make_pdb(api, id, smiles):
    smiles_status = ternviz.check_smiles(smiles)
    if smiles_status:
        se = smiles_error(api, id, smiles)
        api.update_status(status=smiles_status[:179],
                          in_reply_to_status_id=se.id, auto_populate_reply_metadata=True)
        return None, None

    print('Determined SMILES are', smiles)
    p = multiprocessing.Process(
        target=ternviz.gen_coords, args=(smiles, str(id)))
    p.start()
    p.join(TIMEOUT_COORDS)
    if p.is_alive():
        p.terminate()
        p.join()
        smiles_error(api, id, smiles)
        return None, None
    out = os.path.join('/var/tmp', f'{id}.pdb')
    if not os.path.exists(out):
        smiles_error(api, id, smiles)
        return None, None
    return out, ternviz.get_name(smiles)


def do_render(api, pdb_file, id, short_name, width):
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
        target=ternviz.movie, args=(str(id), short_name))
    p2.start()
    p2.join(TIMEOUT_RENDER)
    out = os.path.join('/var/tmp', f'{id}.mp4')
    if not os.path.exists(out):
        render_error(api, id)
        return None
    return out


def render(api, status, smiles, width=800):
    pdb_file, short_name = make_pdb(
        api, status.id, smiles)
    if not pdb_file:
        return None, None
    try:
        init_status = api.update_status(status=f"These SMILES look great 😁 Working on rendering {short_name} 🎨... This can take a few minutes. I will delete this message when I'm done.",
                                        in_reply_to_status_id=status.id, auto_populate_reply_metadata=True)
    except tweepy.os.error.HTTPException as e:
        return None, None
    movie = do_render(api, pdb_file, str(status.id) +
                      short_name, short_name, width=width)
    api.destroy_status(init_status.id)
    return movie, short_name


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
            call_type = ''
            if 'render' in status.text:
                call_type = 'render'
            elif 'compare' in status.text:
                call_type = 'compare'
            if not call_type:
                print('Determined to be reply with no render in body')
                return
            if status.user.id == my_user_id:
                print('Determined to be self-reply')
                return
            start = time.time_ns()
            if call_type == 'render':
                smiles = full_text(api, status.id).split('render')[-1].strip()
                movie, short_name = render(api, status, smiles)
            elif call_type == 'compare':
                smiles = list(full_text(api, status.id).split(
                    'compare')[-1].split())
                m1, s1 = render(api, status, smiles[0], width=800 // 2)
                m2, s2 = render(api, status, smiles[1], width=800 // 2)
                movie = ternviz.multiplex([m1, m2], status.id)
                short_name = s1 + ' vs. ' + s2
            if not movie:
                return
            print('Completed Rendering')
            media = api.media_upload(movie)
            duration = (time.time_ns() - start) / (10**9)
            try:
                api.update_status(status=f'Rendered {short_name} in {duration:.2f} seconds.', media_ids=[media.media_id],
                                  in_reply_to_status_id=status.id, auto_populate_reply_metadata=True)
            except tweepy.os.error.HTTPException as e:
                return

    # Initialize instance of the subclass
    printer = IDPrinter(
        consumer_key, consumer_secret,
        access_token, access_token_secret
    )

    # Filter realtime Tweets by keyword
    print('Set-up done, starting bot')

    while True:
        try:
            p = printer.filter(track=['@ternviz'],  threaded=True)
            p.join()
        except:
            time.sleep(1.0)
            print('Died, restarting')
            continue
