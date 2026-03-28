# Tern Viz
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/whitead/ternviz/blob/main/colab/ternviz.ipynb)

Convert SMILES to a 3D rotating video. Check-out the Twitter bot [@ternviz](https://twitter.com/ternviz)!

## Install

You need to install both [VMD](https://www.ks.uiuc.edu/Research/vmd/) and `ffmpeg` as system dependencies. After that, sync the project environment with `uv`:

```sh
uv sync
```

To include the optional Twitter bot dependencies too:

```sh
uv sync --extra bot
```

## Usage

From a repo checkout, run commands through `uv run`:

It can be used to render stills, rotating videos, comparisons of molecules, and compilation videos of multiple proteins.

### Render a molecule

To render a molecule:

```sh
uv run ternviz CCCO
```

You may need to specify the path to VMD or ffmpeg:

```sh
uv run ternviz CO --vmd /path/to/vmd/executable
```

On macOS, VMD 2 is often installed as an app bundle. `ternviz` can resolve either
the app bundle itself or the nested executable:

```sh
uv run ternviz CO --vmd /Applications/VMD2b1.app
uv run ternviz CO --vmd /Applications/VMD2b1.app/Contents/vmd2b1/vmd_MACOSXARM64
```

## Deploy To Modal

The Modal deployment lives in `ternviz-modal.py`.

Even if you deploy from a Mac, the Modal image needs the Linux VMD archive, not
the macOS app bundle. Point `TERNVIZ_VMD_ARCHIVE_PATH` at a local copy of
`vmd-2.0.0.bin.LINUXAMD64.tar` before deploying. If you leave it unset,
`ternviz-modal.py` will look for `~/Downloads/vmd-2.0.0.bin.LINUXAMD64.tar`.

If you have not configured Modal on this machine yet:

```sh
uv run python -m modal setup
```

Create the bearer-token secret used by the FastAPI endpoint:

```sh
uv run python -m modal secret create web-auth-token AUTH_TOKEN=replace-with-a-random-token
```

Deploy the app:

```sh
export TERNVIZ_VMD_ARCHIVE_PATH=/absolute/path/to/vmd-2.0.0.bin.LINUXAMD64.tar
uv run python -m modal deploy ternviz-modal.py
```

Notes:

- The first deploy is relatively slow because Modal uploads the VMD tarball,
  installs system packages, and builds VMD inside the image.
- The render function requests an `L4` GPU and uses VMD 2's
  `TachyonLOptiXInternal` renderer.
- The deploy output prints the app URL. The render endpoint is `POST /render/`.

Example smoke test with a real PDB:

```sh
curl -X POST "https://<your-modal-url>/render/" \
  -H "Authorization: Bearer $AUTH_TOKEN" \
  -F "pdb_file=@/tmp/1CRN.pdb" \
  -F "message=1CRN" \
  -F "target_total_frames=60" \
  -o ternviz.mp4
```

The endpoint accepts `.pdb`, `.gz`, `.tar`, and `.tar.gz` uploads.

You can modify a few things too:

```sh
uv run ternviz CCCO --name "my molecule" --color white
```

try `--low-quality` to render quickly.

### Render a protein

This will do a still frame
```sh
uv run ternviz-pdb 1A1L
```

To do a rotating movie, specify number of frames (at 60 fps)

```sh
uv run ternviz-pdb 1A1L --frames 60
```

You can pass in multiple PDBs

```sh
uv run ternviz-pdb my.pdb other.pdb
```

You can specify how the structure is colored and other details

```sh
uv run ternviz-pdb 1A1L --frames 60 --color white --scolor Chain
```

### Aligning PDBs

You can align structures before rendering

```sh
uv run ternviz-align ref.pdb *.pdb
```

by default it aligns on the protein. You can also change selection string, using [MDAnalysis Selection String Syntax](https://docs.mdanalysis.org/stable/documentation_pages/selections.html)

```sh
uv run ternviz-align ref.pdb --sel "chainid B" *.pdb
```
## Example Render

https://user-images.githubusercontent.com/908389/147892964-c15ac0fd-44be-4473-b82a-799260e0f373.mp4
