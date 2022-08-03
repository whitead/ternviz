# Tern Viz
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/whitead/ternviz/blob/main/colab/ternviz.ipynb)

Convert SMILES to a 3D rotating video. Check-out the Twitter bot [@ternviz](https://twitter.com/ternviz)!

## Install

You need to install both [VMD](https://www.ks.uiuc.edu/Research/vmd/) and ffmpeg. After that:

```py
pip install ternviz
```

## Usage

It can be used to render stills, rotating videos, comparisons of molecules, and compilation videos of multiple proteins.

### Render a molecule

To render a molecule:

```sh
ternviz CCCO
```

You may need to specify the path to VMD or ffmpeg:

```sh
ternviz CO --vmd /path/to/vmd/executable
```

You can modify a few things too:

```sh
ternviz CCCO --name "my molecule" --color white
```

try `--low-quality` to render quickly.

### Render a protein

This will do a still frame
```sh
ternviz-pdb 1A1L
```

To do a rotating movie, specify number of frames (at 60 fps)

```sh
ternviz-pdb 1A1L --frames 60
```

You can pass in multiple PDBs

```sh
ternviz-pdb my.pdb other.pdb
```

You can specify how the structure is colored and other details

```sh
ternviz-pdb 1A1L --frames 60 --color white --scolor Chain
```

### Aligning PDBs

You can align structures before rendering

```sh
ternviz-align ref.pdb *.pdb
```

by default it aligns on the protein. You can also change selection string, using [MDAnalysis Selection String Syntax](https://docs.mdanalysis.org/stable/documentation_pages/selections.html)

```sh
ternviz-align ref.pdb --sel "chain B" *.pdb
```
## Example Render

https://user-images.githubusercontent.com/908389/147892964-c15ac0fd-44be-4473-b82a-799260e0f373.mp4
