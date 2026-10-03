"""Draw the home page's TALE-on-DNA figure, pkgdown/assets/crd_dna.png.

The crystal structure of PthXo1 bound to its DNA target (PDB 3UGM; Mak et
al. 2012, https://doi.org/10.1126/science.1216211), seen from the side and
down the DNA axis. PDB entries are in the public domain (CC0).

Colours follow R/palette.R (Paul Tol's muted scheme), as in plot.tales():
N-terminal region wine, successive repeats alternately sand and olive so the
modules stand out, RVD residues indigo spheres, DNA grey (the strand the
RVDs read is the darker one).

Needs PyMOL (open source) and Pillow, for instance in a scratch environment:

    micromamba create -p /tmp/pymol-env -c conda-forge pymol-open-source pillow
    /tmp/pymol-env/bin/pymol -cq pkgdown/crd_figure.py

Run from the package root; the structure is downloaded from the PDB.
"""

import os
import re
import tempfile

from PIL import Image, ImageChops
from pymol import cmd

OUT = "pkgdown/assets/crd_dna.png"
HEIGHT = 450  # final height of both panels, in pixels

cmd.fetch("3ugm", "tale", type="cif", path=tempfile.mkdtemp())
cmd.remove("solvent or inorganic or hydro")

# Chain A is PthXo1, chains B and C the DNA duplex; B is the strand the RVDs
# read, 5' to 3' from the N-terminal end.
three = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
         "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
         "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
         "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
ca = cmd.get_model("chain A and name CA").atom
resi = [int(a.resi) for a in ca]
seq = "".join(three[a.resn] for a in ca)

# Each repeat starts L-T/P-x-x-Q-V-V-A-I-A-S; its RVD is residues 12-13, or
# residue 12 alone when 13 is missing (N*, read off the shorter G-G-K-Q
# stretch). The structure stops inside repeat 23: the final half repeat and
# the C-terminal region are not resolved.
starts = [m.start() for m in re.finditer(r"(?=L[TP]..QVVAIAS)", seq)]

colours = {"wine": "0x882255", "sand": "0xDDCC77", "olive": "0x999933",
           "indigo": "0x332288", "grey": "0xCCCCCC", "dark_grey": "0x999999"}

cmd.hide("everything")
# residues before 191 are a few short, disordered fragments
cmd.show("cartoon", "chain A and resi 191-")
cmd.color(colours["wine"], f"chain A and resi {resi[0]}-{resi[starts[0]] - 1}")
rvds = []
for k, s in enumerate(starts):
    first = resi[s]
    last = resi[starts[k + 1]] - 1 if k + 1 < len(starts) else resi[-1]
    cmd.color(colours["sand" if k % 2 == 0 else "olive"],
              f"chain A and resi {first}-{last}")
    rvd = re.match(r"L.{4}VVAIAS(\w{1,2})G(?:GKQ|$)", seq[s:s + 20]).group(1)
    rvds.append(f"{first + 11}-{first + 10 + len(rvd)}")
cmd.select("rvd", "chain A and resi " + "+".join(rvds))
cmd.show("spheres", "rvd and not name N+C+O")
cmd.color(colours["indigo"], "rvd")
cmd.show("cartoon", "chain B+C")
cmd.color(colours["grey"], "chain C")
cmd.color(colours["dark_grey"], "chain B")

cmd.bg_color("white")
settings = {"ray_opaque_background": 1, "ray_trace_mode": 1,
            "ray_trace_gain": 0.05, "antialias": 2, "ambient": 0.5,
            "specular": 0.2, "ray_shadows": 0, "sphere_scale": 0.9,
            "cartoon_gap_cutoff": 0, "cartoon_ring_mode": 3,
            "cartoon_ladder_mode": 0, "cartoon_side_chain_helper": 1}
for name, value in settings.items():
    cmd.set(name, value)

# DNA axis horizontal, N-terminal end on the left
cmd.orient("chain B+C")
cmd.turn("z", 180)
side_view = cmd.get_view()


def panel(turns, path):
    cmd.set_view(side_view)
    for axis, angle in turns:
        cmd.turn(axis, angle)
    cmd.zoom("all", 2)
    cmd.png(path, width=2400, height=1600, ray=1)  # drawn large, then scaled down
    im = Image.open(path).convert("RGB")
    im = im.crop(ImageChops.difference(im, Image.new("RGB", im.size, "white")).getbbox())
    return im.resize((round(im.width * HEIGHT / im.height), HEIGHT), Image.LANCZOS)


tmp = tempfile.mkdtemp()
side = panel([], f"{tmp}/side.png")
end = panel([("y", 90)], f"{tmp}/end.png")  # looking from the N-terminal end
gap = HEIGHT // 8
figure = Image.new("RGB", (side.width + gap + end.width, HEIGHT), "white")
figure.paste(side, (0, 0))
figure.paste(end, (side.width + gap, 0))
os.makedirs(os.path.dirname(OUT), exist_ok=True)
# 256 colours: no visible banding at this size, a third of the file size
figure = figure.quantize(colors=256, method=Image.Quantize.MEDIANCUT,
                         dither=Image.Dither.FLOYDSTEINBERG)
figure.save(OUT, optimize=True)
