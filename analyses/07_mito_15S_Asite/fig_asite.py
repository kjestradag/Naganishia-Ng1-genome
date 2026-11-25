#!/usr/bin/env python3
"""Supplementary figure: aminoglycoside-binding A site (h44 decoding site) of the Ng1 mitochondrial SSU rRNA.

Structural alignment: Infernal 1.1.5 cmalign with Rfam RF00177 (SSU_rRNA_bacteria); E. coli numbering
(J01859 position + 1). Input: all.sto (cmalign output). All A-site columns have posterior probability >=0.95.
Usage: python3 fig_asite.py all.sto <out_prefix>
"""
import sys, glob
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.patches import FancyBboxPatch, Rectangle
for _f in glob.glob("/usr/share/fonts/truetype/msttcorefonts/[Aa]rial*.ttf"):
    fm.fontManager.addfont(_f)
plt.rcParams.update({"font.family": "Arial", "font.size": 7, "svg.fonttype": "none", "mathtext.fontset": "custom",
                     "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold"})
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#ecebe7"
KEY_BG, NG1_BG = "#fdebd9", "#eef4fc"
ORANGE, BLUE, RED = "#eb6834", "#2a78d6", "#e34948"

sto, out = sys.argv[1], sys.argv[2]
seqs = {}
for l in open(sto):
    if l.startswith(("#", "//")) or not l.strip(): continue
    n, s = l.split(); seqs[n] = seqs.get(n, "") + s
e = seqs["Ecoli_16S_rrsB"]; col = {}; k = 0
for i, ch in enumerate(e):
    if ch not in ".-": k += 1; col[k + 1] = i          # E. coli standard numbering = J01859 + 1
def res(n, p): return seqs[n][col[p]].upper().replace("T", "U")

ROWS = [("Ecoli_16S_rrsB", "$\\it{E.~coli}$ 16S", "Bacteria"),
        ("Human_cyto_18S", "Human 18S", "Cytosol"),
        ("Scer_mito_15S", "$\\it{S.~cerevisiae}$ 15S", "Mitochondria"),
        ("Cneoformans_mito_rns", "$\\it{C.~neoformans}$ 15S", "Mitochondria"),
        ("Ng1_rns", "$\\it{N.~adeliensis}$ Ng1 15S", "Mitochondria")]
TOP = list(range(1404, 1413)); BOT = list(range(1488, 1498))
KEY = {1408, 1409, 1491}
PAIRS = [(1404, 1497), (1405, 1496), (1406, 1495), (1407, 1494), (1409, 1491), (1410, 1490), (1411, 1489), (1412, 1488)]
WC = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G")}; GU = {("G", "U"), ("U", "G")}
def pair_type(a, b): return "wc" if (a, b) in WC else "gu" if (a, b) in GU else "mm"

fig = plt.figure(figsize=(7.09, 3.5))
# ---------------- panel A: alignment ----------------
ax = fig.add_axes([0.0, 0.0, 0.64, 1.0]); ax.axis("off")
x0, dx, y0, dy = 0.30, 0.0255, 0.78, 0.125
cols = TOP + [None] + BOT
xs = {}
for j, p in enumerate(cols):
    if p is not None: xs[p] = x0 + j * dx
for p in KEY:
    ax.add_patch(Rectangle((xs[p] - dx / 2, y0 - (len(ROWS) - 0.5) * dy), dx, len(ROWS) * dy + 0.035,
                           fc=KEY_BG, ec="none", zorder=0, transform=ax.transAxes))
iy = len(ROWS) - 1
ax.add_patch(Rectangle((x0 - dx / 2 - 0.27, y0 - iy * dy - dy / 2), 0.27 + len(cols) * dx, dy,
                       fc=NG1_BG, ec="none", zorder=0, transform=ax.transAxes))
for p in TOP + BOT:
    ax.text(xs[p], y0 + 0.07, str(p), rotation=90, ha="center", va="bottom", fontsize=5.6,
            color=INK if p in KEY else MUTED, fontweight="bold" if p in KEY else "normal", transform=ax.transAxes)
ax.text((xs[1404] + xs[1412]) / 2, y0 + 0.165, "5′ strand", ha="center", fontsize=6.2, color=INK2, transform=ax.transAxes)
ax.text((xs[1488] + xs[1497]) / 2, y0 + 0.165, "3′ strand", ha="center", fontsize=6.2, color=INK2, transform=ax.transAxes)
for i, (n, lab, comp) in enumerate(ROWS):
    y = y0 - i * dy
    ax.text(0.005, y + 0.012, lab, ha="left", va="center", fontsize=6.6, color=INK, transform=ax.transAxes,
            fontweight="bold" if n == "Ng1_rns" else "normal")
    ax.text(0.005, y - 0.03, comp, ha="left", va="center", fontsize=5.6, color=MUTED, transform=ax.transAxes)
    for p in TOP + BOT:
        r = res(n, p); ref = res("Ecoli_16S_rrsB", p)
        c = INK if r == ref else RED
        ax.text(xs[p], y, r, ha="center", va="center", fontsize=7.2, family="DejaVu Sans Mono", color=c,
                fontweight="bold" if p in KEY else "normal", transform=ax.transAxes)
    a, b = res(n, 1409), res(n, 1491)
    pt = pair_type(a, b)
    summary = f"{res(n,1408)}1408   {a}1409{'–' if pt == 'wc' else '·'}{b}1491"
    ax.text(xs[1497] + 0.03, y, summary, ha="left", va="center", fontsize=6.4, color=INK, transform=ax.transAxes)
ax.text(xs[1497] + 0.03, y0 + 0.07, "Key positions", ha="left", va="bottom", fontsize=6.2, color=INK2,
        transform=ax.transAxes)
ax.text(0.005, y0 - len(ROWS) * dy + 0.01, "Red: differs from $\\it{E.~coli}$; shaded: positions that determine "
        "aminoglycoside binding (E. coli numbering)", fontsize=5.6, color=INK2, transform=ax.transAxes, va="top")
fig.text(0.008, 0.985, "A", fontsize=11, fontweight="bold", va="top")

# ---------------- panel B: A-site ladders ----------------
LADDER = [("Ecoli_16S_rrsB", "$\\it{E.~coli}$ 16S", "bacterial (sensitive)"),
          ("Human_cyto_18S", "Human 18S", "cytosolic (resistant)"),
          ("Ng1_rns", "$\\it{N.~adeliensis}$ Ng1 15S", "mitochondrial")]
axb = fig.add_axes([0.67, 0.04, 0.32, 0.86]); axb.axis("off")
axb.set_xlim(0, 3); axb.set_ylim(-0.6, 10.6)
fig.text(0.655, 0.985, "B", fontsize=11, fontweight="bold", va="top")
LEFT = TOP                                       # 1404 (top) .. 1412 (bottom)
RIGHT = [1497, 1496, 1495, 1494, 1493, 1492, 1491, 1490, 1489, 1488]
yl = {p: 9.6 - i * 1.05 for i, p in enumerate(LEFT)}
yl[1408] = (yl[1407] + yl[1409]) / 2 + 0.0
for i, p in enumerate(LEFT): yl[p] = 9.6 - i * 1.05 - (0.5 if p >= 1409 else 0)
yr = {1497: yl[1404], 1496: yl[1405], 1495: yl[1406], 1494: yl[1407],
      1493: yl[1407] - 0.75, 1492: yl[1407] - 1.5, 1491: yl[1409], 1490: yl[1410], 1489: yl[1411], 1488: yl[1412]}
for k, (n, lab, sub) in enumerate(LADDER):
    xc = 0.5 + k; xL, xR = xc - 0.22, xc + 0.22
    if n == "Ng1_rns":
        axb.add_patch(FancyBboxPatch((xc - 0.47, -0.45), 0.94, 10.95, boxstyle="round,pad=0,rounding_size=0.12",
                                     fc=NG1_BG, ec="none", zorder=0))
    axb.text(xc, 10.35, lab, ha="center", va="bottom", fontsize=6.3, color=INK,
             fontweight="bold" if n == "Ng1_rns" else "normal")
    axb.text(xc, 10.05, sub, ha="center", va="bottom", fontsize=5.6, color=MUTED)
    for a, b in PAIRS:
        ra, rb = res(n, a), res(n, b); t = pair_type(ra, rb); y = yl[a]
        key = (a, b) == (1409, 1491)
        lc = (ORANGE if key else INK2)
        if t == "wc":
            axb.plot([xL + 0.07, xR - 0.07], [y, y], color=lc, lw=1.6 if key else 0.9, solid_capstyle="round", zorder=1)
        elif t == "gu":
            axb.plot([xc], [y], "o", ms=2.2, color=lc, zorder=1)
        else:
            axb.plot([xL + 0.07, xR - 0.07], [y, y], color="#c9c8c2", lw=0.7, ls=(0, (1.2, 1.4)), zorder=1)
    for p in LEFT:
        r = res(n, p); key = p in (1408, 1409)
        axb.text(xL, yl[p], r, ha="center", va="center", fontsize=7.2, family="DejaVu Sans Mono",
                 color=(ORANGE if p == 1408 and r == "A" else BLUE if p == 1408 else INK),
                 fontweight="bold" if key else "normal", zorder=2)
    for p in RIGHT:
        r = res(n, p); xx = xR + (0.16 if p in (1492, 1493) else 0)
        axb.text(xx, yr[p], r, ha="center", va="center", fontsize=7.2, family="DejaVu Sans Mono", color=INK,
                 fontweight="bold" if p == 1491 else "normal", zorder=2)
    if k == 0:
        for p in (1404, 1408, 1409, 1412):
            axb.text(xL - 0.13, yl[p], str(p), ha="right", va="center", fontsize=5.2, color=MUTED)
        for p in (1497, 1492, 1491, 1488):
            axb.text(xR + (0.3 if p in (1492, 1493) else 0.14), yr[p], str(p), ha="left", va="center",
                     fontsize=5.2, color=MUTED)
axb.text(1.5, -0.35, "line: Watson–Crick pair;  dot: G·U pair;  dotted: non-canonical",
         ha="center", va="top", fontsize=5.4, color=INK2)
fig.savefig(out + ".svg"); fig.savefig(out + ".png", dpi=300)
print("written", out + ".svg", out + ".png")
