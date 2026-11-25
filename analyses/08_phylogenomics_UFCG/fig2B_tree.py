#!/usr/bin/env python3
"""Figure 2B (revised): UFCG phylogenomic tree of Naganishia, with UFBoot support at nodes.

Tree: IQ-TREE 3.1.1 on the UFCG concatenated protein alignment (58 genes, 43,635 aa),
VT+F+I+R3 (ModelFinder, BIC), 1,000 UFBoot / 1,000 SH-aLRT; rooted on Aspergillus clavatus.
Node labels in the Newick are 'SH-aLRT/UFBoot'; only UFBoot is drawn.
The outgroup branch is drawn truncated (break mark) and its true length is printed.

Usage: python3 fig2B_tree.py ng1.treefile <out_prefix>
"""
import sys
from io import StringIO
from Bio import Phylo
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import glob
import matplotlib.font_manager as fm
for _f in glob.glob("/usr/share/fonts/truetype/msttcorefonts/[Aa]rial*.ttf"):
    fm.fontManager.addfont(_f)

tree_file, out = sys.argv[1], sys.argv[2]
OUTGROUP = "GCA_000002715.1_ASM271v1_genomic"
NAMES = {  # (italic part, plain part)
    "GCA_000002715.1_ASM271v1_genomic": ("Aspergillus clavatus", "NRRL 1 (GCA_000002715.1)"),
    "GCA_022813745.1_ASM2281374v1_genomic": ("Naganishia kalamii", "FJI-L2-BK-P3 (GCA_022813745.1)"),
    "GCA_001599735.1_JCM_2334_assembly_v001_genomic": ("Naganishia albida", "JCM 2334 (GCA_001599735.1)"),
    "GCA_015708705.1_Nagvi1_genomic": ("Naganishia vishniacii", "ANT03-052 (GCA_015708705.1)"),
    "GCA_012922715.1_ASM1292271v1_genomic": ("Naganishia adeliensis", "IF1SW-F1 (GCA_012922715.1)"),
    "GCA_012922685.1_ASM1292268v1_genomic": ("Naganishia", "sp. IF7SW-B1 (GCA_012922685.1)"),
    "Naganishia_sp_asm.scaffolds": ("Naganishia adeliensis", "Ng1 (this study)"),
    "GCA_013423385.1_N_liquefaciens_N6_1.0_genomic": ("Naganishia liquefaciens", "N6 (GCA_013423385.1)"),
    "GCA_037832395.1_ASM3783239v1_genomic": ("Naganishia diffluens", "bongo (GCA_037832395.1)"),
}
FOCAL = "Naganishia_sp_asm.scaffolds"

tree = Phylo.read(tree_file, "newick")
# IQ-TREE writes the tree with a basal trifurcation (outgroup, clade, clade); keep that root
# and just put the outgroup on top and ladderize the ingroup.
root = tree.root
og = [c for c in root.clades if c.name == OUTGROUP][0]
ing = [c for c in root.clades if c is not og]
for c in ing:
    c.ladderize(reverse=True)
ing.sort(key=lambda c: c.count_terminals(), reverse=True)
root.clades = [og] + ing

def ufboot(cl):
    lab = cl.name if cl.name and "/" in str(cl.name) else None
    return float(lab.split("/")[1]) if lab else None

OUT_DRAW = 0.035                      # drawn length of the truncated outgroup branch
terms = tree.get_terminals()
ypos = {t: i for i, t in enumerate(terms)}
xpos = {}
def place(cl, x):
    bl = cl.branch_length or 0.0
    if cl is og:
        bl = OUT_DRAW
    xpos[cl] = x + bl
    for c in cl.clades:
        place(c, xpos[cl])
place(root, 0.0)
def y_of(cl):
    if cl.is_terminal():
        return ypos[cl]
    ys = [y_of(c) for c in cl.clades]
    ypos[cl] = (min(ys) + max(ys)) / 2
    return ypos[cl]
y_of(root)

plt.rcParams.update({"font.family": "Arial", "font.size": 8, "svg.fonttype": "none"})
INK, INK2, MUTED, LINE = "#0b0b0b", "#52514e", "#8a8984", "#3d3c39"
fig, ax = plt.subplots(figsize=(6.0, 2.7))
BOX = dict(boxstyle="square,pad=0.08", fc="white", ec="none")

def draw(cl, parent_x):
    x, y = xpos[cl], ypos[cl]
    ax.plot([parent_x, x], [y, y], color=LINE, lw=1.0, solid_capstyle="butt", zorder=1)
    if cl is og:                                   # break mark on the truncated branch
        bx = parent_x + (x - parent_x) * 0.55
        for dx in (-0.0011, 0.0011):
            ax.plot([bx + dx - 0.0007, bx + dx + 0.0007], [y + 0.2, y - 0.2], color=LINE, lw=0.9, zorder=3)
        ax.add_patch(plt.Rectangle((bx - 0.0010, y - 0.05), 0.0020, 0.1, fc="white", ec="none", zorder=2))
        ax.text(bx, y - 0.22, f"{cl.branch_length:.2f}", ha="center", va="bottom", fontsize=6.5, color=INK2)
    if not cl.is_terminal():
        ys = [ypos[c] for c in cl.clades]
        ax.plot([x, x], [min(ys), max(ys)], color=LINE, lw=1.0, zorder=1)
        s = ufboot(cl)
        if s is not None:
            ax.text(x - 0.0012, y, f"{s:.0f}", ha="right", va="center", fontsize=6.5,
                    color=INK2, bbox=BOX, zorder=4)
        for c in cl.clades:
            draw(c, x)
draw(root, 0.0)

# tip labels in one aligned column, dotted leaders from each tip
xmax = max(xpos[t] for t in terms)
xlab = xmax + 0.006
for t in terms:
    it, pl = NAMES[t.name]
    if xlab - xpos[t] > 0.002:
        ax.plot([xpos[t] + 0.001, xlab - 0.001], [ypos[t], ypos[t]], color="#b5b4ae", lw=0.6,
                ls=(0, (1, 1.5)), zorder=0)
    focal = t.name == FOCAL
    ax.text(xlab, ypos[t], f"$\\it{{{it.replace(' ', '~')}}}$ " + pl, va="center", ha="left",
            fontsize=7.5, color=INK, fontweight="bold" if focal else "normal")
    if focal:
        ax.plot(xpos[t], ypos[t], marker="o", ms=4.5, color="#e34948", mec="white", mew=0.6, zorder=5)

# scale bar
sb = 0.02
y0 = len(terms) - 0.15
ax.plot([0.0, sb], [y0, y0], color=LINE, lw=1.0)
ax.text(sb / 2, y0 + 0.12, f"{sb} substitutions/site", ha="left" if False else "center", va="top",
        fontsize=6.5, color=INK2)

ax.set_xlim(-0.003, xlab + 0.12)
ax.set_ylim(len(terms) + 0.45, -0.6)
ax.axis("off")
fig.tight_layout()
fig.savefig(out + ".svg")
fig.savefig(out + ".png", dpi=300)
print("written", out + ".svg", out + ".png")
