#!/usr/bin/env python3
"""Composite Figures 1 and 2 for Vazquez-Tapia et al. (N. adeliensis Ng1), as SVG (editable text) and 300-dpi PNG.

Figure 1  A telomeres per nuclear scaffold (Fig1A_telomeres_table.tsv, from fig1A_telomeres.py with ONT evidence)
          -> written as Figure1_panelA; compose_figure1.py adds the original mitochondrial map (panel B)
Figure 2  A ANIm identity heatmap (pyANI, anioutput7), genomes in tree order
          B UFCG maximum-likelihood tree (IQ-TREE 3.1.1; ng1.treefile), UFBoot at nodes
          C dot plot Ng1 vs IF1SW-F1 (NUCmer 3.1 coords) + detail of contig JAAZPZ010000023.1
Usage (from ncbi/figuras_revision/): python3 make_composites.py compuestas/
"""
import sys, os, re, csv, glob, math, collections
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.colors import LinearSegmentedColormap
from Bio import Phylo
for _f in glob.glob("/usr/share/fonts/truetype/msttcorefonts/[Aa]rial*.ttf"):
    fm.fontManager.addfont(_f)

OUT = sys.argv[1] if len(sys.argv) > 1 else "compuestas/"
os.makedirs(OUT, exist_ok=True)
BASE = "/scratch2/Project_ClauMartinez_2024_04_25_12_22_58"
NC = BASE + "/ncbi"
plt.rcParams.update({"font.family": "Arial", "font.size": 7, "svg.fonttype": "none",
                     "axes.linewidth": 0.6, "pdf.fonttype": 42, "mathtext.fontset": "custom",
                     "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold"})
INK, INK2, MUTED, LINE, GRID = "#0b0b0b", "#52514e", "#8a8984", "#3d3c39", "#ecebe7"
BLUE, ORANGE, VIOLET, RED = "#2a78d6", "#eb6834", "#6250d6", "#e34948"
BAR_T2T, BAR_ONE, BAR_NONE = BLUE, "#9cc0ec", "#dcdbd6"
SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]

def letter(fig, x, y, s):
    fig.text(x, y, s, fontsize=11, fontweight="bold", color=INK, ha="left", va="top")

def it(s):                         # italic run for mathtext (spaces kept)
    return "$\\it{" + s.replace(" ", "~") + "}$"

# ======================================================================================
# Figure 1A — telomeres
# ======================================================================================
MINPLOT = 100000          # scaffolds shown in Fig. 1A (>=100 kb: 25 scaffolds, 98.5% of the nuclear assembly)
def panel_telomeres(ax):
    rows = []
    for r in csv.DictReader(open("Fig1A_telomeres_table.tsv"), delimiter="\t"):
        rows.append(dict(n=int(r["scaffold"].split("_")[-1]), L=int(r["length_bp"]) / 1e6,
                         e5=r["tel5_evidence"], e3=r["tel3_evidence"], cls=r["class"]))
    rows = [r for r in rows if r["L"] * 1e6 >= MINPLOT]
    rows.sort(key=lambda r: r["n"])
    cnt = collections.Counter(r["cls"] for r in rows)
    n_asm = sum((r["e5"] == "assembly") + (r["e3"] == "assembly") for r in rows)
    n_ont = sum((r["e5"] == "ONT_reads_only") + (r["e3"] == "ONT_reads_only") for r in rows)
    for i, r in enumerate(rows):
        col = {"T2T": BAR_T2T, "one_end": BAR_ONE, "none": BAR_NONE}[r["cls"]]
        ax.plot([0, r["L"]], [i, i], color=col, lw=2.6, solid_capstyle="round", zorder=1)
        for ev, x, c in ((r["e5"], 0, ORANGE), (r["e3"], r["L"], VIOLET)):
            if ev == "assembly":
                ax.plot(x, i, "o", ms=4.3, color=c, mec="white", mew=0.6, zorder=3)
            elif ev == "ONT_reads_only":
                ax.plot(x, i, "o", ms=3.6, mfc="white", mec=c, mew=1.1, zorder=3)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([str(r["n"]) for r in rows], fontsize=6.2, color=INK2)
    for lab, r in zip(ax.get_yticklabels(), rows):
        if r["cls"] == "T2T":
            lab.set_fontweight("bold"); lab.set_color(INK)
    ax.set_ylim(len(rows) - 0.3, -0.8)
    ax.set_xlim(-0.06, 2.36)
    ax.set_xlabel("Scaffold length (Mb)", color=INK)
    ax.set_ylabel("Nuclear scaffold", color=INK)
    ax.tick_params(axis="both", colors=INK2, width=0.6, length=2.5, labelsize=6.5)
    ax.tick_params(axis="y", length=0, labelsize=6.2)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(MUTED)
    ax.grid(axis="x", color=GRID, lw=0.6, zorder=0); ax.set_axisbelow(True)
    h = [Line2D([], [], marker="o", ls="", ms=4.3, color=ORANGE, mec="white", label="5′ telomere (CCCCTAA)$_n$"),
         Line2D([], [], marker="o", ls="", ms=4.3, color=VIOLET, mec="white", label="3′ telomere (TTAGGGG)$_n$"),
         Line2D([], [], marker="o", ls="", ms=4.3, color=INK2, mec="white", label=f"In assembly ({n_asm} ends)"),
         Line2D([], [], marker="o", ls="", ms=3.6, mfc="white", mec=INK2, mew=1.1,
                label=f"In ONT reads only ({n_ont} ends)"),
         Line2D([], [], color=BAR_T2T, lw=2.6, label=f"Telomere at both ends ({cnt['T2T']})"),
         Line2D([], [], color=BAR_ONE, lw=2.6, label=f"Telomere at one end ({cnt['one_end']})"),
         Line2D([], [], color=BAR_NONE, lw=2.6, label=f"No terminal telomere ({cnt['none']})")]
    ax.legend(handles=h, loc="lower right", frameon=False, fontsize=6.2, handlelength=1.4,
              borderaxespad=0.1, labelcolor=INK)

# ======================================================================================
# Figure 2B — tree
# ======================================================================================
OUTGROUP = "GCA_000002715.1_ASM271v1_genomic"
FOCAL = "Naganishia_sp_asm.scaffolds"
TIPS = {  # tree tip -> (italic, plain, ANI id)
    "GCA_000002715.1_ASM271v1_genomic": ("Aspergillus clavatus", "NRRL 1 (GCA_000002715.1)", None),
    "GCA_022813745.1_ASM2281374v1_genomic": ("N. kalamii", "FJI-L2-BK-P3 (GCA_022813745.1)", "N_kalamii"),
    "GCA_001599735.1_JCM_2334_assembly_v001_genomic": ("N. albida", "JCM 2334 (GCA_001599735.1)", "N_albida"),
    "GCA_015708705.1_Nagvi1_genomic": ("N. vishniacii", "ANT03-052 (GCA_015708705.1)", "N_vishniacii"),
    "GCA_012922715.1_ASM1292271v1_genomic": ("N. adeliensis", "IF1SW-F1 (GCA_012922715.1)", "N_adeliensis"),
    "GCA_012922685.1_ASM1292268v1_genomic": ("Naganishia", "sp. IF7SW-B1 (GCA_012922685.1)", "N_sp.IF7SW-B1"),
    "Naganishia_sp_asm.scaffolds": ("N. adeliensis", "Ng1 (this study)", "N_sp"),
    "GCA_013423385.1_N_liquefaciens_N6_1.0_genomic": ("N. liquefaciens", "N6 (GCA_013423385.1)", "N_liquefaciens"),
    "GCA_037832395.1_ASM3783239v1_genomic": ("N. diffluens", "bongo (GCA_037832395.1)", "N_diffluens"),
}

def load_tree():
    tree = Phylo.read("ng1.treefile", "newick")
    root = tree.root
    og = [c for c in root.clades if c.name == OUTGROUP][0]
    ing = [c for c in root.clades if c is not og]
    for c in ing: c.ladderize(reverse=True)
    ing.sort(key=lambda c: c.count_terminals(), reverse=True)
    root.clades = [og] + ing
    return tree, og

def panel_tree(ax):
    tree, og = load_tree()
    OUT_DRAW = 0.035
    terms = tree.get_terminals()
    ypos = {t: i for i, t in enumerate(terms)}; xpos = {}
    def place(cl, x):
        bl = OUT_DRAW if cl is og else (cl.branch_length or 0.0)
        xpos[cl] = x + bl
        for c in cl.clades: place(c, xpos[cl])
    place(tree.root, 0.0)
    def y_of(cl):
        if cl.is_terminal(): return ypos[cl]
        ys = [y_of(c) for c in cl.clades]; ypos[cl] = (min(ys) + max(ys)) / 2; return ypos[cl]
    y_of(tree.root)
    BOX = dict(boxstyle="square,pad=0.06", fc="white", ec="none")
    def draw(cl, px):
        x, y = xpos[cl], ypos[cl]
        ax.plot([px, x], [y, y], color=LINE, lw=0.9, solid_capstyle="butt", zorder=1)
        if cl is og:
            bx = px + (x - px) * 0.55
            ax.add_patch(plt.Rectangle((bx - 0.0010, y - 0.06), 0.0020, 0.12, fc="white", ec="none", zorder=2))
            for dx in (-0.0011, 0.0011):
                ax.plot([bx + dx - 0.0007, bx + dx + 0.0007], [y + 0.22, y - 0.22], color=LINE, lw=0.8, zorder=3)
            ax.text(bx, y - 0.25, f"{cl.branch_length:.2f}", ha="center", va="bottom", fontsize=5.8, color=INK2)
        if not cl.is_terminal():
            ys = [ypos[c] for c in cl.clades]
            ax.plot([x, x], [min(ys), max(ys)], color=LINE, lw=0.9, zorder=1)
            lab = cl.name if cl.name and "/" in str(cl.name) else None
            if lab:
                ax.text(x - 0.0012, y, f"{float(lab.split('/')[1]):.0f}", ha="right", va="center",
                        fontsize=5.8, color=INK2, bbox=BOX, zorder=4)
            for c in cl.clades: draw(c, x)
    draw(tree.root, 0.0)
    xmax = max(xpos[t] for t in terms); xlab = xmax + 0.005
    for t in terms:
        ital, plain, _ = TIPS[t.name]
        if xlab - xpos[t] > 0.002:
            ax.plot([xpos[t] + 0.001, xlab - 0.001], [ypos[t]] * 2, color="#b5b4ae", lw=0.5, ls=(0, (1, 1.5)), zorder=0)
        foc = t.name == FOCAL
        ax.text(xlab, ypos[t], it(ital) + " " + plain.split(" (")[0] + (" (this study)" if "this study" in plain else ""),
                va="center", ha="left", fontsize=6.5, color=INK,
                fontweight="bold" if foc else "normal")
        if foc:
            ax.plot(xpos[t], ypos[t], "o", ms=3.8, color=RED, mec="white", mew=0.5, zorder=5)
    y0 = len(terms) - 0.1
    ax.plot([0, 0.02], [y0, y0], color=LINE, lw=0.9)
    ax.text(0.01, y0 + 0.15, "0.02 substitutions/site", ha="center", va="top", fontsize=5.8, color=INK2)
    ax.set_xlim(-0.003, xlab + 0.112); ax.set_ylim(len(terms) + 0.55, -0.6); ax.axis("off")
    return [TIPS[t.name][2] for t in terms if TIPS[t.name][2]]

# ======================================================================================
# Figure 2A — ANI heatmap
# ======================================================================================
def panel_ani(ax, cax, order):
    tab = BASE + "/genomescanner/JolyTree_in/anioutput7/ANIm_percentage_identity.tab"
    lines = [l.rstrip("\n").split("\t") for l in open(tab)]
    cols = lines[0][1:]; M = {r[0]: dict(zip(cols, map(float, r[1:]))) for r in lines[1:]}
    lab = {v[2]: (v[0], v[1].split(" (")[0]) for v in TIPS.values() if v[2]}
    A = np.array([[M[r][c] * 100 for c in order] for r in order])
    cmap = LinearSegmentedColormap.from_list("seq", SEQ)
    im = ax.imshow(A, cmap=cmap, vmin=84, vmax=100, aspect="equal")
    for i in range(len(order)):
        for j in range(len(order)):
            v = A[i, j]
            ax.text(j, i, f"{v:.2f}" if v < 99.995 else "100", ha="center", va="center", fontsize=5.0,
                    color="white" if v > 95 else INK)
    names = [it(lab[o][0]) + " " + lab[o][1] for o in order]
    ax.set_xticks(range(len(order))); ax.set_yticks(range(len(order)))
    ax.set_yticklabels(names, fontsize=6.0, color=INK)
    ax.set_xticklabels([lab[o][1].replace("sp. ", "") for o in order], fontsize=6.0, color=INK, rotation=45,
                       ha="right", rotation_mode="anchor")
    for tl in list(ax.get_yticklabels()) + list(ax.get_xticklabels()):
        if "Ng1" in tl.get_text(): tl.set_fontweight("bold"); tl.set_color(INK)
    ax.tick_params(length=0, pad=2)
    for s in ax.spines.values(): s.set_visible(False)
    ax.set_xticks(np.arange(-0.5, len(order)), minor=True); ax.set_yticks(np.arange(-0.5, len(order)), minor=True)
    ax.grid(which="minor", color="white", lw=0.8); ax.tick_params(which="minor", length=0)
    cb = plt.colorbar(im, cax=cax, orientation="horizontal")
    cb.set_label("ANIm identity (%)", fontsize=6.2, color=INK, labelpad=2)
    cb.ax.tick_params(labelsize=5.8, colors=INK2, width=0.5, length=2); cb.outline.set_visible(False)
    return A

# ======================================================================================
# Figure 2C — dot plot
# ======================================================================================
MISJOIN = "JAAZPZ010000023.1"
def read_coords():
    C = BASE + "/ref_tree/genomes/Naganishia_sp_asm.scaffolds.fastaVSGCA_012922715.1_ASM1292271v1_genomic.fna.delta.filter.coords"
    al = []
    for l in open(C):
        if "\t" not in l or not l.strip()[:1].isdigit(): continue
        t = l.rstrip("\n").split("\t"); w = t[0].split(); n = [x for x in w[:-1] if x != "|"]
        al.append(dict(r=w[-1], q=t[1], s1=int(n[0]), e1=int(n[1]), s2=int(n[2]), e2=int(n[3]),
                       L=int(n[4]), LR=int(n[7]), LQ=int(n[8])))
    return al

def panel_dotplot(ax, axd):
    al = read_coords()
    lens = {l.split()[0]: int(l.split()[1]) for l in open(NC + "/Naganishia_sp_asm.scaffolds.fa.fai")}
    scafs = sorted([k for k in lens if not k.endswith("_27")], key=lambda k: int(k.split("_")[-1]))
    xoff = {}; acc = 0
    for k in scafs: xoff[k] = acc; acc += lens[k]
    XT = acc
    byq = collections.defaultdict(list)
    for a in al:
        if a["r"] in xoff: byq[a["q"]].append(a)
    qinfo = {}
    for q, v in byq.items():
        b = max(v, key=lambda a: a["L"])
        qinfo[q] = (xoff[b["r"]] + min(b["s1"], b["e1"]), b["e2"] < b["s2"], v[0]["LQ"])
    order = sorted(qinfo, key=lambda q: qinfo[q][0])
    yoff = {}; acc = 0
    for q in order: yoff[q] = acc; acc += qinfo[q][2]
    YT = acc
    for q in order:
        flip, LQ = qinfo[q][1], qinfo[q][2]
        for a in byq[q]:
            ys, ye = (LQ - a["s2"], LQ - a["e2"]) if flip else (a["s2"], a["e2"])
            x = [(xoff[a["r"]] + a["s1"]) / 1e6, (xoff[a["r"]] + a["e1"]) / 1e6]
            y = [(yoff[q] + ys) / 1e6, (yoff[q] + ye) / 1e6]
            hl = q == MISJOIN
            ax.plot(x, y, color=RED if hl else LINE, lw=1.6 if hl else 0.9, solid_capstyle="round",
                    zorder=3 if hl else 2)
    for k in scafs:
        if lens[k] > 300000:
            ax.axvline(xoff[k] / 1e6, color=GRID, lw=0.4, zorder=0)
    ax.set_xlim(0, XT / 1e6); ax.set_ylim(0, YT / 1e6)
    ax.set_aspect("equal")
    ax.set_xlabel(it("N. adeliensis") + " Ng1 (Mb)", color=INK)
    ax.set_ylabel(it("N. adeliensis") + " IF1SW-F1 (Mb)", color=INK)
    ax.tick_params(colors=INK2, width=0.6, length=2.5, labelsize=6.2)
    for s in ax.spines.values(): s.set_color(MUTED)
    # annotate the off-diagonal segment (contig 23 vs Ng1_4)
    a4 = [a for a in byq[MISJOIN] if a["r"].endswith("_4")][0]
    flip, LQ = qinfo[MISJOIN][1], qinfo[MISJOIN][2]
    xm = (xoff[a4["r"]] + (a4["s1"] + a4["e1"]) / 2) / 1e6
    ym = (yoff[MISJOIN] + ((LQ - (a4["s2"] + a4["e2"]) / 2) if flip else (a4["s2"] + a4["e2"]) / 2)) / 1e6
    ax.annotate("contig JAAZPZ010000023.1", xy=(xm, ym), xytext=(xm + 1.2, ym + 2.6), fontsize=6.0, color=INK,
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.5), ha="left", va="bottom")
    # ---- detail: contig 23 against Ng1_4 and Ng1_13
    s4, s13 = "Naganishia_adeliensis_Ng1_4", "Naganishia_adeliensis_Ng1_13"
    gap = 0.10
    ox = {s4: 0.0, s13: lens[s4] / 1e6 + gap}
    for a in byq[MISJOIN]:
        if a["r"] not in ox: continue
        axd.plot([ox[a["r"]] + a["s1"] / 1e6, ox[a["r"]] + a["e1"] / 1e6], [a["s2"] / 1e6, a["e2"] / 1e6],
                 color=RED, lw=1.6, solid_capstyle="round", zorder=3)
    for k in (s4, s13):
        axd.axvspan(ox[k], ox[k] + lens[k] / 1e6, color="#f6f5f2", zorder=0, lw=0)
    rep4, rep13, rep23 = 423.24e-3, 255.81e-3, 360.79e-3            # repeat midpoints (Mb)
    for x in (ox[s4] + rep4, ox[s13] + rep13):
        axd.axvline(x, color=MUTED, lw=0.6, ls=(0, (2, 2)), zorder=1)
    axd.axhline(rep23, color=MUTED, lw=0.6, ls=(0, (2, 2)), zorder=1)
    axd.text(ox[s4] + 0.72, rep23 + 0.015, "junction: 0.8-kb repeat present\nin Ng1 scaffolds 4 and 13 (dashed)",
             fontsize=5.8, color=INK2, ha="center", va="bottom")
    xt = [0, 0.5, 1.0]; xt13 = [0, 0.5]
    axd.set_xticks([ox[s4] + v for v in xt] + [ox[s13] + v for v in xt13])
    axd.set_xticklabels([f"{v:g}" for v in xt] + [f"{v:g}" for v in xt13])
    for k, nm in ((s4, "Ng1 scaffold 4"), (s13, "Ng1 scaffold 13")):
        axd.text(ox[k] + lens[k] / 2e6, -0.17, nm, ha="center", va="top", fontsize=6.2, color=INK,
                 transform=axd.get_xaxis_transform())
    axd.set_xlim(-0.02, ox[s13] + lens[s13] / 1e6 + 0.02); axd.set_ylim(0, qinfo[MISJOIN][2] / 1e6)
    axd.set_ylabel("IF1SW-F1 contig\nJAAZPZ010000023.1 (Mb)", color=INK, fontsize=6.4)
    axd.tick_params(colors=INK2, width=0.6, length=2.5, labelsize=6.2)
    for s in ("top", "right"): axd.spines[s].set_visible(False)
    for s in ("left", "bottom"): axd.spines[s].set_color(MUTED)
    return dict(contigs=len(order), aligned=sum(a["L"] for a in al))

# ======================================================================================
def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".svg"))
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=300)
    print("written", os.path.join(OUT, name + ".svg"), os.path.join(OUT, name + ".png"))

# ---- Figure 1, panel A only (panel B = original mitochondrial map, assembled by compose_figure1.py)
fig = plt.figure(figsize=(3.80, 4.0))
axA = fig.add_axes([0.14, 0.105, 0.82, 0.86])
panel_telomeres(axA)
letter(fig, 0.015, 0.99, "A")
save(fig, "Figure1_panelA")

# ---- Figure 2 (180 mm wide)
fig = plt.figure(figsize=(7.09, 6.6))
axT = fig.add_axes([0.545, 0.585, 0.45, 0.39])
order = panel_tree(axT)
axH = fig.add_axes([0.195, 0.615, 0.33, 0.36]); cax = fig.add_axes([0.26, 0.535, 0.20, 0.011])
panel_ani(axH, cax, order)
axC = fig.add_axes([0.075, 0.07, 0.40, 0.40]); axD = fig.add_axes([0.62, 0.13, 0.35, 0.29])
dinfo = panel_dotplot(axC, axD)
letter(fig, 0.008, 0.995, "A"); letter(fig, 0.535, 0.995, "B"); letter(fig, 0.008, 0.49, "C")
print("dotplot:", dinfo)
save(fig, "Figure2")
