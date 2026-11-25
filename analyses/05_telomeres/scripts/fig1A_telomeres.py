#!/usr/bin/env python3
"""Figure 1A (revised): telomeric arrays at the ends of the N. adeliensis Ng1 nuclear scaffolds.

Strict criterion: a scaffold end carries a telomere when a tandem array of >=4 telomeric
units (C-rich (C2-5TAA)n at the 5' end, G-rich (TTAG2-5)n at the 3' end; units <=3 bp apart,
runs <=50 bp apart merged) starts within MAXDIST bp of that end. Arrays further inside are
interstitial and are not counted. The mitochondrial contig (27) is excluded.
ONT evidence (3rd and 4th arguments):
  ont_end_support.tsv     soft-clip analysis; defines 'trim' ends, where the terminal segment beyond a telomeric
                          clip junction lacks ONT support (the boundary is then the junction, not the scaffold end)
  audit_read_termini.tsv  where ONT reads terminate. An end is ONT-confirmed when >=5 reads end in telomeric
                          repeat, >=3 of them uniquely placed (MAPQ 60, >=3 kb aligned), they outnumber reads that
                          continue >=500 bp past the boundary without telomere, and their median tip lies
                          <=1 kb beyond the boundary.
Classification: 'asm' = assembled array at the end (not inside a trim segment) and not contradicted by reads;
'ont' = ONT-confirmed without an assembled array; assembled arrays contradicted by reads are not counted
('conflicting'); assembled arrays without informative reads (<5) are kept as 'unverifiable'.

Usage: python3 fig1A_telomeres.py <scaffolds.fa> <out_prefix> [ont_end_support.tsv] [audit_read_termini.tsv]
"""
import sys, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import glob
import matplotlib.font_manager as fm
for _f in glob.glob("/usr/share/fonts/truetype/msttcorefonts/[Aa]rial*.ttf"):
    fm.fontManager.addfont(_f)
from matplotlib.lines import Line2D
from telo_scan import scan

MAXDIST = 200
MITO = "Naganishia_adeliensis_Ng1_27"

fasta, out = sys.argv[1], sys.argv[2]
ont = {}
if len(sys.argv) > 3:
    for d in csv.DictReader(open(sys.argv[3]), delimiter="\t"):
        ont[(int(d["scaffold"].split("_")[-1]), d["end"])] = d
aud = {}
if len(sys.argv) > 4:
    for d in csv.DictReader((l for l in open(sys.argv[4]) if not l.startswith("#")), delimiter="\t"):
        aud[(int(d["scaffold"]), d["end"])] = d
def verify(n, e):
    d = aud.get((n, e))
    if d is None: return ""
    tel, cont, uq = int(d["tel_term"]), int(d["continue"]), int(d["tel_unique_mapq60_aln3kb"])
    ext = float(d["tel_term_ext_beyond_boundary_median"]) if d["tel_term_ext_beyond_boundary_median"] else None
    if tel + cont < 5: return "unverifiable"
    if tel >= 5 and uq >= 3 and tel >= cont and ext is not None and ext <= 1000: return "confirmed"
    return "conflicting"
rows = []
for r in scan(fasta):
    if r["id"] == MITO:
        continue
    t5 = [a for a in r["a5"] if a[0] <= MAXDIST]
    t3 = [a for a in r["a3"] if a[0] <= MAXDIST]
    rows.append(dict(n=int(r["id"].split("_")[-1]), length=r["len"],
                     t5=t5[0] if t5 else None, t3=t3[0] if t3 else None))
for r in rows:   # evidence per end: 'asm' (in assembly), 'ont' (long reads only) or None
    for e, k in (("5", "t5"), ("3", "t3")):
        o = ont.get((r["n"], e), {}).get("ont_support", "none")
        r["o" + e] = ont.get((r["n"], e), {})
        v = verify(r["n"], e) if aud else ""
        asm = bool(r[k]) and o != "trim"
        r["x" + e] = None
        if asm:
            if v == "conflicting":
                r["e" + e], r["x" + e] = None, "asm"          # assembled array contradicted by reads
            else:
                r["e" + e] = "asm"
        elif v == "confirmed" or (not aud and o in ("extend", "trim")):
            r["e" + e] = "ont"
        else:
            r["e" + e] = None
            if o in ("extend", "trim"):
                r["x" + e] = "ont"                             # clip-based call not confirmed by read termini
        r["v" + e] = v if (r["e" + e] or r["x" + e]) else ""

# table for the supplement / text
with open(out + "_table.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["scaffold", "length_bp", "tel5_start_from_end", "tel5_units", "tel5_bp",
                "tel3_start_from_end", "tel3_units", "tel3_bp",
                "tel5_evidence", "tel5_ont_reads", "tel5_ont_mode", "tel5_unsupported_terminal_bp",
                "tel3_evidence", "tel3_ont_reads", "tel3_ont_mode", "tel3_unsupported_terminal_bp",
                "tel5_verification", "tel3_verification", "class"])
    for r in rows:
        cls = ("T2T" if r["e5"] and r["e3"] else "one_end" if (r["e5"] or r["e3"]) else "none")
        r["cls"] = cls
        f = lambda a: [a[0], a[2], a[1] - a[0]] if a else ["", "", ""]
        def g(e):
            o = r["o" + e]; m = o.get("ont_support", "none")
            ev = {"asm": "assembly", "ont": "ONT_reads_only", None: ""}[r["e" + e]]
            if r["x" + e]:
                ev = "not_counted_conflicting"
            return [ev, o.get("reads", "") if m != "none" else "", m if m != "none" else "",
                    o.get("clip_dist_from_end", "") if m == "trim" else ""]
        w.writerow([f"Naganishia_adeliensis_Ng1_{r['n']}", r["length"], *f(r["t5"]), *f(r["t3"]),
                    *g("5"), *g("3"), r["v5"], r["v3"], cls])

n_asm = sum((r["e5"] == "asm") + (r["e3"] == "asm") for r in rows)
n_ont = sum((r["e5"] == "ont") + (r["e3"] == "ont") for r in rows)
n_t2t = sum(r["cls"] == "T2T" for r in rows)
n_t2t_asm = sum(r["e5"] == "asm" and r["e3"] == "asm" for r in rows)
n_one = sum(r["cls"] == "one_end" for r in rows)
print("not counted (conflicting):", [(r["n"], e) for r in rows for e in "53" if r["x" + e]],
      "| unverifiable:", [(r["n"], e) for r in rows for e in "53" if r["v" + e] == "unverifiable"])
print(f"scaffolds={len(rows)} ends_assembly={n_asm} ends_ont_only={n_ont} T2T={n_t2t} "
      f"(both in assembly={n_t2t_asm}) one_end={n_one}")

# ---- plot ----
plt.rcParams.update({"font.family": "Arial", "font.size": 8, "svg.fonttype": "none",
                     "axes.linewidth": 0.6})
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8984"
BAR_T2T, BAR_ONE, BAR_NONE = "#2a78d6", "#9cc0ec", "#dcdbd6"
C5, C3 = "#eb6834", "#6250d6"          # categorical slots (orange, violet): 5' C-rich, 3' G-rich

fig, ax = plt.subplots(figsize=(6.0, 6.2))
rows.sort(key=lambda r: r["n"])
for i, r in enumerate(rows):
    y = i
    L = r["length"] / 1e6
    col = {"T2T": BAR_T2T, "one_end": BAR_ONE, "none": BAR_NONE}[r["cls"]]
    ax.plot([0, L], [y, y], color=col, lw=3.2, solid_capstyle="round", zorder=1)
    for e, x, c in (("5", 0, C5), ("3", L, C3)):
        if r["e" + e] == "asm":
            ax.plot(x, y, marker="o", ms=5.5, color=c, mec="white", mew=0.8, zorder=3)
        elif r["e" + e] == "ont":
            ax.plot(x, y, marker="o", ms=4.6, mfc="white", mec=c, mew=1.3, zorder=3)

ax.set_yticks(range(len(rows)))
ax.set_yticklabels([str(r["n"]) for r in rows], fontsize=6.5, color=INK2)
for lab, r in zip(ax.get_yticklabels(), rows):
    if r["cls"] == "T2T":
        lab.set_fontweight("bold"); lab.set_color(INK)
ax.invert_yaxis()
ax.set_ylim(len(rows) - 0.4, -0.8)
ax.set_xlim(-0.05, 2.35)
ax.set_xlabel("Scaffold length (Mb)", color=INK)
ax.set_ylabel("Scaffold", color=INK)
ax.tick_params(axis="both", colors=INK2, width=0.6, length=3)
ax.tick_params(axis="y", length=0)
for s in ("top", "right", "left"):
    ax.spines[s].set_visible(False)
ax.spines["bottom"].set_color(MUTED)
ax.grid(axis="x", color="#ecebe7", lw=0.6, zorder=0)
ax.set_axisbelow(True)

handles = [
    Line2D([], [], marker="o", ls="", ms=5.5, color=C5, mec="white", label="5' telomere (CCCCTAA)$_n$"),
    Line2D([], [], marker="o", ls="", ms=5.5, color=C3, mec="white", label="3' telomere (TTAGGGG)$_n$"),
    Line2D([], [], marker="o", ls="", ms=5.5, color=INK2, mec="white", label=f"In assembly ({n_asm} ends)"),
    Line2D([], [], marker="o", ls="", ms=4.6, mfc="white", mec=INK2, mew=1.3,
           label=f"In ONT reads only ({n_ont} ends)"),
    Line2D([], [], color=BAR_T2T, lw=3.2, label=f"Telomere at both ends ({n_t2t})"),
    Line2D([], [], color=BAR_ONE, lw=3.2, label=f"Telomere at one end ({n_one})"),
    Line2D([], [], color=BAR_NONE, lw=3.2, label=f"No terminal telomere ({len(rows) - n_t2t - n_one})"),
]
leg = ax.legend(handles=handles, loc="lower right", frameon=False, fontsize=7,
                handlelength=1.6, borderaxespad=0.2, labelcolor=INK)
fig.tight_layout()
fig.savefig(out + ".svg")
fig.savefig(out + ".png", dpi=300)
print("written", out + ".svg", out + ".png")
