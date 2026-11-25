#!/usr/bin/env python3
"""Supplementary S1 Table: telomeric evidence for each end of the N. adeliensis Ng1 nuclear scaffolds.

Inputs (run from ncbi/figuras_revision/):
  Fig1A_telomeres_table.tsv              classification per scaffold (fig1A_telomeres.py, with ONT arguments)
  telomeros_ont/ont_end_support.tsv      'trim' ends: terminal segment lacking ONT support (ont_end_support.py)
  telomeros_ont/audit_read_termini.tsv   reads ending in telomere / continuing past each end (audit_read_termini.py)
  telomeros_ont/endclips.txt             ONT depth over the terminal 200 bp (endclips.py)
  ../Naganishia_sp_asm_with_annot_fixidsmito.gff3
Usage: python3 telomeros_ont/make_table_telomeres.py S1_Table_telomeres
"""
import csv, re, sys, statistics as st
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, Border, Side
from openpyxl.utils import get_column_letter

out = sys.argv[1]
P = "Naganishia_adeliensis_Ng1_"
GFF = "../Naganishia_sp_asm_with_annot_fixidsmito.gff3"
sc = list(csv.DictReader(open("Fig1A_telomeres_table.tsv"), delimiter="\t"))
ont = {(d["scaffold"], d["end"]): d for d in csv.DictReader(open("telomeros_ont/ont_end_support.tsv"), delimiter="\t")}
aud = {(P + d["scaffold"], d["end"]): d for d in
       csv.DictReader((l for l in open("telomeros_ont/audit_read_termini.tsv") if not l.startswith("#")), delimiter="\t")}
null_line = [l for l in open("telomeros_ont/audit_read_termini.tsv") if l.startswith("#")][0]
null_hits = sum(int(x) for x in re.findall(r"(\d+)/\d+", null_line))
null_tot = sum(int(x) for x in re.findall(r"\d+/(\d+)", null_line))
dep = {}
for l in open("telomeros_ont/endclips.txt"):
    m = re.match(r"\s*(\d+) (\d)'\s+depth_end=\s*(\d+)", l)
    dep[(P + m.group(1), m.group(2))] = int(m.group(3))
mrna = [l.rstrip("\n").split("\t") for l in open(GFF) if "\tmRNA\t" in l]

def genes_in(c, a, b):          # genes with >=50% of their span inside [a, b]
    ids = []
    for g in mrna:
        if g[0] != c: continue
        s, e = int(g[3]), int(g[4]); ov = min(e, b) - max(s, a) + 1
        if ov > 0.5 * (e - s + 1):
            ids.append(re.search(r"ID=([^;-]+)", g[8]).group(1))
    return ids

num = lambda x, f=int: f(float(x)) if x not in ("", None) else ""
CLASS = {"T2T": "Telomere-to-telomere", "one_end": "Telomere at one end", "none": "No terminal telomere"}
STATUS = {"assembly": "Assembly", "ONT_reads_only": "ONT reads only", "": "Absent",
          "not_counted_conflicting": "Not counted (conflicting)"}
VER = {"confirmed": "Confirmed", "unverifiable": "Unverifiable (<5 informative reads)",
       "conflicting": "Conflicting", "": ""}
head = ["Scaffold", "Scaffold length (bp)", "Scaffold class", "End", "Telomere evidence",
        "Assembled array: start from end (bp)", "Assembled array: units", "Assembled array: length (bp)",
        "ONT depth at terminal 200 bp", "ONT reads ending in telomere", "…of which uniquely placed",
        "ONT reads continuing past the end", "Telomere tip, median (bp beyond end)",
        "Terminal segment without ONT support (bp)", "Annotated genes in that segment", "Read-based verification"]
rows = []
for r in sc:
    c = r["scaffold"]; L = int(r["length_bp"])
    for e in ("5", "3"):
        o = ont[(c, e)]; a = aud[(c, e)]
        ev = r[f"tel{e}_evidence"]
        trim = o["ont_support"] == "trim"
        unsup = int(o["clip_dist_from_end"]) if trim else ""
        g = ""
        if trim:
            s_, e_ = (1, unsup) if e == "5" else (L - unsup + 1, L)
            g = ", ".join(genes_in(c, s_, e_)) or "none"
        rows.append([c.replace(P, "Ng1_"), L, CLASS[r["class"]], f"{e}′", STATUS[ev],
                     num(r[f"tel{e}_start_from_end"]), num(r[f"tel{e}_units"]), num(r[f"tel{e}_bp"]),
                     dep[(c, e)], int(a["tel_term"]), int(a["tel_unique_mapq60_aln3kb"]), int(a["continue"]),
                     num(a["tel_term_ext_beyond_boundary_median"]), unsup, g,
                     VER[r[f"tel{e}_verification"]]])

# ---- numbers quoted in the text ----
ontonly = [x for x in rows if x[4] == "ONT reads only"]
n_reads = [x[9] for x in ontonly]
ext = [x for x in ontonly if x[13] == ""]                       # ONT-only ends that are not 'trim'
tips = [max(x[12], 0) for x in ext]
print(f"ends={len(rows)} assembly={sum(x[4]=='Assembly' for x in rows)} ont_only={len(ontonly)} "
      f"not_counted={sum(x[4].startswith('Not') for x in rows)} absent={sum(x[4]=='Absent' for x in rows)}")
print(f"ONT-only ends: reads ending in telomere {min(n_reads)}-{max(n_reads)} (median {st.median(n_reads)}); "
      f"uniquely placed {min(x[10] for x in ontonly)}-{max(x[10] for x in ontonly)}")
print(f"ONT-only ends beyond the scaffold end: n={len(ext)}, tip median per end {min(tips)}-{max(tips)} bp, "
      f"total {sum(tips)} bp | trim ends: n={len(ontonly)-len(ext)}, "
      f"unsupported {sum(x[13] for x in ontonly if x[13] != '')} bp, genes: "
      f"{[x[14] for x in rows if x[14] not in ('', 'none')]}")
print(f"null: {null_hits}/{null_tot} internal read ends")

with open(out + ".tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(head); w.writerows(rows)

TITLE = "S1 Table. Telomeric repeats at the ends of the 39 nuclear scaffolds of Naganishia adeliensis Ng1."
NOTES = [
    "Assembled arrays: tandem arrays of ≥4 telomeric units, (C2–5TAA)n at 5′ ends and (TTAG2–5)n at 3′ ends "
    "(units ≤3 bp apart, runs ≤50 bp apart merged), starting ≤200 bp from the scaffold end.",
    "ONT reads were mapped with minimap2 v2.11 (-ax map-ont -Y). For each end, reads aligned (MAPQ ≥1) across a point "
    "1 kb inside the end were classified as ending in telomeric repeat in the expected orientation (≥40% of their "
    f"outermost 150 bp) or as continuing ≥500 bp past the end without telomere. Only {null_hits} of {null_tot:,} read "
    "ends located >20 kb from scaffold ends met the telomeric test. Uniquely placed: MAPQ 60 and ≥3 kb aligned.",
    "A telomere was considered supported by ONT reads when ≥5 reads (≥3 uniquely placed) ended in telomere, they "
    "outnumbered continuing reads, and their median tip lay ≤1 kb beyond the end. Telomere tip: distance between the "
    "end of the reads and the scaffold end (negative values, inside the scaffold).",
    "Terminal segment without ONT support: soft-clipped telomeric overhangs start inside the scaffold, and the "
    "sequence between them and the scaffold end has <25% of the ONT depth of the adjacent 3 kb; for these ends, read "
    "counts and tip are measured from that junction. The Ng1_11 3′ segment, which contains a 4-unit array, is 99.8% "
    "identical to scaffold Ng1_39. These segments also lack uniquely mapped Illumina coverage (depth at MAPQ ≥20: "
    "0–8×, vs 4–99× over the adjacent 3 kb; minimap2 -ax sr).",
    "'ONT reads only' with a tip ≤20 bp from the end (Ng1_4 5′, Ng1_18 3′, Ng1_22 5′ junction): degenerate telomeric "
    "repeats are present in the assembly but do not meet the array criterion.",
    "Not counted: Ng1_24 has twice the genome-wide ONT depth (108× vs 61×); reads cross both ends, and at the 3′ end "
    "they continue ~2.3 kb into scaffold Ng1_34, whose 3′ end carries the telomere. Its assembled arrays are therefore "
    "interstitial. At Ng1_21 5′, reads continuing past the junction (27) outnumber those ending in telomere (20).",
    "Absent ends with reads ending in telomere >1 kb beyond the scaffold end (Ng1_8 5′, Ng1_12 3′, Ng1_20 5′, "
    "Ng1_39 5′) lack subtelomeric sequence in the assembly and were not counted; at Ng1_16 3′ none of the reads "
    "ending in telomere was uniquely placed. Ng1_10 3′ ends in a collapsed "
    "high-copy repeat (ONT depth ~290×). The assembled array at Ng1_35 3′ has no informative ONT reads.",
    "The deposited assembly was not modified. The mitochondrial contig (Ng1_27) is excluded.",
]
wb = Workbook(); ws = wb.active; ws.title = "S1 Table"
ws.append([TITLE]); ws["A1"].font = Font(bold=True, size=11)
ws.append([])
ws.append(head)
thin = Side(style="thin", color="000000")
for j in range(1, len(head) + 1):
    cl = ws.cell(row=3, column=j)
    cl.font = Font(bold=True); cl.alignment = Alignment(wrap_text=True, vertical="bottom")
    cl.border = Border(top=thin, bottom=thin)
for x in rows:
    ws.append(x)
last = ws.max_row
for j in range(1, len(head) + 1):
    ws.cell(row=last, column=j).border = Border(bottom=thin)
for i in range(4, last + 1):
    ws.cell(row=i, column=2).number_format = "#,##0"
ws.append([])
for n in NOTES:
    ws.append([n])
widths = [9, 12, 20, 5, 17, 12, 10, 11, 11, 11, 11, 11, 12, 13, 20, 18]
for j, wdt in enumerate(widths, 1):
    ws.column_dimensions[get_column_letter(j)].width = wdt
ws.row_dimensions[3].height = 58
ws.freeze_panes = "A4"
wb.save(out + ".xlsx")
print("written", out + ".xlsx", out + ".tsv")
