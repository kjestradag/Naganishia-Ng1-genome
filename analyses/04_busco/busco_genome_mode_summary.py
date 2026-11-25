#!/usr/bin/env python3
"""Collect BUSCO v5.8.0 genome-mode (MetaEuk) results for Ng1 and the Table 2 genomes, next to the earlier
AUGUSTUS genome-mode values. Reads metaeuk_*/short_summary*.txt in the BUSCO working folder.
Usage: python3 busco_metaeuk_summary.py <busco_folder> > busco_metaeuk_vs_augustus.tsv
"""
import sys, glob, re, os
folder = sys.argv[1]
# Table 2 genomes: accession -> (name, AUGUSTUS genome-mode C%, or protein mode where marked)
T2 = {"ng1": ("N. adeliensis Ng1", "90.6"),
      "GCA_022813745.1": ("N. kalamii", "92.5"), "GCA_042257065.1": ("N. uzbekistanensis", "92.1"),
      "GCA_001599735.1": ("N. albida", "91.7"), "GCA_037832395.1": ("N. diffluens", "91.2"),
      "GCA_023628075.1": ("Naganishia sp. PG3902-1", "91.0"), "GCA_013423385.1": ("N. liquefaciens", "90.7"),
      "GCA_012922715.1": ("N. adeliensis IF1SW-F1", "90.4"), "GCA_012922685.1": ("Naganishia sp. IF7SW-B1", "90.1"),
      "GCA_025758675.1": ("N. antarctica", "90.0"), "GCA_015708705.1": ("N. vishniacii ANT03-052", "89.6"),
      "GCA_030039245.1": ("N. vaughanmartiniae", "89.0"), "GCA_023628115.1": ("Naganishia sp. PT2105A-2", "88.9"),
      "GCA_030039215.1": ("N. onofrii", "88.6"), "GCA_013461525.1": ("N. randhawae", "88.6"),
      "GCA_023628055.1": ("Naganishia sp. PT2105A", "87.6"), "GCA_030039255.1": ("N. friedmannii", "83.2 (prot)"),
      "GCA_030039265.1": ("N. cerealis", "66.8 (prot)")}
rows = []
for f in glob.glob(os.path.join(folder, "metaeuk_*", "short_summary*.txt")):
    run = os.path.basename(os.path.dirname(f))
    acc = "ng1" if run == "metaeuk_ng1" else (re.search(r"GCA_\d+\.\d", run) or [None])[0]
    if acc not in T2: continue
    txt = open(f).read()
    mode = re.search(r"run in mode: (\S+)", txt).group(1)
    m = re.search(r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%", txt)
    rows.append((acc, T2[acc][0], mode, *m.groups(), T2[acc][1]))
rows.sort(key=lambda r: -float(r[3]))
print("rank\taccession\torganism\tmode\tC\tS\tD\tF\tM\tC_augustus_previous")
for i, r in enumerate(rows, 1):
    print(i, *r, sep="\t")
missing = sorted(set(T2) - {r[0] for r in rows})
if missing:
    print("# not finished / not found:", ", ".join(missing), file=sys.stderr)
