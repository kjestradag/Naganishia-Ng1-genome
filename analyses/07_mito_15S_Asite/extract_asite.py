#!/usr/bin/env python3
"""Print the A-site residues (E. coli numbering) and their alignment posterior probabilities from a cmalign
Stockholm file. E. coli standard numbering = position in J01859.1 + 1.
Usage: python3 extract_asite.py all.sto
"""
import sys
seqs, pp = {}, {}
for l in open(sys.argv[1]):
    if l.startswith("#=GR"):
        _, n, t, s = l.split()
        if t == "PP": pp[n] = pp.get(n, "") + s
    elif l.startswith(("#", "//")) or not l.strip():
        continue
    else:
        n, s = l.split(); seqs[n] = seqs.get(n, "") + s
e = seqs["Ecoli_16S_rrsB"]; col = {}; k = 0
for i, ch in enumerate(e):
    if ch not in ".-":
        k += 1; col[k + 1] = i                     # J01859 residue k = E. coli position k+1
pos = list(range(1404, 1413)) + list(range(1488, 1498))
print("position             " + " ".join(f"{p:>4}" for p in pos))
for n in seqs:
    print(f"{n[:20]:20s} " + " ".join(f"{seqs[n][col[p]].upper():>4}" for p in pos))
print("\nposterior probability ('*' = >=0.95)")
for n in seqs:
    print(f"{n[:20]:20s} " + " ".join(f"{pp[n][col[p]]:>4}" for p in pos))
