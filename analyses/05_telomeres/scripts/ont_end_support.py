#!/usr/bin/env python3
"""Classify ONT read support for telomeres at each nuclear scaffold end (from endclips.py output).

An end is read-supported when >=MINREADS ONT reads carry a telomeric soft-clipped overhang
(C-rich at 5' ends, G-rich at 3' ends; >=30% of the 400 bp next to the clip junction):
  extend : the clip cluster lies within 200 bp of the scaffold end (telomere missing from the assembly)
  trim   : the clip cluster lies further inside AND ONT depth over the terminal segment beyond it is
           <25% of the depth over the next 3 kb inward (terminal segment unsupported by long reads);
           if several clusters qualify, the one with most reads is kept
Usage: python3 ont_end_support.py endclips.txt ont.bam scaffolds.fa.fai > ont_end_support.tsv
"""
import re, sys, subprocess
clips_txt, bam, fai = sys.argv[1:4]
MINREADS, NEAR, RATIO = 5, 200, 0.25
P = "Naganishia_adeliensis_Ng1_"
L = {l.split()[0]: int(l.split()[1]) for l in open(fai)}

def depth(c, a, b):
    out = subprocess.run(["samtools", "depth", "-a", "-r", f"{c}:{a}-{b}", bam],
                         capture_output=True, text=True).stdout.splitlines()
    d = [int(x.split()[2]) for x in out]
    return sum(d) / max(len(d), 1)

print("scaffold\tend\tont_support\treads\tclip_dist_from_end\tterminal_depth\tinner_depth")
for line in open(clips_txt):
    m = re.match(r"\s*(\d+) (\d)'\s+depth_end=\s*(\d+)\s+telomeric_clip_clusters\(dist_from_end:reads\)=(\[.*?\])", line)
    n, e, clusters = m.group(1), m.group(2), eval(m.group(4))
    c = P + n; Ln = L[c]
    near = [x for x in clusters if x[0] <= NEAR and x[1] >= MINREADS]
    trim = []
    for p, nr in sorted(clusters):
        if p > NEAR and nr >= MINREADS:
            term = (1, p) if e == "5" else (Ln - p, Ln)
            inner = (p + 1, p + 3000) if e == "5" else (Ln - p - 3000, Ln - p - 1)
            dt, di = depth(c, *term), depth(c, *inner)
            if dt < RATIO * di:
                trim.append((p, nr, dt, di))
    if trim:
        p, nr, dt, di = max(trim, key=lambda t: (t[1], -t[0]))   # best-supported junction
        print(f"{c}\t{e}\ttrim\t{nr}\t{p}\t{dt:.1f}\t{di:.1f}")
    elif near:
        print(f"{c}\t{e}\textend\t{near[0][1]}\t{near[0][0]}\t\t")
    else:
        print(f"{c}\t{e}\tnone\t\t\t\t")
