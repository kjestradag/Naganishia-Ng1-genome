#!/usr/bin/env python3
"""Test which junction configuration the IF1SW-F1 Illumina reads support.

Ng1 configuration:       4a-[R4]-4b      and 13a-[R13]-13b
IF1SW-F1 contig 23:      4a-[R]-13b      (reciprocal 13a-4b would be expected for a true translocation)
R = ~0.8 kb repeat (94.7% identity between copies):
  Ng1_4:422,821-423,651   Ng1_13:255,412-256,212   JAAZPZ010000023.1:360,382-361,207

For a BAM of read pairs mapped to a reference, count for each repeat interval [a, b]:
  spanning   : properly oriented pairs (same contig, MAPQ >= Q both mates) whose fragment starts before a-F and
               ends after b+F (the fragment physically links both flanks across the repeat)
  flank-to-R : pairs with one mate in the left/right unique flank (<= 1 kb away) and the other mate inside R,
               MAPQ >= Q both (links each flank to this specific repeat copy)
And discordant pairs linking the 4a/4b flanks with the 13a/13b flanks (Ng1 coordinates only).
Usage: python3 junction_pairs.py <bam> ng1|if1sw [MAPQ]
"""
import sys, subprocess, collections
bam, mode = sys.argv[1], sys.argv[2]
Q = int(sys.argv[3]) if len(sys.argv) > 3 else 20
F, W = 20, 1000
P = "Naganishia_adeliensis_Ng1_"
REP = {"ng1": {"R4 (Ng1_4)": (P + "4", 422821, 423651), "R13 (Ng1_13)": (P + "13", 255412, 256212)},
       "if1sw": {"R (contig 23, 4a-13b)": ("JAAZPZ010000023.1", 360382, 361207)}}[mode]

def reads(c, a, b, extra=()):
    out = subprocess.run(["samtools", "view", "-F", "0xF0C", *extra, bam, f"{c}:{a}-{b}"],
                         capture_output=True, text=True).stdout.splitlines()
    for r in out:
        f = r.split("\t")
        yield f[0], int(f[1]), f[2], int(f[3]), int(f[4]), f[6], int(f[7]), int(f[8])

print(f"# {bam}  mode={mode}  MAPQ>={Q}")
for lab, (c, a, b) in REP.items():
    span = set(); fl = collections.Counter()
    for name, flag, ch, pos, mq, mch, mpos, tlen in reads(c, a - W - 600, b + W + 600):
        if mq < Q: continue
        # proper orientation: this read forward and leftmost, mate on same contig
        if mch == "=" and tlen > 0 and not flag & 16 and flag & 32:
            s, e = pos, pos + tlen - 1
            if s < a - F and e > b + F and tlen < 3000:
                span.add(name)
    # flank-to-R links (mate MAPQ via MQ tag)
    out = subprocess.run(["samtools", "view", "-F", "0xF0C", "-q", str(Q), bam, f"{c}:{a}-{b}"],
                         capture_output=True, text=True).stdout.splitlines()
    for r in out:
        f = r.split("\t"); pos = int(f[3]); mpos = int(f[7])
        mq_m = [t for t in f[11:] if t.startswith("MQ:i:")]
        if not (a <= pos <= b - 100) or f[6] != "=" or (mq_m and int(mq_m[0][5:]) < Q): continue
        if a - W <= mpos < a - 150: fl["left flank <-> R"] += 1
        elif b < mpos <= b + W: fl["right flank <-> R"] += 1
    print(f"{lab:24s} spanning fragments: {len(span):4d}   " + "  ".join(f"{k}: {v}" for k, v in sorted(fl.items())))

if mode == "ng1":
    # discordant pairs between the flanks of the two Ng1 junctions
    zones = {"4a": (P + "4", 421821, 422820), "4b": (P + "4", 423652, 424651),
             "13a": (P + "13", 254412, 255411), "13b": (P + "13", 256213, 257212)}
    def zone(c, p):
        for z, (zc, za, zb) in zones.items():
            if c == zc and za <= p <= zb: return z
    links = collections.Counter()
    for z, (c, a, b) in zones.items():
        for name, flag, ch, pos, mq, mch, mpos, tlen in reads(c, a, b, ("-q", str(Q))):
            mc = ch if mch == "=" else mch
            z2 = zone(mc, mpos)
            if z2 and z2 != z: links[tuple(sorted((z, z2)))] += 1
    print("pairs linking flanks (each pair counted from both mates):", dict(links))
