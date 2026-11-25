#!/usr/bin/env python3
"""Independent audit of telomere calls, based on where ONT reads END rather than on soft-clips.

For each nuclear scaffold end, take an anchor 1,000 bp inside the end (or inside the junction for
'trim' ends) and every read (primary or supplementary, MAPQ >= 1) aligned across that anchor. Look at the
outermost 150 bp of the read in the outward direction:
  tel_term : the read terminates in telomeric repeat in the expected orientation (>=40% of 150 bp)
  continue : the read extends >=500 bp beyond the scaffold end/junction and does NOT end in telomere
             (=> the chromosome continues past this point: assembly break, not a chromosome end)
  inside   : the read ends before reaching the end/junction (uninformative)
Unique support: reads ending in telomere with MAPQ 60 and >=3 kb aligned.
Null: the same telomeric test applied to read ends whose alignment terminates >20 kb from any scaffold end.
Usage: python3 audit_read_termini.py ont.bam scaffolds.fa.fai ont_end_support.tsv
"""
import sys, re, csv, subprocess, statistics as st, collections
bam, fai, sup = sys.argv[1:4]
L = {l.split()[0]: int(l.split()[1]) for l in open(fai)}
TELC = re.compile(r'C{2,6}TA{1,2}'); TELG = re.compile(r'T{1,2}AG{2,6}')
cig = re.compile(r'(\d+)([MIDNSHP=X])')
def frac(s, R): return sum(len(m.group()) for m in R.finditer(s)) / max(len(s), 1)
TERM, THR, ANCH, CONT = 150, 0.40, 1000, 500

def parse(r):
    f = r.split("\t"); ops = [(int(n), o) for n, o in cig.findall(f[5])]
    pos = int(f[3]); reflen = sum(n for n, o in ops if o in "MDN=X")
    return f[0], int(f[1]), int(f[4]), pos, pos + reflen - 1, ops, f[9]

# ---------- null distribution: read ends terminating internally ----------
null = collections.Counter()
p = subprocess.Popen(["samtools", "view", "-F", "0x904", "-q", "1", bam], stdout=subprocess.PIPE, text=True)
for r in p.stdout:
    name, flag, mq, s, e, ops, seq = parse(r)
    c = r.split("\t", 3)[2]
    if c.endswith("_27") or len(seq) < 2 * TERM: continue
    Ln = L[c]
    if s > 20000 and (ops[0][1] != "S" or ops[0][0] < 50):                     # left read end inside
        null["left_n"] += 1; null["left_tel"] += frac(seq[:TERM], TELC) >= THR
    if e < Ln - 20000 and (ops[-1][1] != "S" or ops[-1][0] < 50):              # right read end inside
        null["right_n"] += 1; null["right_tel"] += frac(seq[-TERM:], TELG) >= THR
print(f"# NULL internal read ends: left C-rich {null['left_tel']}/{null['left_n']} "
      f"({100*null['left_tel']/null['left_n']:.3f}%), right G-rich {null['right_tel']}/{null['right_n']} "
      f"({100*null['right_tel']/null['right_n']:.3f}%)")

# ---------- per end ----------
junction = {}
for d in csv.DictReader(open(sup), delimiter="\t"):
    if d["ont_support"] == "trim":
        junction[(d["scaffold"], d["end"])] = int(d["clip_dist_from_end"])
print("scaffold\tend\tboundary_from_end\treads_at_anchor\ttel_term\tcontinue\tinside\tunique_tel_reads\t"
      "tel_reads_median_aln_kb\ttel_reads_mapq60_frac\ttel_term_ext_beyond_boundary_median\ttel_unique_mapq60_aln3kb")
for c, Ln in L.items():
    if c.endswith("_27"): continue
    for e in ("5", "3"):
        b = junction.get((c, e), 0)                      # boundary = scaffold end, or the trim junction
        bpos = b + 1 if e == "5" else Ln - b               # boundary coordinate
        apos = bpos + ANCH if e == "5" else bpos - ANCH    # anchor inside
        if apos < 1 or apos > Ln: apos = (1 + Ln) // 2
        out = subprocess.run(["samtools", "view", "-F", "0x104", "-q", "1", bam, f"{c}:{apos}-{apos}"],
                             capture_output=True, text=True).stdout.splitlines()
        tel = cont = ins = 0; names = set(); alen = []; mq60 = 0; ext = []; uq = 0
        for r in out:
            name, flag, mq, s, en, ops, seq = parse(r)
            if e == "5":
                lead = ops[0][0] if ops[0][1] == "S" else 0
                beyond = lead - (s - bpos) if s >= bpos else lead + (bpos - s)
                is_tel = frac(seq[:TERM], TELC) >= THR
            else:
                trail = ops[-1][0] if ops[-1][1] == "S" else 0
                beyond = trail - (bpos - en) if en <= bpos else trail + (en - bpos)
                is_tel = frac(seq[-TERM:], TELG) >= THR
            if is_tel:
                tel += 1; names.add(name); alen.append((en - s) / 1000); mq60 += mq >= 60; ext.append(beyond)
                uq += (mq >= 60 and en - s >= 3000)
            elif beyond >= CONT:
                cont += 1
            else:
                ins += 1
        print(f"{c.split('_')[-1]}\t{e}\t{b}\t{len(out)}\t{tel}\t{cont}\t{ins}\t{len(names)}\t"
              f"{round(st.median(alen),1) if alen else ''}\t{round(mq60/tel,2) if tel else ''}\t"
              f"{st.median(ext) if ext else ''}\t{uq}")
