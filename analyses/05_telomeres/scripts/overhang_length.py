"""Length of the telomeric overhang carried by ONT reads at each read-supported scaffold end,
and genes in terminal segments not supported by ONT reads ('trim' ends).
Usage: python3 overhang_length.py ont_end_support.tsv ont.bam scaffolds.fa.fai annotation.gff3"""
import sys, re, subprocess, statistics as st
tsv, bam, fai, gff = sys.argv[1:5]
L = {l.split()[0]: int(l.split()[1]) for l in open(fai)}
TELC = re.compile(r'C{2,6}TA{1,2}'); TELG = re.compile(r'T{1,2}AG{2,6}')
cig = re.compile(r'(\d+)([MIDNSHP=X])')
def frac(s, R): return sum(len(m.group()) for m in R.finditer(s)) / max(len(s), 1)
def tract(s, R, w=50):            # s oriented from the junction outward
    n = 0
    for i in range(0, len(s) - w + 1, w):
        if frac(s[i:i + w], R) < 0.3: break
        n += w
    return n
genes = [l.split("\t") for l in open(gff) if "\tgene\t" in l]
print("scaffold\tend\tmode\treads\tmedian_telomeric_tract_bp\tmax_tract_bp\tmedian_clip_bp\tmedian_clip_tel_frac\tunsupported_bp\tgenes_in_unsupported")
for d in list(open(tsv))[1:]:
    c, e, mode, nr, dist = d.rstrip("\n").split("\t")[:5]
    if mode == "none": continue
    dist = int(dist); Ln = L[c]
    pos = dist if e == "5" else Ln - dist
    reg = f"{c}:{max(1, pos - 300)}-{min(Ln, pos + 300)}"
    R = TELC if e == "5" else TELG
    tr, cl, fr = [], [], []
    for r in subprocess.run(["samtools", "view", "-F", "0x104", "-q", "1", bam, reg], capture_output=True, text=True).stdout.splitlines():
        f = r.split("\t"); p = int(f[3]); ops = cig.findall(f[5]); seq = f[9]
        rl = sum(int(n) for n, o in ops if o in "MDN=X")
        if e == "5" and ops[0][1] == "S" and int(ops[0][0]) >= 80 and abs((p - 1) - dist) <= 150:
            s = seq[:int(ops[0][0])][::-1]
            tel = tract(s[::-1][-len(s):][::-1], R)   # outward from junction = reversed string
            # use window scan on reversed clip: windows read right-to-left; motif check on forward windows
            ws = [seq[:int(ops[0][0])][max(0, i - 50):i] for i in range(int(ops[0][0]), 0, -50)]
            n = 0
            for w in ws:
                if len(w) < 50 or frac(w, R) < 0.3: break
                n += 50
            if frac(seq[:int(ops[0][0])][-400:], R) >= 0.3:
                tr.append(n); cl.append(int(ops[0][0])); fr.append(frac(seq[:int(ops[0][0])], R))
        if e == "3" and ops[-1][1] == "S" and int(ops[-1][0]) >= 80 and abs((Ln - (p + rl - 1)) - dist) <= 150:
            s = seq[-int(ops[-1][0]):]
            if frac(s[:400], R) >= 0.3:
                tr.append(tract(s, R)); cl.append(len(s)); fr.append(frac(s, R))
    ung = ""
    if mode == "trim":
        a, b = (1, dist) if e == "5" else (Ln - dist + 1, Ln)
        ung = sum(1 for g in genes if g[0] == c and int(g[3]) <= b and int(g[4]) >= a)
    print(f"{c.split('_')[-1]}\t{e}\t{mode}\t{len(tr)}\t{st.median(tr) if tr else ''}\t{max(tr) if tr else ''}\t"
          f"{st.median(cl) if cl else ''}\t{round(st.median(fr), 2) if fr else ''}\t{dist if mode == 'trim' else ''}\t{ung}")
