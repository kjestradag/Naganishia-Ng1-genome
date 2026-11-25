"""teloclip-like check: for each nuclear scaffold end, find ONT reads soft-clipped within W bp of
the end (clip pointing outward) and score whether the clipped overhang is telomeric.
Usage: python3 endclips.py ont.bam scaffolds.fa"""
import sys, re, subprocess, collections
bam, fa = sys.argv[1], sys.argv[2]
W, MINCLIP = 12000, 80
TELC = re.compile(r'C{2,6}TA{1,2}'); TELG = re.compile(r'T{1,2}AG{2,6}')
lens = {}
for l in open(fa + ".fai") if False else []: pass
name = None
for line in open(fa):
    if line.startswith(">"):
        name = line[1:].split()[0]; lens[name] = 0
    else:
        lens[name] += len(line.strip())
def telfrac(s, R):
    return sum(len(m.group()) for m in R.finditer(s)) / max(len(s), 1)
cig = re.compile(r'(\d+)([MIDNSHP=X])')
res = collections.defaultdict(lambda: collections.Counter())
for ctg, L in lens.items():
    if ctg.endswith("_27"):
        continue
    for end, reg in (("5", f"{ctg}:1-{W}"), ("3", f"{ctg}:{max(1, L - W)}-{L}")):
        p = subprocess.run(["samtools", "view", "-F", "0x104", "-q", "1", bam, reg],
                           capture_output=True, text=True)
        for r in p.stdout.splitlines():
            f = r.split("\t"); pos = int(f[3]); ops = cig.findall(f[5]); seq = f[9]
            reflen = sum(int(n) for n, o in ops if o in "MDN=X")
            if end == "5" and ops[0][1] == "S" and int(ops[0][0]) >= MINCLIP and pos <= W:
                clip = seq[:int(ops[0][0])][-400:]           # 400 bp next to the junction, C-rich expected
                key = (pos - 1) // 100 * 100                 # bin of clip position (from 5' end)
                res[(ctg, end)][(key, telfrac(clip, TELC) >= 0.3)] += 1
            if end == "3" and ops[-1][1] == "S" and int(ops[-1][0]) >= MINCLIP:
                e = pos + reflen - 1
                if e >= L - W:
                    clip = seq[-int(ops[-1][0]):][:400]
                    key = (L - e) // 100 * 100               # distance of clip from 3' end
                    res[(ctg, end)][(key, telfrac(clip, TELG) >= 0.3)] += 1
        # depth right at the terminal 200 bp
        a, b = (1, 200) if end == "5" else (L - 199, L)
        d = subprocess.run(["samtools", "depth", "-a", "-r", f"{ctg}:{a}-{b}", bam],
                           capture_output=True, text=True).stdout.split("\n")
        dv = [int(x.split("\t")[2]) for x in d if x]
        res[(ctg, end)][("depth", None)] = round(sum(dv) / max(len(dv), 1))
for (ctg, end), c in sorted(res.items(), key=lambda kv: (int(kv[0][0].split("_")[-1]), kv[0][1])):
    tel = sorted(((k, n) for (k, t), n in c.items() if t is True and n >= 3), key=lambda x: -x[1])
    non = sum(n for (k, t), n in c.items() if t is False)
    print(f"{ctg.split('_')[-1]:>3} {end}'  depth_end={c[('depth', None)]:>3}  "
          f"telomeric_clip_clusters(dist_from_end:reads)={tel[:4]}  nontel_clips={non}")
