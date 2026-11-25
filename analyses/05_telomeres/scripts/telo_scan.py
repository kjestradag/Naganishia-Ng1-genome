"""Strict telomere scan for N. adeliensis Ng1.
Telomeric unit: C-rich (C{2,5}TAA) at 5' ends, G-rich (TTAG{2,5}) at 3' ends.
Units separated by <=GAP bp are merged into one array."""
import re, sys, json
GAP_UNIT=3; GAP_ARRAY=50
UC=re.compile(r'C{2,5}TAA'); UG=re.compile(r'TTAG{2,5}')
def load(f):
    d={};n=None
    for l in open(f):
        if l[0]=='>': n=l[1:].split()[0]; d[n]=[]
        else: d[n].append(l.strip().upper())
    return {k:''.join(v) for k,v in d.items()}
def runs(s,U,gap):
    out=[];cur=None
    for m in U.finditer(s):
        a,b=m.start(),m.end()
        if cur and a-cur[1]<=gap: cur[1]=b; cur[2]+=1
        else:
            if cur: out.append(cur)
            cur=[a,b,1]
    if cur: out.append(cur)
    return out
def arrays(s,U):
    r=[x for x in runs(s,U,GAP_UNIT) if x[2]>=4]          # tandem runs >=4 units
    out=[]
    for x in r:                                            # merge runs of the same array
        if out and x[0]-out[-1][1]<=GAP_ARRAY: out[-1]=[out[-1][0],x[1],out[-1][2]+x[2]]
        else: out.append(list(x))
    return out
COMP=str.maketrans('ACGTN','TGCAN')
def revcomp(x): return x.translate(COMP)[::-1]
def scan(fasta,W=10000):
    res=[]
    for k,s in load(fasta).items():
        L=len(s)
        a5=arrays(s[:W],UC)                     # C-rich strand read 5'->3' from the left end
        a3=arrays(revcomp(s[-W:]),UC)           # G-rich 3' end == C-rich on the reverse complement
        # positions are distances (bp) from the respective scaffold end
        res.append(dict(id=k,len=L,a5=a5,a3=a3))
    return res
if __name__=='__main__':
    for r in scan(sys.argv[1]):
        f=lambda a:'; '.join(f"{x[0]}-{x[1]}:{x[2]}u" for x in a[:3])
        print(r['id'].split('_')[-1],r['len'],"| 5':",f(r['a5'])," | 3':",f(r['a3']))
