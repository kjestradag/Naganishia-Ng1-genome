#!/usr/bin/env python3
"""Extrae los CDS nucleotídicos del envío (.sqn/.gbf) a un .ffn,
con encabezados idénticos a los del .faa, y valida traduciéndolos."""
import warnings, re
warnings.filterwarnings("ignore")
from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter

D   = "/scratch2/Project_ClauMartinez_2024_04_25_12_22_58/ncbi/"
GBF = D + "output_file.gbf"
FA  = D + "Naganishia_sp_asm.scaffolds.fa"
FAA = D + "Naganishia_adeliensis_Ng1_proteins.faa"
FFN = D + "Naganishia_adeliensis_Ng1_cds.ffn"

genome = {r.id: str(r.seq).upper() for r in SeqIO.parse(FA, "fasta")}

# encabezados del .faa, para reutilizarlos tal cual
hdr_faa = {}
for line in open(FAA):
    if line.startswith(">"):
        hdr_faa[line[1:].split()[0]] = line[1:].rstrip("\n")

feats = []
for rec in SeqIO.parse(GBF, "genbank"):
    for f in rec.features:
        if f.type != "CDS":
            continue
        feats.append(dict(
            contig=rec.name,
            locus_tag=f.qualifiers.get("locus_tag", [None])[0],
            strand=f.location.strand,
            parts=[(int(p.start), int(p.end)) for p in f.location.parts],
            start=int(f.location.start) + 1, end=int(f.location.end),
            codon_start=int(f.qualifiers.get("codon_start", ["1"])[0]),
            transl_table=f.qualifiers.get("transl_table", [None])[0],
            translation=f.qualifiers.get("translation", [None])[0],
            product=f.qualifiers.get("product", ["hypothetical protein"])[0],
        ))
print(f"CDS en el .gbf: {len(feats)}")
print("transl_table declarados:", Counter(f["transl_table"] for f in feats))

def extraer(f):
    g = genome[f["contig"]]
    partes = sorted(f["parts"])
    s = "".join(g[a:b] for a, b in partes)
    if f["strand"] == -1:
        s = str(Seq(s).reverse_complement())
    return s

st = Counter()
mito_code1, mito_code4 = 0, 0
with open(FFN, "w") as out:
    for f in feats:
        nt = extraer(f)
        off = f["codon_start"] - 1
        prot_esperada = f["translation"]

        # traducción con código 1 (lo que declara el archivo)
        aa = str(Seq(nt[off:]).translate(table=1)).rstrip("*")
        ok1 = (aa == prot_esperada)
        # y con código 4, por si el contig es mitocondrial
        aa4 = str(Seq(nt[off:]).translate(table=4)).rstrip("*")
        ok4 = (aa4 == prot_esperada)

        es_mito = f["contig"].endswith("_27")
        if es_mito:
            mito_code1 += ok1; mito_code4 += ok4
        st["coincide_code1" if ok1 else ("coincide_solo_code4" if ok4 else "NO_COINCIDE")] += 1
        if not ok1 and not ok4:
            st["_ejemplo_" + f["locus_tag"]] = 1 if st.get("_n_ej", 0) < 5 else 0

        hdr = hdr_faa.get(f["locus_tag"])
        if hdr is None:
            hdr = (f'{f["locus_tag"]} location={f["contig"]}:{f["start"]}-{f["end"]}'
                   f'({"+" if f["strand"]==1 else "-"}) product={f["product"]}')
        out.write(f">{hdr} length={len(nt)}\n")
        for i in range(0, len(nt), 60):
            out.write(nt[i:i+60] + "\n")

print(f"\n[{FFN}]  {len(feats)} CDS escritos")
for k, v in st.most_common():
    if not k.startswith("_"):
        print(f"  {v:6}  {k}")
n_mito = sum(1 for f in feats if f["contig"].endswith("_27"))
print(f"\nContig 27 (mitocondria): {n_mito} CDS")
print(f"  traducen bien con código 1: {mito_code1}")
print(f"  traducen bien con código 4: {mito_code4}")
