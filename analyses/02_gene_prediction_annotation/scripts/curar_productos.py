#!/usr/bin/env python3
"""Cura los productos de eggNOG segun las reglas de nomenclatura de NCBI.

Regla conservadora:
  - Si el gen YA tenia un producto especifico en el .sqn, se conserva.
  - Solo se rellenan las 'hypothetical protein', y solo si el nombre de eggNOG
    es aceptable (o se puede reescribir de forma segura).
"""
import csv, re
from collections import Counter

D = "/scratch2/Project_ClauMartinez_2024_04_25_12_22_58/ncbi/"
IN = D + "Naganishia_Ng1_eggnog_productos.tsv"
OUT = D + "Naganishia_Ng1_productos_curados.tsv"

# --- patrones que NCBI no acepta como nombre de proteina ---
FRASE = re.compile(r"^(belongs to|mediates|catalyzes|catalyses|involved in|required for|"
                   r"may |might |component of|part of|responsible for|acts as|binds|"
                   r"converts|plays|functions|this |the |it |in the|is a|has |provides|"
                   r"participates|contributes|essential|probable component)", re.I)
DOMINIO = re.compile(r"\b(domain|motif|repeat|fold|superfamily|family)\b", re.I)
CONTIENE_OK = re.compile(r"-containing protein$|domain-containing protein$", re.I)
DUF = re.compile(r"\b(DUF\d+|UPF\d+)\b", re.I)
VAGO = re.compile(r"^(uncharacteri[sz]ed|unknown|hypothetical|protein of unknown function)", re.I)

def limpiar(p):
    p = re.sub(r"\s+", " ", p).strip().rstrip(".")
    p = re.sub(r"\s*\([^)]*\)\s*$", "", p).strip()   # quita parentesico final
    return p

def clasificar(p):
    """Devuelve (producto_final, estado)."""
    if not p:
        return "hypothetical protein", "sin_nombre"
    q = limpiar(p)
    if not q or VAGO.match(q):
        return "hypothetical protein", "sin_nombre"
    if FRASE.match(q):
        return "hypothetical protein", "descartado_es_frase"
    if len(q) > 100:
        return "hypothetical protein", "descartado_muy_largo"
    # Ya termina en 'protein': es un nombre valido para NCBI, no tocar
    # (conserva especificidad tipo 'glycoside hydrolase family 16 protein')
    if re.search(r"\bprotein$", q, re.I):
        return q, "aceptado_ya_es_protein"
    if DUF.search(q):
        m = DUF.search(q)
        return f"{m.group(1).upper()} domain-containing protein", "reescrito_DUF"
    if DOMINIO.search(q) and not CONTIENE_OK.search(q):
        es_familia = re.search(r"\b(superfamily|family)\b", q, re.I) is not None
        base = re.sub(r"\s*\b(domain|motif|repeat|fold|superfamily|family)\b.*$", "", q, flags=re.I).strip()
        base = re.sub(r"[,;:]$", "", base).strip()
        if base and len(base) < 60:
            if es_familia:
                return f"{base} family protein", "reescrito_familia"
            return f"{base} domain-containing protein", "reescrito_dominio"
        return "hypothetical protein", "descartado_dominio_sin_base"
    return q, "aceptado"

rows = list(csv.DictReader(open(IN), delimiter="\t"))
st = Counter(); est = Counter()
COLS = ["locus_tag", "transcript_id", "old_FUN_id", "producto_actual_sqn",
        "producto_eggNOG_crudo", "producto_final_recomendado", "origen",
        "requiere_revision", "COG_category", "EC", "KEGG_ko", "GOs", "PFAMs", "CAZy"]

with open(OUT, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(COLS)
    for r in rows:
        actual = r["producto_actual_sqn"]
        crudo  = r["producto_propuesto_eggNOG"]
        if actual != "hypothetical protein":
            final, origen, rev = actual, "conservado_del_sqn", "no"
            st["conservado"] += 1
        else:
            final, estado = clasificar(crudo)
            est[estado] += 1
            if final != "hypothetical protein":
                origen = f"eggNOG ({estado})"
                rev = "si" if estado.startswith("reescrito") else "no"
                st["rellenado"] += 1
            else:
                origen, rev = "sigue hypothetical", "no"
                st["sigue_hypothetical"] += 1
        w.writerow([r["locus_tag"], r["transcript_id"], r["old_FUN_id"], actual, crudo,
                    final, origen, rev, r["COG_category"], r["EC"], r["KEGG_ko"],
                    r["GOs"], r["PFAMs"], r["CAZy"]])

print(f"[{OUT}]\n")
tot = sum(st.values())
for k, v in st.most_common():
    print(f"  {v:6}  {k}")
print(f"\n  detalle de las 3057 que eran hypothetical:")
for k, v in est.most_common():
    print(f"  {v:6}  {k}")
hyp_final = st["sigue_hypothetical"]
print(f"\n  hypothetical antes: 3057   despues: {hyp_final}   "
      f"reduccion: {3057-hyp_final} ({100*(3057-hyp_final)/3057:.0f}%)")
