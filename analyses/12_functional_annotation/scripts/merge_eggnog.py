#!/usr/bin/env python3
"""Integra la anotación de eggNOG-mapper al catálogo de Naganishia adeliensis Ng1.

Produce:
  1. Naganishia_Ng1_eggnog_productos.tsv  - una fila por proteína, con producto propuesto
  2. ..._NadelNg1_eggnog.xlsx             - el catálogo con columnas de producto eggNOG
"""
import csv, re, pickle, copy, sys
from collections import defaultdict, Counter
import openpyxl

D    = "/scratch2/Project_ClauMartinez_2024_04_25_12_22_58/ncbi/"
ANN  = D + "eggnog/naganishia_ng1.emapper.annotations"
SRC  = D + "Naganishia_catalogo_interaccion_plantas_06082026_NadelNg1.xlsx"
DST  = D + "Naganishia_catalogo_interaccion_plantas_06082026_NadelNg1_eggnog.xlsx"
TSV  = D + "Naganishia_Ng1_eggnog_productos.tsv"

fm       = pickle.load(open("/tmp/claude-1000/-scratch2-Project-ClauMartinez-2024-04-25-12-22-58-ncbi/fc6145ef-3427-4927-868a-08765ead4d9b/scratchpad/final_map.pkl", "rb"))
map_tx   = fm["map_tx"]                     # FUN_x-T1 -> (locus_tag, transcript_id)
products = fm["products"]                   # locus_tag -> producto actual del .sqn
old_of_lt = {v[0]: k for k, v in map_tx.items()}
tx_of_lt  = {v[0]: v[1] for v in map_tx.values()}

# ---------------------------------------------------------------- eggNOG
rows = {}
with open(ANN) as fh:
    hdr = None
    for line in fh:
        if line.startswith("##"):
            continue
        if line.startswith("#"):
            hdr = line.lstrip("#").rstrip("\n").split("\t")
            continue
        if not line.strip() or hdr is None:
            continue
        rec = dict(zip(hdr, line.rstrip("\n").split("\t")))
        rows[rec[hdr[0]]] = rec
print(f"anotaciones eggNOG leidas: {len(rows)}")
if rows:
    print("columnas:", list(next(iter(rows.values())).keys()))

def clean(v):
    return "" if v in (None, "-", "") else v.strip()

# Nombres no informativos que NCBI rechaza o que no aportan sobre 'hypothetical protein'
VAGO = re.compile(r"^(uncharacterized|unknown|hypothetical|putative uncharacterized|"
                  r"protein of unknown function|duf\d+)", re.I)

def propose(r):
    """Devuelve (producto_propuesto, fuente)."""
    desc = clean(r.get("Description"))
    name = clean(r.get("Preferred_name"))
    if desc and not VAGO.match(desc):
        if name:
            return f"{desc}", f"eggNOG Description (gen {name})"
        return desc, "eggNOG Description"
    if name:
        return f"{name} protein", "eggNOG Preferred_name"
    return "", "sin asignacion informativa"

# ---------------------------------------------------------------- TSV
COLS = ["locus_tag", "transcript_id", "old_FUN_id", "producto_actual_sqn",
        "producto_propuesto_eggNOG", "fuente", "cambia",
        "Preferred_name", "COG_category", "EC", "KEGG_ko", "GOs", "PFAMs", "CAZy",
        "eggNOG_OGs", "max_annot_lvl", "evalue", "score"]

# Solo proteinas codificantes: los locus_tag que realmente estan en el .faa
prot_lts = []
for line in open(D + "Naganishia_adeliensis_Ng1_proteins.faa"):
    if line.startswith(">"):
        prot_lts.append(line[1:].split()[0])
print(f"proteinas en el .faa: {len(prot_lts)}")

stats = Counter()
with open(TSV, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(COLS)
    for lt in sorted(prot_lts):
        r = rows.get(lt)
        cur = products.get(lt, "")
        if r is None:
            stats["sin_hit_eggnog"] += 1
            w.writerow([lt, tx_of_lt.get(lt, ""), old_of_lt.get(lt, ""), cur,
                        "", "sin hit en eggNOG", "no"] + [""] * 11)
            continue
        prop, src = propose(r)
        cambia = "si" if prop and prop != cur else "no"
        if cur == "hypothetical protein" and prop:
            stats["hypothetical_rescatado"] += 1
        elif prop and prop != cur:
            stats["producto_cambia"] += 1
        elif not prop:
            stats["eggnog_no_informativo"] += 1
        else:
            stats["igual"] += 1
        w.writerow([lt, tx_of_lt.get(lt, ""), old_of_lt.get(lt, ""), cur, prop, src, cambia,
                    clean(r.get("Preferred_name")), clean(r.get("COG_category")),
                    clean(r.get("EC")), clean(r.get("KEGG_ko")), clean(r.get("GOs")),
                    clean(r.get("PFAMs")), clean(r.get("CAZy")),
                    clean(r.get("eggNOG_OGs"))[:80], clean(r.get("max_annot_lvl")),
                    clean(r.get("evalue")), clean(r.get("score"))])
print(f"\n[1] {TSV}")
for k, v in stats.most_common():
    print(f"    {v:6}  {k}")

# ---------------------------------------------------------------- Excel
best = {}
for lt in products:
    r = rows.get(lt)
    best[lt] = propose(r)[0] if r else ""

wb = openpyxl.load_workbook(SRC)
PAT = re.compile(r"NadelNg1_\d+(?:-T\d+)?")
lt_of_tx = {v[1]: v[0] for v in map_tx.values()}

for ws in wb.worksheets:
    if ws.title == "Tabla_general":
        continue
    hdrs = [c.value for c in ws[1]]
    try:
        cF = hdrs.index("Todos los ID transcrito (NadelNg1)") + 1
    except ValueError:
        print(f"  !! {ws.title}: sin columna de IDs"); continue
    style = ws.cell(row=1, column=cF)
    cNew = ws.max_column + 1
    c = ws.cell(row=1, column=cNew)
    c.value = "Productos eggNOG"
    c.font = copy.copy(style.font); c.fill = copy.copy(style.fill)
    c.border = copy.copy(style.border); c.alignment = copy.copy(style.alignment)
    ws.column_dimensions[c.column_letter].width = 60
    for rr in range(2, ws.max_row + 1):
        v = ws.cell(row=rr, column=cF).value
        if not v:
            continue
        prods = []
        for g in PAT.findall(v):
            p = best.get(lt_of_tx.get(g, ""), "")
            if p and p not in prods:
                prods.append(p)
        ws.cell(row=rr, column=cNew).value = "; ".join(prods) if prods else None
wb.save(DST)
print(f"\n[2] {DST}")
