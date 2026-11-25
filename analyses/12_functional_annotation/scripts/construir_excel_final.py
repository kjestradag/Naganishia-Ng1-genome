#!/usr/bin/env python3
"""Reconstruye el catálogo desde el original:
   - migra los IDs FUN_ -> NadelNg1_ en las columnas F y H
   - conserva los IDs originales al final (trazabilidad, no anotación)
   - agrega UNA sola columna de anotación con lo que eggNOG aporta de NUEVO
"""
import csv, re, copy, pickle
from collections import Counter
import openpyxl

D    = "/scratch2/Project_ClauMartinez_2024_04_25_12_22_58/ncbi/"
SRC  = D + "Naganishia_catalogo_interaccion_plantas_06082026.xlsx"
DST  = D + "Naganishia_catalogo_interaccion_plantas_06082026_FINAL.xlsx"
PKL  = "/tmp/claude-1000/-scratch2-Project-ClauMartinez-2024-04-25-12-22-58-ncbi/fc6145ef-3427-4927-868a-08765ead4d9b/scratchpad/final_map.pkl"

fm      = pickle.load(open(PKL, "rb"))
map_tx  = fm["map_tx"]                                   # FUN_x-T1 -> (locus_tag, tx_id)
cur     = {r["locus_tag"]: r for r in
           csv.DictReader(open(D + "Naganishia_Ng1_productos_curados.tsv"), delimiter="\t")}

PAT_FUN = re.compile(r"FUN_\d+(?:-T\d+)?")

# ------------------------------------------------------------- comparación
STOP = re.compile(r"\b(putative|probable|protein|domain|containing|family|subunit|"
                  r"like|type|related|superfamily|component|chain|complex)\b")

def toks(s):
    s = STOP.sub(" ", (s or "").lower())
    s = re.sub(r"[^a-z0-9 ]", " ", s)
    return {w for w in s.split() if len(w) > 3}

def ya_esta(cand, existentes_tok, existentes_raw):
    """True si 'cand' ya está representado en la anotación previa."""
    t = toks(cand)
    if not t:
        return True
    if cand.lower().strip() in existentes_raw:
        return True
    for e in existentes_tok:
        if not e:
            continue
        inter = len(t & e)
        if inter == len(t) or inter / len(t) >= 0.6:
            return True
    return False

def limpia_lista(v, sep=","):
    return [x.strip() for x in (v or "").split(sep) if x.strip() and x.strip() != "-"]

# --- mapa Pfam: nombre corto <-> accesión (exacto, no heurístico) ---
CLANS = "/scratch/free_old/databases/funannotate/Pfam-A.clans.tsv"
acc_de_nombre, nombre_de_acc = {}, {}
with open(CLANS) as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 4 and p[0].startswith("PF"):
            acc_de_nombre[p[3]] = p[0]
            nombre_de_acc[p[0]] = p[3]
print(f"mapa Pfam cargado: {len(acc_de_nombre)} familias")

# ------------------------------------------------------------- construcción
wb = openpyxl.load_workbook(SRC)
stats = Counter()

for ws in wb.worksheets:
    if ws.title == "Tabla_general":
        for j, c in enumerate(ws[1], 1):
            if c.value == "No. genes FUN únicos":
                c.value = "No. genes únicos (NadelNg1)"
        continue

    hdrs = [c.value for c in ws[1]]
    try:
        cE = hdrs.index("Productos anotados") + 1
        cF = hdrs.index("Todos los ID FUN") + 1
        cH = hdrs.index("ID FUN con SignalP") + 1
        cD = hdrs.index("Descripción Pfam") + 1
    except ValueError:
        print(f"  !! {ws.title}: encabezados inesperados"); continue

    style = ws.cell(row=1, column=cF)
    base  = ws.max_column
    cNew, cOrigF, cOrigH = base + 1, base + 2, base + 3
    for col, txt, wd in ((cNew,   "Anotación nueva (eggNOG)", 72),
                         (cOrigF, "Todos los ID FUN (original)", 40),
                         (cOrigH, "ID FUN con SignalP (original)", 40)):
        c = ws.cell(row=1, column=col)
        c.value = txt
        c.font = copy.copy(style.font); c.fill = copy.copy(style.fill)
        c.border = copy.copy(style.border); c.alignment = copy.copy(style.alignment)
        ws.column_dimensions[c.column_letter].width = wd

    ws.cell(row=1, column=cF).value = "Todos los ID transcrito (NadelNg1)"
    ws.cell(row=1, column=cH).value = "ID transcrito con SignalP (NadelNg1)"

    for r in range(2, ws.max_row + 1):
        vF = ws.cell(row=r, column=cF).value
        vH = ws.cell(row=r, column=cH).value
        if not vF:
            continue

        # --- migración de identificadores
        ids = PAT_FUN.findall(vF)
        ws.cell(row=r, column=cOrigF).value = vF
        ws.cell(row=r, column=cF).value = ", ".join(map_tx[i][1] if i in map_tx else i
                                                    for i in ids)
        if vH:
            ws.cell(row=r, column=cOrigH).value = vH
            ws.cell(row=r, column=cH).value = ", ".join(
                map_tx[i][1] if i in map_tx else i for i in PAT_FUN.findall(vH))

        # --- anotación previa de la fila, para comparar
        prev_raw = ws.cell(row=r, column=cE).value or ""
        pfam_desc = ws.cell(row=r, column=cD).value or ""
        prev_items = [x.strip() for x in prev_raw.split(";") if x.strip()]
        prev_tok = [toks(x) for x in prev_items]
        prev_low = {x.lower() for x in prev_items}
        pfam_tok = toks(pfam_desc)
        # accesión Pfam de la fila, para descartar el mismo dominio por accesión
        m_acc = re.search(r"PF\d{5}", str(ws.cell(row=r, column=2).value or ""))
        pfam_acc = m_acc.group(0) if m_acc else None

        prods, ecs, kos, cogs, cazy, pfams = [], [], [], [], [], []
        for i in ids:
            lt = map_tx.get(i, (None, None))[0]
            row = cur.get(lt)
            if not row:
                continue
            p = row["producto_final_recomendado"]
            if p and p != "hypothetical protein" and p not in prods:
                # se compara contra los productos ya anotados Y contra la
                # descripción del Pfam: un producto que solo reformula el
                # dominio no aporta nada
                if not ya_esta(p, prev_tok + [pfam_tok], prev_low):
                    prods.append(p)
                else:
                    stats["producto_omitido_por_duplicado"] += 1
            for e in limpia_lista(row["EC"]):
                if e not in ecs: ecs.append(e)
            for k in limpia_lista(row["KEGG_ko"]):
                k = k.replace("ko:", "")
                if k not in kos: kos.append(k)
            for ch in (row["COG_category"] or ""):
                if ch.isalpha() and ch not in cogs: cogs.append(ch)
            for z in limpia_lista(row["CAZy"]):
                if z not in cazy: cazy.append(z)
            for pf in limpia_lista(row["PFAMs"]):
                acc = acc_de_nombre.get(pf)
                # 1) mismo dominio por accesión: descarte exacto
                if acc and pfam_acc and acc == pfam_acc:
                    stats["pfam_omitido_exacto"] += 1
                    continue
                # 2) sin accesión conocida: comparación por texto, como respaldo
                if not acc and toks(pf.replace("_", " ")) & pfam_tok:
                    stats["pfam_omitido_por_texto"] += 1
                    continue
                etq = f"{pf} ({acc})" if acc else pf
                if etq not in pfams: pfams.append(etq)

        def campo(etq, items, cap=8):
            if not items:
                return None
            extra = f" (+{len(items)-cap})" if len(items) > cap else ""
            return f"{etq} {', '.join(items[:cap])}{extra}"

        partes = [campo("Prod+", prods, 5), campo("EC", ecs), campo("KO", kos),
                  campo("COG", sorted(cogs), 12), campo("CAZy", cazy),
                  campo("Pfam+", pfams, 6)]
        partes = [p for p in partes if p]
        ws.cell(row=r, column=cNew).value = " · ".join(partes) if partes else None
        ws.cell(row=r, column=cNew).alignment = copy.copy(
            ws.cell(row=r, column=cF).alignment)
        if partes: stats["filas_con_aporte"] += 1
        else:      stats["filas_sin_aporte"] += 1

# ------------------------------------------------------------- hoja por gen
lt_of_tx = {v[1]: v[0] for v in map_tx.values()}
old_of_lt = {v[0]: k for k, v in map_tx.items()}
PAT_NEW = re.compile(r"NadelNg1_\d+(?:-T\d+)?")
from openpyxl.styles import Font, PatternFill, Alignment
from collections import defaultdict

info = defaultdict(lambda: {"cat": [], "pfam": set(), "sp": False})
for ws in wb.worksheets:
    if ws.title == "Tabla_general":
        continue
    for r in range(2, ws.max_row + 1):
        v = ws.cell(row=r, column=6).value
        if not v: continue
        sub = ws.cell(row=r, column=1).value; pf = ws.cell(row=r, column=2).value
        sp = set(PAT_NEW.findall(ws.cell(row=r, column=8).value or ""))
        for g in PAT_NEW.findall(v):
            d = info[g]
            tag = f"{ws.title} > {sub}"
            if tag not in d["cat"]: d["cat"].append(tag)
            d["pfam"].add(pf)
            if g in sp: d["sp"] = True

gs = wb.create_sheet("Genes_anotacion_eggNOG", 1)
HDR = ["Gen (NadelNg1)", "ID FUN original", "Pfam (catálogo)", "N categorías",
       "Categorías > subcategorías", "Producto (eggNOG curado)", "Origen del producto",
       "Requiere revisión", "COG", "EC", "KEGG ko", "GOs", "PFAMs (eggNOG)", "CAZy", "SignalP"]
W = [22, 20, 14, 11, 80, 55, 26, 12, 8, 22, 26, 45, 34, 14, 10]
gs.append(HDR)
for c in gs[1]:
    c.font = Font(bold=True, color="FFFFFF"); c.fill = PatternFill("solid", fgColor="4472C4")
    c.alignment = Alignment(vertical="center", wrap_text=True)
for i, w in enumerate(W, 1):
    gs.column_dimensions[gs.cell(row=1, column=i).column_letter].width = w
gs.freeze_panes = "A2"
for g in sorted(info):
    lt = lt_of_tx.get(g, ""); rr = cur.get(lt, {}); d = info[g]
    gs.append([g, old_of_lt.get(lt, ""), "; ".join(sorted(x for x in d["pfam"] if x)),
               len(d["cat"]), " || ".join(d["cat"]),
               rr.get("producto_final_recomendado", ""), rr.get("origen", ""),
               rr.get("requiere_revision", ""), rr.get("COG_category", ""),
               rr.get("EC", ""), rr.get("KEGG_ko", "")[:200], rr.get("GOs", "")[:200],
               rr.get("PFAMs", ""), rr.get("CAZy", ""), "sí" if d["sp"] else "no"])

wb.save(DST)
print("Guardado:", DST)
for k, v in stats.most_common():
    print(f"  {v:6}  {k}")
print(f"  {gs.max_row-1:6}  genes en la hoja por gen")
