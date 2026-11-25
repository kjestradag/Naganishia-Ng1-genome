#!/usr/bin/env python3
"""Figure 1 = panel A (compuestas/Figure1_panelA.svg, from make_composites.py) + panel B, the ORIGINAL
mitochondrial map used in the manuscript (mito/prokka_mito/fix_mito/prokka_mito/mito_v2.svg), embedded
unchanged as vector graphics. Writes compuestas/Figure1.svg and a 300-dpi PNG (rsvg-convert).
Usage (from ncbi/figuras_revision/): python3 compose_figure1.py
"""
import re, subprocess

PANEL_A = "compuestas/Figure1_panelA.svg"
MITO = "/scratch2/Project_ClauMartinez_2024_04_25_12_22_58/mito/prokka_mito/fix_mito/prokka_mito/mito_v2.svg"
OUT = "compuestas/Figure1"
W_MM, H_MM = 180.0, 101.6
PT = 72 / 25.4
W, H = W_MM * PT, H_MM * PT                      # canvas in pt

def body(svg):
    svg = re.sub(r"<\?xml[^>]*\?>", "", svg)
    svg = re.sub(r"<!DOCTYPE[^>]*>", "", svg, flags=re.S)
    return svg.strip()

def place(svg, x, y, w, h, viewbox=None):
    """Replace the attributes of the root <svg> so it is drawn at (x, y, w, h) on the canvas."""
    m = re.search(r"<svg\b[^>]*>", svg, flags=re.S)
    tag = m.group(0)
    vb = viewbox or re.search(r'viewBox="([^"]+)"', tag).group(1)
    new = (f'<svg x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" viewBox="{vb}" '
           f'preserveAspectRatio="xMidYMid meet" overflow="hidden">')
    return svg[:m.start()] + new + svg[m.end():]

# panel A (matplotlib): 3.80 x 4.0 in
a = body(open(PANEL_A).read())
aw, ah = 3.80 * 72, 4.0 * 72
a = place(a, 0, 0, aw, ah)

# panel B: original map, 20.915 x 17.902 cm, no viewBox -> user units are px at 96 dpi
m = body(open(MITO).read())
m = re.sub(r'\bid="([^"]+)"', r'id="mito_\1"', m)                 # avoid id clashes
m = re.sub(r'(href=")#([^"]+)"', r'\1#mito_\2"', m)
vw, vh = 20.915001 / 2.54 * 96, 17.902 / 2.54 * 96
bx0 = aw + 6
bw = W - bx0 - 2
bh = bw * vh / vw
by = (H - bh) / 2 + 4
m = place(m, bx0, by, bw, bh, viewbox=f"0 0 {vw:.3f} {vh:.3f}")

letters = (f'<text x="{bx0 - 2:.2f}" y="{0.01 * ah + 8:.2f}" font-family="Arial" font-weight="bold" '
           f'font-size="11" fill="#0b0b0b">B</text>')
svg = (f'<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'
       f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" version="1.1" '
       f'width="{W_MM}mm" height="{H_MM}mm" viewBox="0 0 {W:.2f} {H:.2f}">\n'
       f'<rect width="100%" height="100%" fill="#ffffff"/>\n{a}\n{m}\n{letters}\n</svg>\n')
open(OUT + ".svg", "w").write(svg)
subprocess.run(["rsvg-convert", "-d", "300", "-p", "300", "-b", "white", OUT + ".svg", "-o", OUT + ".png"], check=True)
from PIL import Image                                   # rsvg does not write the dpi tag
Image.open(OUT + ".png").convert("RGB").save(OUT + ".png", dpi=(300, 300))
print("written", OUT + ".svg", OUT + ".png", f"(panel B at y={by:.1f} pt, {bw:.1f} x {bh:.1f} pt)")
