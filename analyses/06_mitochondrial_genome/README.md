# 06 — Mitochondrial genome (Fig. 1B)

* The mitochondrial contig of the hybrid assembly (63,121 bp) contained an unresolved repeat. It was resolved with
  Nanopore reads spanning the repeat (reads re-basecalled and corrected with Canu; consensus of the corrected reads)
  and confirmed by PCR primer walking (oligonucleotides 4742, 4743 and 4744; see article Methods).
* Final mitogenome: circular, 54,731 bp (contig `Naganishia_adeliensis_Ng1_27`), 37.8% GC.
* Annotation: MITOS2 (Galaxy) — [`annotation/mito_MITOS2.gff`](annotation/mito_MITOS2.gff) (includes rRNAs *rnl*,
  *rns*); Prokka 1.14.6 in compliant mode, genetic code 4 —
  [`annotation/mito_prokka_compliant.gff`](annotation/mito_prokka_compliant.gff).
  The NCBI annotation of contig 27 (44 CDS + 27 tRNA) is in `data/annotation/`.
* Map used in Fig. 1B (Proksee): [`mito_map_original.svg`](mito_map_original.svg). Note: the map shows Prokka
  (bacterial) gene symbols; the NCBI annotation uses mitochondrial gene names.
