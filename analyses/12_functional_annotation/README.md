# 12 — Functional annotation

* **eggNOG-mapper v2.1.12** (DIAMOND 2.2.5) on the final proteome (6,818 proteins):
  ```
  emapper.py -i Naganishia_adeliensis_Ng1_proteins.faa --itype proteins -m diamond --cpu 32 --data_dir eggnog_db \
    --output_dir eggnog -o naganishia_ng1 --tax_scope auto --go_evidence non-electronic --target_orthologs all
  ```
  Output: `data/functional_annotation/Naganishia_adeliensis_Ng1.emapper.annotations.gz`
  (COG category, EC, KEGG KO, GO, Pfam, CAZy).
* **SwissProt/Pfam annotation** of the funannotate proteome (FUN_ identifiers; map to NadelNg1_ with
  `data/annotation/FUN_to_NadelNg1_locus_map.tsv.gz`):
  `data/functional_annotation/Naganishia_Ng1_SwissProt_Pfam_annotation_FUNids.tsv.gz` — 5,560 proteins
  (5,235 with SwissProt/Pfam hits + 325 with SignalP only).
* **Merged and curated products**: [`scripts/merge_eggnog.py`](scripts/merge_eggnog.py),
  [`scripts/construir_excel_final.py`](scripts/construir_excel_final.py) →
  `data/functional_annotation/Naganishia_adeliensis_Ng1_curated_products.tsv.gz`.
  Proteins with SwissProt/Pfam or eggNOG annotation: 5,909 of 6,774 nuclear proteins (87.2%).
* **Plant-interaction gene catalog** (functional categories discussed in the article):
  `data/functional_annotation/Naganishia_Ng1_plant_interaction_gene_catalog.xlsx`.
