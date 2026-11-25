# 08 — Phylogenomics (Fig. 2B)

* **UFCG** (fungal universal core genes): 9 genomes, protein mode, 50% filter → 58 genes, concatenated alignment of
  43,635 aa ([`tree.log`](tree.log), [`aligned_concatenated.fasta.gz`](aligned_concatenated.fasta.gz),
  original UFCG tree [`concatenated.nwk`](concatenated.nwk); run of 2025-03-18).
* **IQ-TREE 3.1.1** re-inference on the same alignment (same topology as the original UFCG tree):
  ```
  iqtree -s aligned_concatenated.fasta -m MFP -B 1000 -alrt 1000 -T 8 -o GCA_000002715.1_ASM271v1_genomic --prefix ng1
  ```
  Best model (BIC): VT+F+I+R3. Node labels in [`ng1.treefile`](ng1.treefile) are SH-aLRT/UFBoot; the figure shows
  UFBoot. Full report: [`ng1.iqtree`](ng1.iqtree).
* Figure: [`fig2B_tree.py`](fig2B_tree.py) (`python3 fig2B_tree.py ng1.treefile Fig2B_UFCG_tree`); composite in
  `figures/scripts/make_composites.py`.
* Genomes: *A. clavatus* NRRL 1 GCA_000002715.1 (outgroup); *N. adeliensis* IF1SW-F1 GCA_012922715.1;
  *Naganishia* sp. IF7SW-B1 GCA_012922685.1; *N. liquefaciens* N6 GCA_013423385.1; *N. diffluens* bongo
  GCA_037832395.1; *N. kalamii* FJI-L2-BK-P3 GCA_022813745.1; *N. albida* JCM 2334 GCA_001599735.1;
  *N. vishniacii* ANT03-052 GCA_015708705.1; *N. adeliensis* Ng1 (this study).
