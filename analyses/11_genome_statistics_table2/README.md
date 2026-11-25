# 11 — Assembly statistics of Naganishia genomes (Table 2)

* [`stats_17refNaganishiaNCBI.txt`](stats_17refNaganishiaNCBI.txt): size, longest scaffold, N50/L50, N90/L90,
  number of scaffolds, GC and CDS (NCBI annotation) for the 17 public genomes and Ng1. BUSCO values in Table 2 come
  from the homogeneous miniprot runs of `../04_busco` (the BUSCO column of this file is not used).
* Ng1 row updated to the deposited assembly and annotation (19.69 Mb; 6,818 CDS = 6,774 nuclear + 44 mitochondrial;
  BUSCO 92.2% genome mode).
* [`make_table2_reduced.js`](make_table2_reduced.js) builds `figures/tables/Table2_comparative_reduced.docx`
  (`npm install docx`; `node make_table2_reduced.js Table2_comparative_reduced.docx`).
