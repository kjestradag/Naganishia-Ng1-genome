# 03 — RNA-seq support of the gene models

Input: the HISAT2 alignment from `../02_gene_prediction_annotation` (20.2 M read pairs from 28 °C and 4 °C cultures;
96.53% of primary reads mapped) and the final annotation. Nuclear scaffolds of the alignment reference are identical
in length to the final assembly; the mitochondrial contig was excluded.

Coverage was computed with mosdepth (genome-wide thresholds) and samtools bedcov / bedtools (exon intervals of the
6,774 nuclear protein-coding genes).

| Metric | Value |
|---|---|
| Nuclear genome (19,639,661 bp) covered ≥1x / ≥5x / ≥10x | 84.91% / 70.35% / 65.33% |
| Exonic bases covered ≥1x / ≥5x | 98.29% / 95.68% |
| Protein-coding genes with ≥1 read / ≥10 reads | 6,744 (99.6%) / 6,647 (98.1%) |
| Genes with ≥90% of exonic bases ≥1x / ≥5x | 6,430 (94.9%) / 6,165 (91.0%) |

tRNAs show low coverage, as expected for poly(A)-selected libraries. The BAM file (several GB) is not included.
