# 02 — Gene prediction, annotation and NCBI submission

Final files: `data/annotation/` (GFF3, proteins, CDS, GenBank flat file).

## Evidence
* **RNA-seq** (28 °C and 4 °C):
  ```
  hisat2-build genome_masked_sorted.fa genome_masked_hisat
  hisat2 -x genome_masked_hisat -p 24 -1 Ng1optimo_S38_R1_001.fastq.gz,Ng1frio_S39_R1_001.fastq.gz \
         -2 Ng1optimo_S38_R2_001.fastq.gz,Ng1frio_S39_R2_001.fastq.gz -S rna_vs_genome_hisat2.sam
  samtools view -b -@ 24 -T genome_masked_sorted.fa -o rna_vs_genome_hisat2.bam rna_vs_genome_hisat2.sam
  bamtools sort -in rna_vs_genome_hisat2.bam -out rna_vs_genome_hisat2_sorted.bam
  ```
* **Transcript evidence**: proteins of *N. onofrii* DBVPG 5303 (GCA_030039215.1) transferred to the genome (PATT
  workflow), filtered (≥100 aa), annotated (`../01_assembly/workflows/Snakefile_annot`) and kept when
  e-value < 1e-30; their CDS were used as `--transcript_evidence`.

## Prediction
```
funannotate predict --cpus 24 --force --species "Naganishia_sp" --isolate "asm" \
  --transcript_evidence cds_for_funannotate.fasta --rna_bam rna_vs_genome_hisat2_sorted.bam \
  --out funannotate_out -i genome_masked_sorted.fa
```
funannotate 1.8.7 produced 6,948 genes (6,775 mRNA + 173 tRNA) on the pre-curation assembly.

## NCBI submission
* [`table2asn_command.sh`](table2asn_command.sh). The submission template (`template.sbt`) and the assembly
  structured comment (`genome.asm`) are not included.
* **Important**: the `.sqn` that was submitted was built from `Naganishia_sp_asm_with_annot_fixidsmito.gff3`
  (= `data/annotation/Naganishia_adeliensis_Ng1.gff3.gz`), not from the GFF3 named in `table2asn_command.sh`.
* Final annotation: 6,999 genes — nuclear 6,774 mRNA + 154 tRNA; mitochondrion (contig 27, circular, 54,731 bp)
  44 CDS + 27 tRNA.
* Locus renaming `FUN_` → `NadelNg1_`: `data/annotation/FUN_to_NadelNg1_locus_map.tsv.gz`. Numbering is shifted by
  +51 after FUN_006869, and mitochondrial genes have no FUN_ equivalent
  (`data/annotation/NadelNg1_genes_without_FUN_equivalent.tsv`).
* Scripts: [`scripts/curar_productos.py`](scripts/curar_productos.py) (product-name curation for NCBI),
  [`scripts/extraer_cds.py`](scripts/extraer_cds.py) (CDS/protein FASTA from the final GFF3).
