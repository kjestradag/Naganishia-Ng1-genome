# *Naganishia adeliensis* Ng1 — chromosome-level genome assembly and annotation

Bioinformatic analyses supporting **Vázquez-Tapia *et al.*** (in preparation), *"Genome of the psychrotolerant,
plant growth-promoting yeast* Naganishia adeliensis *Ng1"*.

This repository contains the final genome and annotation submitted to NCBI, the scripts and lightweight results of
every bioinformatic analysis reported in the article, and the figures and tables generated from them. Raw reads,
alignments (BAM) and other large intermediate files are not included (see [Data availability](#data-availability)).

| | |
|---|---|
| Organism | *Naganishia adeliensis* strain Ng1 (Basidiomycota, Tremellomycetes, Filobasidiales) |
| BioProject | [PRJNA1274400](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1274400) |
| Locus tag prefix | `NadelNg1` |
| Assembly | 19.69 Mb; 40 sequences (39 nuclear scaffolds + circular mitogenome of 54,731 bp); N50 1.16 Mb (L50 7); GC 53.4% |
| Telomeres | 37 telomeric scaffold ends; 9 telomere-to-telomere scaffolds (9.96 Mb, 50.7% of the nuclear genome) |
| Annotation | 6,999 genes: 6,774 nuclear protein-coding + 154 tRNA; mitochondrion 44 CDS + 27 tRNA |
| BUSCO (basidiomycota_odb10) | 92.2% genome mode (miniprot); 93.4% annotated proteome |

---

## Repository layout

```
data/                      final genome and annotation (NCBI submission versions, gzip) + md5sums.txt
  genome/                  Naganishia_adeliensis_Ng1_genome.fa.gz
  annotation/              GFF3, proteins (.faa), CDS (.ffn), GenBank flat file (.gbf), FUN_ -> NadelNg1_ locus map
  functional_annotation/   eggNOG-mapper, curated products, SwissProt/Pfam table, plant-interaction gene catalog
analyses/                  one folder per analysis: README.md + scripts + small results
  01_assembly/                     hybrid assembly pipeline (Flye, Pilon, SPAdes, GARM2, redundans, BESST, ULTRA)
  02_gene_prediction_annotation/   funannotate, NCBI submission (table2asn), locus renaming
  03_rnaseq_support/               RNA-seq support of the gene models
  04_busco/                        BUSCO genome mode (miniprot) for Ng1 and 17 Naganishia genomes; Ng1 proteome
  05_telomeres/                    telomere detection in the assembly + ONT read evidence (Fig. 1A, S1 Table)
  06_mitochondrial_genome/         mitogenome curation and annotation (Fig. 1B)
  07_mito_15S_Asite/               aminoglycoside-binding A site of the mitochondrial 15S rRNA (S1 Fig)
  08_phylogenomics_UFCG/           UFCG phylogenomic tree (Fig. 2B)
  09_ANI/                          average nucleotide identity, pyANI ANIm (Fig. 2A)
  10_comparison_IF1SW-F1/          dot plot vs N. adeliensis IF1SW-F1 and analysis of the discordant contig (Fig. 2C)
  11_genome_statistics_table2/     assembly statistics of Naganishia genomes (Table 2)
  12_functional_annotation/        eggNOG-mapper annotation and product curation
figures/
  main/                    Figure 1, Figure 2 (PNG 300 dpi + SVG)
  supplementary/           S1 Fig (A site)
  tables/                  Table 2 (docx), S1 Table telomeres (xlsx)
  scripts/                 scripts that build the composite figures
```

Every analysis folder has its own `README.md` with the exact commands, inputs and outputs. Scripts keep the absolute
paths of the original working directory (`/scratch2/Project_ClauMartinez_2024_04_25_12_22_58`); adjust them to rerun.

---

## Methods summary

### 1. Sequencing
Genomic DNA (28 °C culture) was sequenced with Illumina NextSeq 500 (paired-end) and Oxford Nanopore MinION
(R10.4.1 flow cell FLO-MIN114, Kit 14 chemistry, barcodes 52 and 53; live basecalling in MinKNOW with the fast model dna_r10.4.1_e8.2_400bps_fast@v4.3.0). RNA from cultures grown at 28 °C and 4 °C was sequenced
with Illumina (paired-end) and used as gene-prediction evidence.

### 2. Genome assembly — [`analyses/01_assembly`](analyses/01_assembly)
Nanopore reads were assembled with Flye v2.9.4 (--nano-raw) and polished with Pilon v1.24 (7 iterations, Illumina reads).
An independent Illumina assembly (SPAdes) was combined with the Nanopore assembly through simulated error-free reads
(re-assembled with SPAdes) and merged with GARM2; redundancy was removed with redundans, and the result was scaffolded
with BESST v2.2.8 using simulated mate-pair libraries (3 and 10 kb) from both assemblies. Repeats were soft-masked with
ULTRA v0.9.0 and contigs were sorted/renamed with `funannotate sort` (min. 1 kb). The mitochondrial genome was curated
separately (section 6).

### 3. Gene prediction and annotation — [`analyses/02_gene_prediction_annotation`](analyses/02_gene_prediction_annotation)
Genes were predicted with funannotate v1.8.7 using RNA-seq alignments (HISAT2 v2.1.0) and transcript evidence obtained by
transferring and filtering proteins of *N. onofrii* DBVPG 5303 (GCA_030039215.1). The annotation was submitted to NCBI with
table2asn (locus tag `NadelNg1`); funannotate identifiers (`FUN_`) were mapped to the final locus tags
(`data/annotation/FUN_to_NadelNg1_locus_map.tsv.gz`).

### 4. RNA-seq support — [`analyses/03_rnaseq_support`](analyses/03_rnaseq_support)
20.2 M read pairs (28 °C + 4 °C) mapped with HISAT2 (96.5%). 84.9% of the nuclear genome and 98.3% of exonic bases
are covered by ≥1 read; 6,744 of the 6,774 nuclear protein-coding genes (99.6%) have RNA-seq reads.

### 5. Genome completeness — [`analyses/04_busco`](analyses/04_busco)
BUSCO v5.8.0, basidiomycota_odb10 (1,764 genes), genome mode with miniprot for Ng1 and the 17 public Naganishia
genomes of Table 2: Ng1 92.2% (15 of 18 assemblies between 92.1% and 93.9%); 93.4% on the annotated proteome.

### 6. Telomeres and chromosome-level organisation — [`analyses/05_telomeres`](analyses/05_telomeres)
Telomeric arrays — (C₂₋₅TAA)n at 5′ ends, (TTAG₂₋₅)n at 3′ ends — were searched in the terminal 10 kb of each nuclear
scaffold (≥4 units starting ≤200 bp from the end). ONT reads (minimap2 2.11, map-ont) were used to confirm telomeres:
an end was supported when ≥5 reads (≥3 uniquely placed) ended in telomeric repeat in the expected orientation within
1 kb of the end and outnumbered reads continuing past it (null: 1 of 395,880 internal read ends). Result: 37
telomeric ends (16 in the assembly, 21 supported only by ONT reads) and 9 telomere-to-telomere scaffolds.
The deposited assembly was not modified.

### 7. Mitochondrial genome — [`analyses/06_mitochondrial_genome`](analyses/06_mitochondrial_genome)
A repetitive region of the mitogenome was resolved with Nanopore reads and confirmed by primer walking; the final
circular mitogenome (54,731 bp) was annotated with MITOS2 and Prokka (genetic code 4).

### 8. Mitochondrial 15S rRNA A site — [`analyses/07_mito_15S_Asite`](analyses/07_mito_15S_Asite)
The 15S rRNA was aligned with Infernal 1.1.5 (cmalign) to the Rfam bacterial SSU rRNA model (RF00177) with
*E. coli*, human cytosolic 18S, *S. cerevisiae* and *C. neoformans* SSU rRNAs. Ng1 keeps a bacterial-type
aminoglycoside-binding site (A1408 and a U1409–A1491 pair, *E. coli* numbering).

### 9. Phylogenomics — [`analyses/08_phylogenomics_UFCG`](analyses/08_phylogenomics_UFCG)
UFCG: 58 fungal core genes (proteins, 50% filter), concatenated alignment of 43,635 aa, maximum-likelihood tree with
IQ-TREE 3.1.1 (ModelFinder: VT+F+I+R3; 1,000 ultrafast bootstrap and SH-aLRT replicates); outgroup
*Aspergillus clavatus* NRRL 1.

### 10. Average nucleotide identity — [`analyses/09_ANI`](analyses/09_ANI)
pyANI, ANIm method, eight Naganishia genomes.

### 11. Comparison with *N. adeliensis* IF1SW-F1 — [`analyses/10_comparison_IF1SW-F1`](analyses/10_comparison_IF1SW-F1)
NUCmer 3.1 alignment and dot plot. The only discordance (contig JAAZPZ010000023.1 joining Ng1 scaffolds 4 and 13)
was shown to be a misassembly of IF1SW-F1 at a 0.8-kb two-copy repeat, using Ng1 ONT reads and the IF1SW-F1
Illumina reads (SRR11658514).

### 12. Functional annotation — [`analyses/12_functional_annotation`](analyses/12_functional_annotation)
eggNOG-mapper (COG, KEGG, GO, Pfam, CAZy) on the final proteome; products were curated and merged with the
SwissProt/Pfam annotation.

---

## Software

| Step | Software (version) |
|---|---|
| Assembly | Flye 2.9.4, Pilon 1.24, SPAdes, ART (read simulation), GARM2, redundans, BESST 2.2.8, ULTRA 0.9.0 |
| Annotation | funannotate 1.8.7, HISAT2 2.1.0, samtools, table2asn |
| Mitogenome | MITOS2 (Galaxy), Prokka 1.14.6, Canu (read correction) |
| Completeness | BUSCO 5.8.0 (miniprot, HMMER 3.4) |
| Telomeres / ONT | minimap2 2.11, samtools 1.19.2, Python 3 |
| rRNA structure | Infernal 1.1.5, Rfam RF00177 |
| Phylogenomics | UFCG (MAFFT alignment, default), IQ-TREE 3.1.1 |
| ANI / comparison | pyANI (ANIm), MUMmer/NUCmer 3.1, minimap2 2.11 |
| Functional annotation | eggNOG-mapper 2.1.12, DIAMOND 2.2.5 |
| Figures | Python 3, matplotlib 3.11, Biopython, docx (npm) |

---

## Data availability

* Final genome and annotation: this repository (`data/`) and NCBI BioProject
  [PRJNA1274400](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1274400).
* Raw reads (Illumina DNA, Nanopore, RNA-seq): SRA accessions to be added.
* Public genomes used: accessions listed in Table 2 (`figures/tables/Table2_comparative_reduced.docx`) and in the analysis
  READMEs; IF1SW-F1 Illumina reads: SRA SRR11658514.

Checksums of all files in `data/` are in [`data/md5sums.txt`](data/md5sums.txt).

## Contact

Karel Estrada — bioinformatic analyses.
