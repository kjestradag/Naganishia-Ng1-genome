# 01 — Hybrid genome assembly

Final assembly: `data/genome/Naganishia_adeliensis_Ng1_genome.fa.gz` (40 sequences, 19,694,392 bp).
Snakemake workflows called below: [`workflows/`](workflows) (`Snakefile_GARM2`, `Snakefile_bwa`, `Snakefile_annot`).

| Step | Tool | Command / note |
|---|---|---|
| 1 | Flye 2.9.4 | Nanopore assembly (`--nano-raw`) (`NANO_ASM`); independent Illumina assembly with SPAdes (`ILLUM_ASM`) |
| 2 | Pilon 1.24 | 7 iterations of `NANO_ASM` polishing with Illumina reads (~50x) |
| 3 | ART | error-free Illumina reads simulated from both assemblies (60x) |
| 4 | SPAdes | re-assembly of the simulated reads (`ILLSIM_ASM`) |
| 5 | GARM2 | merge of `ILLSIM_ASM` and `NANO_ASM`: `snakemake --cores 24 --config G1=scaffolds.fasta G2=Naganishia_sp_np7pilon.fasta prefix=spades_simul_and_nanop iterations=10 -s Snakefile_GARM2` |
| 6 | redundans | `redundans.py -v -f garm_consensus.fasta -o asm_redundans.fasta` |
| 7 | ART + bwa + BESST 2.2.8 | simulated mate-pair libraries (3 and 10 kb) from both assemblies: `art_illumina -i <asm> -mp -l 250 -f 60 -m 10000 -s 500 -na -ir 0 -ir2 0 -dr 0 -dr2 0 -qs 93 -qs2 93`; mapping with `Snakefile_bwa`; `runBESST -c contigs.reduced.fa -f <4 BAMs> -o besstall -orientation rf rf rf rf` |
| 8 | ULTRA 0.9.0 | `ultra -o ultra_repeatsmask.bed --mask genome_masked.fa -t 24 --min_length 5 --pval -p 200 --min_split_window 10 --at 0.57 Scaffolds_pass4.fa` |
| 9 | funannotate 1.8.7 | `funannotate sort -i genome_masked.fa -o genome_masked_sorted.fa -b Naganishia_sp --minlen 1000` |
| 10 | — | mitochondrial contig replaced by the curated mitogenome (see `../06_mitochondrial_genome`) |

