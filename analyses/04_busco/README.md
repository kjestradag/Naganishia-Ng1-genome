# 04 — Genome completeness (BUSCO)

BUSCO v5.8.0, lineage basidiomycota_odb10 (1,764 genes), offline lineage, Docker image `ezlabgva/busco:v5.8.0_cv1`.

## Genome mode with miniprot, Ng1 and the 17 public Naganishia assemblies (Table 2)
```
for f in ng1.fa GCA_*_genomic.fna; do
  docker run --rm -u $(id -u):$(id -g) -v $PWD:/busco_wd -w /busco_wd ezlabgva/busco:v5.8.0_cv1 \
    busco -i $f -m genome --lineage busco_downloads/lineages/basidiomycota_odb10 --offline --cpu 28 -f --out metaeuk_${f%%.f*}
done
python3 busco_genome_mode_summary.py <busco_folder> > busco_miniprot_vs_augustus.tsv
```
BUSCO 5.8 uses miniprot by default for eukaryotic genomes (`euk_genome_min`); output folders are named `metaeuk_*`
only because of the command above. Summaries: [`all_genomes_miniprot/`](all_genomes_miniprot);
table: [`busco_miniprot_vs_augustus.tsv`](busco_miniprot_vs_augustus.tsv) (column `C` is the value used in Table 2).

**Ng1: C 92.2% [S 92.1%, D 0.2%], F 0.9%, M 6.9%** — 15 of the 18 assemblies fall between 92.1% and 93.9%.

## Annotated Ng1 proteome (Table 1, Table 2 footnote)
Protein mode on the final NCBI proteome (6,818 proteins):
```
docker run --rm -u $(id -u):$(id -g) -v $PWD:/busco_wd -w /busco_wd ezlabgva/busco:v5.8.0_cv1 \
  busco -i prot/ng1_final_proteins.faa -m proteins --lineage busco_downloads/lineages/basidiomycota_odb10 --offline --cpu 28 -f --out prot_ng1_final_proteins
```
**C 93.4% [S 93.2%, D 0.2%], F 1.9%, M 4.8%** —
[`annotated_proteomes/ng1_final_proteins_short_summary.txt`](annotated_proteomes/ng1_final_proteins_short_summary.txt).
