A-site (helix 44 decoding site) analysis of the N. adeliensis Ng1 mitochondrial 15S rRNA
=====================================================================================
Purpose: test whether the Ng1 mitoribosome carries a bacterial-type aminoglycoside-binding site (S1 Fig).

Software: Infernal 1.1.5 (cmsearch, cmalign), samtools 1.19.2, curl, Python 3 + matplotlib.
Input sequences, the structural alignment and the scripts are in this folder.

---------------------------------------------------------------------------
STEP 0. Working folder
---------------------------------------------------------------------------
mkdir -p ~/asite_repro && cd ~/asite_repro
P=/scratch2/Project_ClauMartinez_2024_04_25_12_22_58

---------------------------------------------------------------------------
STEP 1. Sequences
---------------------------------------------------------------------------
# 1a. Ng1 mitochondrial 15S (rrnS). Coordinates from the MITOS2 annotation of the final mitogenome
#     (mito_v2.fna, 54,731 bp, identical to deposited contig 27):
#     $P/mito/prokka_mito/fix_mito/prokka_mito/Galaxy3_MITOS2_on_dataset_1_GFF.gff  (line "rRNA ... rrnS")
samtools faidx $P/mito/prokka_mito/fix_mito/prokka_mito/mito_v2.fna Naganishia_sp_27:53177-54523 \
  | sed 's/>.*/>Ng1_rns/' > ng1_rns.fa                                  # 1,347 bp

# 1b. S. cerevisiae mitochondrial 15S (NC_027264.1:6523-8170)
sed 's/^>.*/>Scer_mito_15S/' $P/ref_tree/scerev/srna.fna > scer_15S.fa

# 1c. Controls from NCBI
E="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta"
curl -s "$E&id=J01859.1"   | sed 's/>.*/>Ecoli_16S_rrsB/' > ecoli_16S.fa   # E. coli 16S, 1,541 bp
curl -s "$E&id=NR_003286.4" | sed 's/>.*/>Human_cyto_18S/' > human_18S.fa  # human cytosolic 18S, 1,869 bp
curl -s "$E&id=NC_012920.1&seq_start=648&seq_stop=1601" | sed 's/>.*/>Human_mito_12S/' > human_12S.fa  # human mito 12S
curl -s "$E&id=NC_018792.1" > crypto_mito.fa                               # C. neoformans H99 mitogenome

# 1d. Rfam covariance model: bacterial SSU rRNA
curl -s -o RF00177.cm "https://rfam.org/family/RF00177/cm"

# 1e. C. neoformans mitochondrial 15S, located with the SSU model
cmsearch --cpu 8 --tblout crypto.tbl RF00177.cm crypto_mito.fa > /dev/null
grep -v '^#' crypto.tbl | head -1        # best hit: NC_018792.1 16994-18302, strand +, E = 4e-241
samtools faidx crypto_mito.fa
samtools faidx crypto_mito.fa NC_018792.1:16994-18302 | sed 's/>.*/>Cneoformans_mito_rns/' > crypto_rns.fa
#   (if a hit is on the - strand use: samtools faidx -i crypto_mito.fa CHR:END-START)

# Optional check: Ng1 15S covers the full model (model 1-1533, sequence 1-1346, E = 6e-229)
cmsearch --tblout ng1.tbl RF00177.cm ng1_rns.fa > /dev/null

---------------------------------------------------------------------------
STEP 2. Structural alignment
---------------------------------------------------------------------------
cat ecoli_16S.fa human_18S.fa human_12S.fa scer_15S.fa crypto_rns.fa ng1_rns.fa > all.fa
cmalign --cpu 8 -o all.sto RF00177.cm all.fa
#   all.sto (Stockholm) contains the consensus structure (#=GC SS_cons) and per-residue posterior probabilities (#=GR PP)

---------------------------------------------------------------------------
STEP 3. Extract the A site (E. coli numbering)
---------------------------------------------------------------------------
python3 $P/ncbi/figuras_revision/asite/extract_asite.py all.sto
#   Prints positions 1404-1412 and 1488-1497 and their posterior probabilities ('*' = >=0.95).
#   Numbering: E. coli standard position = position in J01859.1 + 1. Verified with the conserved A-site motif
#   (G1405-U1406-C1407-A1408 / U1490-G1491-A1492-A1493-G1494) and with the controls:
#   E. coli = A1408, C1409-G1491; human cytosolic 18S = G1408, C1409.A1491 (textbook configurations).

---------------------------------------------------------------------------
STEP 4. Figure
---------------------------------------------------------------------------
python3 $P/ncbi/figuras_revision/asite/fig_asite.py all.sto FigS1_Asite     # -> FigS1_Asite.svg / .png (300 dpi)

---------------------------------------------------------------------------
RESULT (E. coli numbering; posterior probability >=0.95 at all A-site columns)
---------------------------------------------------------------------------
                                        1408   1409-1491
  E. coli 16S (bacterial, sensitive)    A      C-G
  human 18S (cytosolic, resistant)      G      C.A
  S. cerevisiae mito 15S                A      C-G
  C. neoformans mito 15S                A      C-G
  N. adeliensis Ng1 mito 15S            A      U-A   (as human mitochondrial variant m.1494C>T)
  In Ng1 all pairs of this helix region are Watson-Crick or G.U.

Human mito 12S is in the alignment but EXCLUDED from the figure: cmalign places it one position off in this region
(literature: m.1494 = E. coli 1409, m.1555 = 1491; cmalign assigns 1493/1556) because human 12S has two extra
nucleotides (C1493, C1556). STEP 3 shows it with C at 1491; do not use that row.

Interpretation: Ng1 (like S. cerevisiae and C. neoformans) keeps a bacterial-type A site, predicting high-affinity
aminoglycoside binding to the mitoribosome; this is a compatible mechanism for the gentamicin/kanamycin sensitivity,
not proof of it. References: Recht et al. 1999 EMBO J (A1408); Prezant et al. 1993 Nat Genet (m.1555A>G);
Zhao et al. 2004 Am J Hum Genet (m.1494C>T).

Files in this folder: ng1_rns.fa, scer_15S.fa, ecoli_16S.fa, human_18S.fa, human_12S.fa, crypto_rns.fa, all.fa,
all.sto, extract_asite.py, fig_asite.py. RF00177.cm is downloaded in STEP 1d; the figure is in figures/supplementary/.

