Fig 2C "Check rearrangement" (Ng1 vs N. adeliensis IF1SW-F1, GCA_012922715.1) — 2026-11-09
Source alignment: ref_tree/genomes/Naganishia_sp_asm.scaffolds.fastaVSGCA_012922715.1_ASM1292271v1_genomic.fna.delta.filter.coords
(nuclear scaffolds identical to the final assembly; only the mito differs).
- 205 alignments, 19.51 Mb aligned. Only one non-collinear block >=5 kb genome-wide:
  IF1SW-F1 contig JAAZPZ010000023.1 (760,850 bp; SPAdes 3.14, NovaSeq-only assembly, 149 scaffolds):
    1-361,207       = Ng1_4:62,447-423,651 (+)
    360,382-760,850 = Ng1_13:255,412-655,830 (+)      (no N gap at the junction)
  In IF1SW-F1 the reciprocal pieces (Ng1_13 up to 255.6 kb; Ng1_4 from 423.5 kb) are contig ends
  (JAAZPZ010000013.1 ends at 158,375; JAAZPZ010000103.1 starts at 1).
- Junction lies on a ~350-bp segment present only at Ng1_4:423.3-423.65 kb and Ng1_13:255.87-256.21 kb (92-96% identity).
- Ng1 configuration supported by ONT (Ng1.pc.fastq, minimap2 map-ont): primary reads spanning the WHOLE repeat
  (Ng1_4:422,821-423,651; Ng1_13:255,412-256,212) with >=1 / 2 / 3 kb of unique sequence beyond each repeat end:
  Ng1_4 28 / 11 / 3, Ng1_13 36 / 18 / 9; all MAPQ 60. ONT depth ~60x (genome mean 61x).
  (Earlier draft figures "18 and 23 with >=2 kb" were measured from an internal point, not from the repeat ends: superseded.)
  No ONT read joins Ng1_4:~423 kb with Ng1_13:~255 kb (the only split read found is a chimera with the mito).
- Other non-collinear contig: JAAZPZ010000025.1 (70 kb) links small Ng1 scaffolds 31, 32, 35, 36, 37 (fragmentation, not a rearrangement).
Conclusion: most likely a misjoin in the short-read IF1SW-F1 assembly across a 2-copy repeat; a genuine reciprocal
translocation in IF1SW-F1 cannot be excluded without its reads (SRA, BioProject PRJNA623412).

== Test with IF1SW-F1's own Illumina reads (SRR11658514, ENA; 7.4 M pairs NovaSeq, insert 410 +- 132 bp) — 2026-11-09
Mapped with minimap2 -ax sr to Ng1 (99.0% mapped, 96.8% properly paired) and to IF1SW-F1's own assembly (94.3%, 90.9%).
Repeat at the junction: Ng1_4:422,821-423,651 / Ng1_13:255,412-256,212 (~0.8 kb, 94.7% identity) / contig 23:360,382-361,207.
1. No pair spans the whole repeat in either assembly (insert too short) -> direct linkage test not possible (junction_pairs.py).
2. Contig 23's repeat copy is a MOSAIC of the two Ng1 copies (alternating Ng1_4-like and Ng1_13-like blocks at the
   ~200 variant columns). Reads inside the repeat body aligned to each version (H4.fa, H13.fa, H23.fa = Ng1_4,
   Ng1_13 and contig-23 windows): IF1SW-F1 reads -> 309 fit a pure Ng1 copy better than the mosaic, 0 fit the
   mosaic better than both pure copies (Ng1 Illumina control: 748 vs 0). Both pure copies are present in IF1SW-F1.
3. Read pairs (MAPQ>=20 both) around the repeat linking Ng1_4 with Ng1_13: IF1SW-F1 1 pair vs ~1,550-1,600 pairs
   within each junction (Ng1 control: 1-2 vs ~3,700).
Conclusion: IF1SW-F1 reads carry the two distinct repeat copies and no evidence of the 4a-13b junction; contig 23
joins the regions through a collapsed (mosaic) consensus of the repeat -> misassembly of IF1SW-F1, not a rearrangement.
