# 05 — Telomeres and chromosome-level organisation (Fig. 1A, S1 Table)

Run order and the commands that generate the ONT alignment: [`run_order.txt`](run_order.txt).

1. **Telomeric arrays in the assembly** — [`scripts/telo_scan.py`](scripts/telo_scan.py): terminal 10 kb of each
   nuclear scaffold; units (C₂₋₅TAA) at 5′ ends and (TTAG₂₋₅) at 3′ ends (3′ end scanned as reverse complement);
   units ≤3 bp apart contiguous, runs ≤50 bp apart merged; an end is capped when an array of ≥4 units starts
   ≤200 bp from it.
2. **ONT alignment** — `minimap2 -t 28 -ax map-ont -Y Naganishia_sp_asm.scaffolds.fa Ng1.pc.fastq | samtools sort`.
3. **Telomeric soft-clips** — [`scripts/endclips.py`](scripts/endclips.py) and
   [`scripts/ont_end_support.py`](scripts/ont_end_support.py): reads soft-clipped within 12 kb of each end whose
   overhang is ≥30% telomeric (expected orientation); defines terminal segments without ONT support ('trim' ends:
   <25% of the ONT depth of the adjacent 3 kb).
4. **Read-terminus verification (ONT criterion)** — [`scripts/audit_read_termini.py`](scripts/audit_read_termini.py):
   reads crossing a point 1 kb inside each end are classified as ending in telomeric repeat (≥40% of their outermost
   150 bp) or continuing ≥500 bp past the end. Confirmed telomere: ≥5 reads ending in telomere, ≥3 uniquely placed
   (MAPQ 60, ≥3 kb aligned), more than continuing reads, median tip ≤1 kb beyond the end.
   Null: 1 of 395,880 read ends >20 kb from scaffold ends is telomeric.
5. **Classification and figure** — [`scripts/fig1A_telomeres.py`](scripts/fig1A_telomeres.py);
   **S1 Table** — [`scripts/make_table_telomeres.py`](scripts/make_table_telomeres.py);
   overhang lengths — [`scripts/overhang_length.py`](scripts/overhang_length.py).

## Results (`results/`)
* 37 telomeric ends: 16 in the assembly, 21 supported only by ONT reads; 9 telomere-to-telomere scaffolds
  (1, 2, 4, 5, 7, 11, 15, 17, 18; 9.96 Mb, 50.7% of the nuclear genome); 19 scaffolds with one telomere (95.2% of the
  nuclear assembly in scaffolds with ≥1 telomere).
* Not counted: scaffold 24 (2x ONT depth, reads cross both ends; its arrays are interstitial) and the 5′ end of
  scaffold 21 (more continuing than telomere-terminated reads).
* Missing telomeres lie 0–464 bp beyond the scaffold ends (~3.0 kb in total, 16 ends). The deposited assembly was not
  modified.
* `S1_Table_telomeres.tsv` (per end), `Fig1A_telomeres_table.tsv` (per scaffold), `audit_read_termini.tsv`,
  `ont_end_support.tsv`, `endclips.txt`, `overhang_length.tsv`.

Relaxing the distance to the end (e.g. ≤5 kb) inflates T2T counts with interstitial subtelomeric arrays (only 4
arrays of ≥4 units occur in the 18.4 Mb scaffold interiors), so it was not used.
