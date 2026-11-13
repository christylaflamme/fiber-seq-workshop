#!/usr/bin/env python
"""Haplotype-selective accessibility test for one or more FIRE peaks.

For each peak, take the ft pileup row at the peak summit and build the 2x2 table
    hap1: fibers with a FIRE | fibers without
    hap2: fibers with a FIRE | fibers without
then run a two-sided Fisher's exact test (as described in the FIRE docs:
fiberseq.github.io/fire/methods/haplotype-selective.html), and rank the peaks by p-value.
p_adjust_locus is a Benjamini-Hochberg correction across the peaks tested here. The FIRE
pipeline corrects across all peaks genome-wide instead, so compare our p_value with the
pipeline's p_value column; its p_adjust will differ from ours.

This is our re-implementation of the FIRE pipeline's test, run on the tutorial's own pileup.

usage: fisher_hap_test.py pileup.bed peaks.bed
  pileup.bed : output of `ft pileup --haps`
  peaks.bed  : FIRE peaks (columns: chrom, peak_start, peak_end, start, end, ...);
               the summit is taken as the midpoint of start/end
"""
import sys
import pandas as pd
from scipy.stats import fisher_exact, false_discovery_control

# --- read both tables; ft and the FIRE pipeline start the header with '#chrom', so strip the '#' ---
pileup = pd.read_csv(sys.argv[1], sep="\t")
pileup.columns = [c.lstrip("#") for c in pileup.columns]
peaks = pd.read_csv(sys.argv[2], sep="\t")
peaks.columns = [c.lstrip("#") for c in peaks.columns]

rows = []
for _, pk in peaks.iterrows():
    # the summit: midpoint of the peak's start/end columns
    summit = (int(pk["start"]) + int(pk["end"])) // 2
    # the pileup row (an interval of constant counts) that contains the summit
    hit = pileup[(pileup["chrom"] == pk["chrom"]) & (pileup["start"] <= summit) & (pileup["end"] > summit)]
    if hit.empty:
        # peak lies outside the pileup region (e.g. the NAPA peaks; step 7 piles up SNRPN only)
        continue
    r = hit.iloc[0]
    # fibers with a FIRE at the summit, and all fibers covering it, per haplotype
    h1_fire, h1_cov = int(r["fire_coverage_H1"]), int(r["coverage_H1"])
    h2_fire, h2_cov = int(r["fire_coverage_H2"]), int(r["coverage_H2"])
    if min(h1_cov, h2_cov) < 10:  # FIRE tests only peaks with >= 10 fibers per haplotype
        continue
    # 2x2 table: rows = haplotype, columns = FIRE / no FIRE
    table = [[h1_fire, h1_cov - h1_fire], [h2_fire, h2_cov - h2_fire]]
    _, p = fisher_exact(table, alternative="two-sided")
    rows.append(dict(chrom=pk["chrom"], peak_start=pk["peak_start"], peak_end=pk["peak_end"],
                     hap1_fire=h1_fire, hap1_cov=h1_cov, hap2_fire=h2_fire, hap2_cov=h2_cov,
                     hap1_frac=round(h1_fire / h1_cov, 3), hap2_frac=round(h2_fire / h2_cov, 3), p_value=p))

# multiple-testing correction across the peaks tested here only (see docstring), then most significant first
out = pd.DataFrame(rows)
out["p_adjust_locus"] = false_discovery_control(out["p_value"], method="bh")
out = out.sort_values("p_value")
out.to_csv(sys.stdout, sep="\t", index=False)
