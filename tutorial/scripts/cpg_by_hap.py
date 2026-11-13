#!/usr/bin/env python
"""Percent CpG methylation (5mC) per CpG, for all fibers and for each haplotype.

Input: the output of `ft pileup --cpg --haps --per-base` and the reference FASTA.
ft counts, at each base, the fibers with a 5mC call there (`cpg_coverage`, calls above ft's
default probability cutoff). PacBio reports a CpG's call on the C for reads aligned to + and on
the G for reads aligned to -, so for each CpG in the reference:
    methylated fibers = cpg_coverage[C] + cpg_coverage[G]
    % methylated      = 100 * methylated fibers / fibers covering the CpG
This per-CpG summary is our own calculation on top of ft's counts, not an ft output.

Output columns: chrom, start, end (the 2-bp CpG), then for all / H1 / H2:
cov_* (fibers covering the C), meth_* (methylated fibers), pct_* (% methylated).

usage: cpg_by_hap.py pileup.bed ref.fa > cpg-by-hap.tsv
"""
import sys
import pandas as pd

# --- the per-base pileup; strip the '#' from ft's '#chrom' header ---
pile = pd.read_csv(sys.argv[1], sep="\t")
pile.columns = [c.lstrip("#") for c in pile.columns]
chrom = pile["chrom"].iloc[0]          # the pileup covers one region, so one chromosome

# --- read that chromosome's sequence from the FASTA (to find where the CpGs are) ---
seq = []
with open(sys.argv[2]) as fh:
    keep = False
    for line in fh:
        if line.startswith(">"):
            # header line: keep the sequence lines that follow only if this is our chromosome
            keep = line[1:].split()[0] == chrom
            continue
        if keep:
            seq.append(line.strip())
seq = "".join(seq).upper()             # one string; index 0 = BED position 0

# --- one output row per CpG ---
p = pile.set_index("start")            # look up pileup rows by base position
rows = []
for pos in p.index:
    # a CpG starts here only if the reference reads "CG" at pos, pos+1 and the pileup has the G too
    if seq[pos:pos + 2] != "CG" or pos + 1 not in p.index:
        continue
    a, b = p.loc[pos], p.loc[pos + 1]  # a = the C (+ strand calls), b = the G (- strand calls)
    r = dict(chrom=chrom, start=pos, end=pos + 2)
    # same calculation for all fibers ("" columns), hap1 ("_H1") and hap2 ("_H2")
    for suf, lab in [("", "all"), ("_H1", "H1"), ("_H2", "H2")]:
        cov = int(a["coverage" + suf])                                   # fibers covering the CpG
        meth = int(a["cpg_coverage" + suf]) + int(b["cpg_coverage" + suf])  # methylated: + strand on C, - strand on G
        r[f"cov_{lab}"] = cov
        r[f"meth_{lab}"] = meth
        r[f"pct_{lab}"] = round(100 * meth / cov, 1) if cov else float("nan")  # nan where no fibers
    rows.append(r)
pd.DataFrame(rows).to_csv(sys.stdout, sep="\t", index=False)
