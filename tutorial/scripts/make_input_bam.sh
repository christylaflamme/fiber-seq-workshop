#!/usr/bin/env bash
# Instructor only: build the tutorial's starting BAM (run once, then upload it for students).
#   conda activate fiberseq-workshop
#   bash scripts/make_input_bam.sh <out.bam>
#
# From the public GM12878 FIRE CRAM, for two regions (one BAM):
#   SNRPN  chr15:24850000-25050000  imprinting: haplotype-selective FIREs + 5mC
#   NAPA   chr19:47480000-47525000  TF footprints (FiberHMM); one known site checked: a CTCF motif
#   keeps   alignment (hg38), haplotype tags (HP, PS), m6A + 5mC calls (MM, ML), PacBio read tags
#   removes every Fiber-seq feature call: nucleosomes (ns, nl), MSPs (as, al), FIRE scores (aq)
#   removes the producer's processing history (@PG) and reference paths (UR:) from the header
# so students start from "aligned, phased reads with m6A" and make every Fiber-seq feature themselves.

# stop at the first error, on unset variables, and on failures inside pipes
set -euo pipefail
OUT=${1:?usage: make_input_bam.sh <out.bam>}
# the lab's public GM12878 FIRE results; the CRAM holds every read, genome-wide
HUB=https://s3.kopah.uw.edu/stergachis/public/FIRE/broadly-consented/GM12878
CRAM=$HUB/GM12878-fire-v0.1-filtered.cram
REGIONS="chr15:24850000-25050000 chr19:47480000-47525000"   # SNRPN, NAPA (hg38)

# scratch folder next to the output ($$ = this run's process ID, so runs don't collide)
WORK=$(dirname "$OUT")/.make_input.$$
mkdir -p "$WORK"
# a CRAM stores reads as differences from the reference, so reading it needs hg38 chr15 + chr19
for c in chr15 chr19; do
  curl -sSf https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/$c.fa.gz | gunzip
done > "$WORK/ref.fa"
samtools faidx "$WORK/ref.fa"

# pull the two regions from the remote CRAM as BAM, dropping (-x) the old feature tags
samtools view -T "$WORK/ref.fa" -b --no-PG -x ns,nl,as,al,aq -o "$WORK/raw.bam" "$CRAM" $REGIONS
# clean header: drop @PG lines (processing history) and UR: fields (paths to the producer's reference)
samtools view -H "$WORK/raw.bam" | grep -v '^@PG' | sed $'s/\tUR:[^\t]*//' > "$WORK/clean.hdr"
samtools reheader --no-PG "$WORK/clean.hdr" "$WORK/raw.bam" > "$OUT"
samtools index "$OUT"

# checks: no feature tags left, m6A + haplotype tags present, no paths in the header
# awk counts the first 3 characters ("ns:", "MM:", ...) of every tag (SAM columns 12+) across all reads
samtools view "$OUT" | awk -F'\t' '
  { for (i = 12; i <= NF; i++) { t = substr($i, 1, 3); seen[t]++ } }
  END { bad = seen["ns:"] + seen["nl:"] + seen["as:"] + seen["al:"] + seen["aq:"]
        if (bad > 0)        { print "ERROR: feature tags remain" > "/dev/stderr"; exit 1 }
        if (!seen["MM:"])   { print "ERROR: no MM tags" > "/dev/stderr"; exit 1 }
        if (!seen["HP:"])   { print "ERROR: no HP tags" > "/dev/stderr"; exit 1 } }'
# refuse if the header still mentions a local/cluster path
if samtools view -H --no-PG "$OUT" | grep -q -E '/Users/|/home/|/mmfs1/|/grid/'; then echo "ERROR: path in header" >&2; exit 1; fi

# htslib caches the remote CRAM index in the working directory
rm -rf "$WORK" GM12878-fire-v0.1-filtered.cram.crai
echo "wrote $OUT ($(samtools view -c "$OUT") reads)"
