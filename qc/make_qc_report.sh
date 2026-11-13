#!/usr/bin/env bash
# Instructor only: rebuild the GM12878 fiberseq-qc report (qc/GM12878_report/, open GM12878.overview.html).
#   conda env create -f qc/qc-environment.yml && conda activate fiberseq-qc
#   bash qc/make_qc_report.sh <work dir>           (needs ~4 GB disk for the hg38 autosomes)
# Sample: reads in 80 random 100 kb autosomal windows (windows.bed) from the public GM12878 FIRE CRAM,
# old feature tags removed, nucleosomes called with ft add-nucleosomes (as in the tutorial), then
# github.com/fiberseq/fiberseq-qc (commit 4b0aa3b). It reports on this sample, not the whole genome.
#
# Steps: 1. hg38 autosomes (to read the CRAM)  2. sample reads from the windows  3. call nucleosomes
#        4. run fiberseq-qc  5. check for failed plots / local paths  6. copy the report into qc/

# stop at the first error, on unset variables, and on failures inside pipes
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)        # qc/
W=${1:?usage: make_qc_report.sh <work dir>}
mkdir -p "$W" && cd "$W"
CRAM=https://s3.kopah.uw.edu/stergachis/public/FIRE/broadly-consented/GM12878/GM12878-fire-v0.1-filtered.cram

# 1. reference; skipped if already downloaded into the work dir
if [ ! -s ref.fa.fai ]; then   # hg38 autosomes, to decode the CRAM
  for i in $(seq 1 22); do curl -sSf https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr$i.fa.gz | gunzip; done > ref.fa
  samtools faidx ref.fa
fi
# 2. windows.bed (0-based BED) -> samtools regions "chr:start+1-end"; pull those reads, drop old feature tags
REGIONS=$(grep -v '^#' "$HERE/windows.bed" | awk '{printf "%s:%d-%d ", $1, $2+1, $3}')
samtools view -T ref.fa -b --no-PG -x ns,nl,as,al,aq -o raw.bam "$CRAM" $REGIONS
# reads from many regions come out unsorted: sort + index
samtools sort -o raw.sorted.bam raw.bam && samtools index raw.sorted.bam
# 3. nucleosomes + MSPs, as in tutorial step 4 (the QC plots need them)
ft add-nucleosomes -t 8 raw.sorted.bam GM12878.qc-sample.bam
samtools index GM12878.qc-sample.bam
# clean up intermediates and the cached remote CRAM index
rm -f raw.bam raw.sorted.bam* GM12878-fire-v0.1-filtered.cram.crai

# 4. fiberseq-qc at a fixed commit, so the report is reproducible
[ -d fiberseq-qc ] || git clone -q https://github.com/fiberseq/fiberseq-qc.git
git -C fiberseq-qc checkout -q 4b0aa3b
mkdir -p tmp
# the tool insists on an environment named fiberseq-qc
# (runall-qc.tcsh <output dir> <sample name> <bam>; writes report/; all messages go to qc.log)
CONDA_DEFAULT_ENV=fiberseq-qc TMPDIR=$PWD/tmp tcsh fiberseq-qc/src/runall-qc.tcsh report GM12878 GM12878.qc-sample.bam > qc.log 2>&1
# 5. the tool does not always exit with an error when a plot fails: look for R / ImageMagick errors in the log
if grep -q -E "Execution halted|unable to open image" qc.log; then echo "ERROR: some plots failed, see $W/qc.log" >&2; exit 1; fi
# and refuse to publish a report that mentions a local path
if grep -l -E '/Users/|/home/|/mmfs1/' report/* > /dev/null 2>&1; then echo "ERROR: local path in report" >&2; exit 1; fi
# 6. replace qc/GM12878_report/ with the new report
rm -rf "$HERE/GM12878_report" && cp -R report "$HERE/GM12878_report"
echo "wrote $HERE/GM12878_report (open GM12878.overview.html)"
