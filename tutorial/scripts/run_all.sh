#!/usr/bin/env bash
# Runs every command-line step of the tutorial end to end (steps 1-11). Used to (re)generate the backup outputs.
#   conda activate fiberseq-workshop
#   INPUT_BAM=<URL or path of GM12878.workshop.bam> bash scripts/run_all.sh <outdir> [threads]
# The starting BAM is made once by scripts/make_input_bam.sh. Each step's commands match README.md.
# Step 12 (FiberBrowser) runs on a laptop and has no outputs.
# Every output file is named <step number>.<name>, so you can see which step made it.
# If you change a command in README.md, change it here too (and rebuild: make_backup.sh, make_bucket.sh).

# stop at the first error, on unset variables, and on failures inside pipes
set -euo pipefail
OUT=${1:?usage: run_all.sh <outdir> [threads]}
T=${2:-4}                                  # threads (default 4)
SCRIPTS=$(cd "$(dirname "$0")" && pwd)     # tutorial/scripts/ (our Python helpers)
mkdir -p "$OUT" && cd "$OUT"               # everything below is written into <outdir>

# the lab's public GM12878 FIRE results (genome-wide peaks, hap-differences), read over the internet
HUB=https://s3.kopah.uw.edu/stergachis/public/FIRE/broadly-consented/GM12878
SNRPN=chr15:24850000-25050000     # imprinting center (hg38): open on one allele, methylated on the other
NAPA=chr19:47480000-47525000      # NAPA promoter (hg38): TF footprints
INPUT=${INPUT_BAM:?set INPUT_BAM to the starting BAM (URL or local path)}

echo "== Step 1: starting BAM + reference"
# download the starting BAM (URL) or copy it (local path)
case "$INPUT" in
  http*) curl -sSf -o 01.gm12878.bam "$INPUT" ;;
  *)     cp "$INPUT" 01.gm12878.bam ;;
esac
samtools index 01.gm12878.bam
# hg38 chr15 from UCSC; only step 10 needs it (to find CpGs). Skipped if already here.
if [ ! -s ref.fa.fai ]; then
  curl -sSf https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr15.fa.gz | gunzip > ref.fa
  samtools faidx ref.fa
fi

echo "== Step 2: look at one fiber"
# saves what students print to the screen: the first read's 9 standard SAM fields, then its tags one per line
set +o pipefail   # head closes the pipe early; that's expected here
samtools view 01.gm12878.bam | head -1 | cut -f 1-9 > 02.one-fiber.fields.txt
samtools view 01.gm12878.bam | head -1 | cut -f 12- | tr '\t' '\n' > 02.one-fiber.tags.txt
set -o pipefail

echo "== Step 3: QC"
# QC summary: statistic / value / count (fiber length, m6A per fiber, phased reads, ...), across all fibers
ft qc -t $T 01.gm12878.bam 03.qc.tsv

echo "== Step 4: per molecule: nucleosomes + MSPs"
# call nucleosomes and MSPs on every fiber from its m6A pattern; written into the Ma tag
ft add-nucleosomes -t $T 01.gm12878.bam 04.nuc.bam
samtools index 04.nuc.bam
# QC again, now that nucleosome stats exist (instructor reference only)
ft qc -t $T 04.nuc.bam 04.qc.tsv

echo "== Step 5: per molecule: FIRE"
# score every MSP; MSPs that pass become FIREs (Ma fire block + Aq scores). Flags match the FIRE pipeline.
ft fire -t $T -F 260 --min-msp 10 --min-ave-msp-size 10 --skip-no-m6a 04.nuc.bam 05.fire.bam
samtools index 05.fire.bam
# the same calls as a BED (one row per MSP; FIREs are rows with score < 100)
ft fire -t $T --extract 05.fire.bam 05.fire-elements.bed

echo "== Step 6: tags -> BED / tables"
# MSPs, nucleosomes, and a one-row-per-fiber table of everything
ft extract -t $T 05.fire.bam --msp 06.msp.bed.gz --nuc 06.nuc.bed.gz --all 06.all.tsv.gz

echo "== Step 7: pileup (per locus, SNRPN)"
# per-base counts across fibers (coverage, FIREs, ...), also split by haplotype (--haps)
ft pileup -t $T --haps --fiber-coverage -r $SNRPN -o 07.pileup.bed 05.fire.bam

echo "== Step 8: genome-wide FIRE peaks for our regions"
# peaks are called genome-wide by the FIRE pipeline; tabix pulls just our two regions from the public file
tabix -h $HUB/GM12878-fire-v0.1-peaks.bed.gz $SNRPN $NAPA > 08.peaks.bed

echo "== Step 9: haplotype-selective peaks"
# our Fisher's exact test, hap1 vs hap2, at each peak summit (scripts/fisher_hap_test.py)
python "$SCRIPTS/fisher_hap_test.py" 07.pileup.bed 08.peaks.bed > 09.hap-test.tsv
# the FIRE pipeline's own hap-differences calls for SNRPN, to compare against (instructor reference)
curl -sSf $HUB/GM12878-fire-v0.1-hap-differences.bed.gz \
  | gunzip | awk -F'\t' 'NR==1 || ($1=="chr15" && $2>=24850000 && $3<=25050000)' \
  > 09.hap-differences.pipeline.bed

echo "== Step 10: DNA methylation (5mC) on the same fibers"
# per-base 5mC counts (--cpg), by haplotype; then % methylated per CpG (scripts/cpg_by_hap.py)
ft pileup -t $T --cpg --haps --fiber-coverage --per-base -r $SNRPN -o 10.cpg.pileup.bed 05.fire.bam
python "$SCRIPTS/cpg_by_hap.py" 10.cpg.pileup.bed ref.fa > 10.cpg-by-hap.tsv

echo "== Step 11: TF footprints with FiberHMM"
# FiberHMM branch: starts again from the step 1 BAM (not from ft's output)
# call nucleosomes, MSPs and TF footprints per fiber (MA tag)
fiberhmm-call -i 01.gm12878.bam -o 11.fiberhmm.bam --enzyme hia5 --seq pacbio -c $T --no-qc
samtools index 11.fiberhmm.bam
# the calls as BED12 files in 11.tracks/ (one line per fiber, one block per feature)
fiberhmm-extract -i 11.fiberhmm.bam -o 11.tracks --tf --nucleosome --msp --bed-only -c $T

# htslib caches remote index files in the working directory; remove them
rm -f GM12878-fire-v0.1-filtered.cram.crai GM12878-fire-v0.1-peaks.bed.gz.tbi
echo "done: $(pwd)"
