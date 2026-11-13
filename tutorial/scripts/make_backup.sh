#!/usr/bin/env bash
# Instructor only: regenerate tutorial/backup/ (the pick-up-here copy of every step's output).
#   conda activate fiberseq-workshop
#   bash scripts/make_backup.sh [threads]        (run from tutorial/)
# Builds the starting BAM (GM12878.workshop.bam, also the file to upload for students), runs
# run_all.sh on it in a scratch folder, removes machine-specific paths, then copies the
# outputs into backup/. The reference (ref.fa) is not kept: students download it in step 1.
#
# Full rebuild order after any tutorial change:
#   make_backup.sh  ->  make_bucket.sh  ->  upload .bucket/ (README "Instructors")
# Only backup/01.gm12878.bam(.bai) goes into git (check_setup.sh uses it); the rest is served from the bucket.

# stop at the first error, on unset variables, and on failures inside pipes
set -euo pipefail
T=${1:-4}                                  # threads
HERE=$(cd "$(dirname "$0")/.." && pwd)     # tutorial/
WORK="$HERE/.work"   # scratch folder (not published)
rm -rf "$WORK"
mkdir -p "$WORK"
# 1. the starting BAM (scripts/make_input_bam.sh)
bash "$HERE/scripts/make_input_bam.sh" "$WORK/GM12878.workshop.bam"
# 2. every tutorial step on it (scripts/run_all.sh)
INPUT_BAM="$WORK/GM12878.workshop.bam" bash "$HERE/scripts/run_all.sh" "$WORK" "$T"
# the UCSC hub students open in step 12a
bash "$HERE/scripts/make_hub.sh" "$WORK" "$WORK/hub" "$T"
cd "$WORK"

rm -f ref.fa ref.fa.fai   # ref.fa is ~100 MB; students download it in step 1

# FiberHMM records the full path of the fiberhmm-call program in its @PG line: keep just the name
# (sed keeps "CL:fiberhmm-call" and drops the folder path before it; then the header is swapped in)
samtools view -H --no-PG 11.fiberhmm.bam | sed -E 's#CL:[^[:space:]]*/(fiberhmm-[a-z-]+)#CL:\1#' > 11.hdr
samtools reheader --no-PG 11.hdr 11.fiberhmm.bam > 11.tmp.bam && mv 11.tmp.bam 11.fiberhmm.bam
rm 11.hdr && samtools index 11.fiberhmm.bam

# refuse to publish anything that still holds a local or internal path
# BAM headers first, then the text outputs
for b in *.bam; do
  if samtools view -H --no-PG "$b" | grep -q -E '/Users/|/home/|/mmfs1/|/grid/'; then
    echo "ERROR: $b header still contains a local path" >&2; exit 1
  fi
done
if grep -l -E '/Users/|/home/|/mmfs1/|/grid/' *.txt *.tsv *.bed 11.tracks/*.bed 2>/dev/null; then
  echo "ERROR: the files above contain a local path" >&2; exit 1
fi

# 3. replace backup/ with the new outputs
rm -rf "$HERE/backup" && mkdir -p "$HERE/backup"
mv "$WORK"/* "$HERE/backup/"
rmdir "$WORK"
du -sh "$HERE/backup"
