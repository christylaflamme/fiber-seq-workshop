#!/usr/bin/env bash
# Instructor only: package the starting BAM and every step's outputs for the public bucket.
#   bash scripts/make_bucket.sh            (run from tutorial/, after make_backup.sh)
# Writes tutorial/.bucket/ (not tracked in git):
#   GM12878.workshop.bam(.bai)   the starting BAM (step 1 downloads it)
#   chr15.fa.gz               hg38 chr15, in case UCSC's download fails in step 1
#   hub/                      the SNRPN UCSC track hub (step 12a)
#   stepNN.tar.gz             step NN's outputs; students recover with
#                             curl -sSf $BUCKET/stepNN.tar.gz | tar xz
# Then upload the whole folder (see README "Instructors").

# stop at the first error, on unset variables, and on failures inside pipes
set -euo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd)     # tutorial/
B="$HERE/backup"                           # input: make_backup.sh's outputs
OUT="$HERE/.bucket"                        # output: the folder to upload
[ -s "$B/GM12878.workshop.bam" ] || { echo "run scripts/make_backup.sh first" >&2; exit 1; }
rm -rf "$OUT" && mkdir -p "$OUT"

# the starting BAM + index
cp "$B/GM12878.workshop.bam" "$B/GM12878.workshop.bam.bai" "$OUT/"
# hg38 chr15 (UCSC), as a fallback for step 1's reference download
curl -sSf -o "$OUT/chr15.fa.gz" https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr15.fa.gz
cd "$B"
# one archive per step: every file/folder whose name starts with the step number (e.g. 05.*)
for n in 01 02 03 04 05 06 07 08 09 10 11; do
  files=$(ls -d ${n}.* 2>/dev/null || true)
  # skip steps with no outputs
  [ -n "$files" ] || continue
  # COPYFILE_DISABLE stops macOS tar from adding ._ metadata files
  COPYFILE_DISABLE=1 tar -czf "$OUT/step$n.tar.gz" $files
done
cp -r "$B/hub" "$OUT/hub"   # step 12a: students load $BUCKET/hub/hub.txt in UCSC
# list sizes + names, and the total
ls -la "$OUT" | awk 'NR>1 {print $5, $9}'
du -sh "$OUT"
