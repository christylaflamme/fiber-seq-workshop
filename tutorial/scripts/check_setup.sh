#!/usr/bin/env bash
# Check that your computer is ready for the Fiber-seq tutorial. Run it BEFORE the workshop:
#   conda activate fiberseq-workshop
#   bash scripts/check_setup.sh            (from the tutorial/ folder)
# It runs every tool on 20 real fibers in a temporary folder (~10 s) and prints PASS / FAIL for each.
# All PASS: you're set. Any FAIL: send the whole output to the instructors before the workshop.
#
# What it checks, in order:
#   1. your operating system / chip is one fibertools has a build for
#   2. the fiberseq-workshop conda environment is the one active
#   3. every tool the tutorial uses is installed
#   4. a mini tutorial on ~20 fibers: ft add-nucleosomes, ft fire, fiberhmm-call, FiberBrowser import
#   5. the internet sites the tutorial downloads from can be reached
# Nothing is written outside a temporary folder, which is deleted at the end.

# --- settings ---
HERE=$(cd "$(dirname "$0")/.." && pwd)    # the tutorial/ folder (this script lives in tutorial/scripts/)
TINY="$HERE/backup/01.gm12878.bam"        # the starting BAM, shipped in the repo for this check
# a small public file in the lab's bucket: if we can download it, the bucket is reachable
BUCKET_TEST=https://s3.kopah.uw.edu/stergachis/public/FIRE/broadly-consented/GM12878/GM12878-fire-v0.1-peaks.bed.gz.tbi
fails=0                                   # number of failed checks; becomes the exit code

# --- helpers: print one PASS or FAIL line (FAIL adds a hint) ---
ok()   { printf '  PASS  %s\n' "$1"; }
bad()  { printf '  FAIL  %s\n        -> %s\n' "$1" "$2"; fails=$((fails + 1)); }
# run a command quietly (output saved to $W/log.txt); PASS if it succeeds, FAIL with the log's last line if not
check() { # check "<label>" "<hint if it fails>" command...
  local label=$1 hint=$2; shift 2
  if "$@" > "$W/log.txt" 2>&1; then ok "$label"; else bad "$label" "$hint (log: $(tail -1 "$W/log.txt"))"; fi
}

# --- 1. operating system + chip (fibertools builds: Mac Intel/Apple silicon, Linux x86-64) ---
echo "System: $(uname -s) $(uname -m)"
case "$(uname -s)-$(uname -m)" in
  Darwin-arm64|Darwin-x86_64|Linux-x86_64) ok "operating system supported" ;;
  Linux-aarch64) bad "operating system" "fibertools has no ARM Linux build (e.g. Windows on ARM / Snapdragon). Use the cluster instead (README step 0)." ;;
  *) bad "operating system" "not supported here. On Windows, run this inside WSL 2 (Ubuntu); see README step 0." ;;
esac
# --- 2. the conda environment: is the `ft` on your PATH the one inside the active environment? ---
if [ -n "${CONDA_PREFIX:-}" ] && [ "$(command -v ft)" = "$CONDA_PREFIX/bin/ft" ]; then
  ok "conda environment active ($CONDA_PREFIX)"
else
  bad "conda environment active" "run: conda activate fiberseq-workshop"
fi

# temporary working folder, deleted automatically when the script exits (trap ... EXIT)
W=$(mktemp -d "${TMPDIR:-/tmp}/fiberseq-check.XXXXXX")
trap 'rm -rf "$W"' EXIT

# --- 3. every program the tutorial calls is on the PATH ---
echo "Tools:"
for t in ft samtools tabix fiberhmm-call fiberhmm-extract fiberbrowser python; do
  if command -v "$t" > /dev/null 2>&1; then ok "$t found"; else bad "$t found" "missing: recreate the environment from environment.yml"; fi
done
# the Python libraries our step 9-11 scripts import
check "python packages (pandas, scipy, pysam)" "recreate the environment" python -c "import pandas, scipy, pysam"

# --- 4. a mini tutorial on real data ---
echo "A small run (20 fibers):"
if [ -s "$TINY" ]; then
  # 20 fibers from the NAPA region of the starting BAM: take the reads over chr19:47515000-47516000,
  # keep a random 20% of them (-s 0.2), index, and make sure at least one read came out
  check "samtools: read the BAM" "is backup/01.gm12878.bam complete? re-clone the repo" \
    bash -c "samtools view -b '$TINY' chr19:47515000-47516000 | samtools view -b -s 0.2 -o '$W/tiny.bam' - && samtools index '$W/tiny.bam' && [ \$(samtools view -c '$W/tiny.bam') -gt 0 ]"
  # the fibertools branch (tutorial steps 4 and 5)
  check "ft add-nucleosomes" "fibertools problem" ft add-nucleosomes "$W/tiny.bam" "$W/nuc.bam"
  check "ft fire" "fibertools problem" ft fire "$W/nuc.bam" "$W/fire.bam"
  # the FiberHMM branch (tutorial step 11), 2 cores
  check "fiberhmm-call (TF footprints)" "FiberHMM problem" fiberhmm-call -i "$W/tiny.bam" -o "$W/fhmm.bam" --enzyme hia5 --seq pacbio -c 2 --no-qc
  # the tools must also have written features: FIREs (fire.Q: in the Ma tag) and footprints (tf.QQQ: in MA)
  check "FIREs and footprints written" "tools ran but wrote no features" \
    bash -c "samtools view '$W/fire.bam' | grep -q 'fire.Q:' && samtools view '$W/fhmm.bam' | grep -q 'tf.QQQ:'"
  # FiberBrowser (step 12b): just check it imports; we don't open a browser here
  check "fiberbrowser starts" "FiberBrowser problem" python -c "import fiberbrowser.browser_client"
else
  bad "tutorial data" "backup/01.gm12878.bam not found: run this from the cloned repo's tutorial/ folder"
fi

# --- 5. internet: the lab bucket (data, hub) and UCSC (reference genome) ---
echo "Internet:"
check "lab bucket reachable" "can't reach s3.kopah.uw.edu (firewall or proxy?)" curl -sSf -o /dev/null --max-time 20 "$BUCKET_TEST"
# -r 0-99: download only the first 100 bytes of chr15, enough to prove the site answers
check "UCSC reachable" "can't reach hgdownload.soe.ucsc.edu" curl -sSf -o /dev/null --max-time 20 -r 0-99 https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr15.fa.gz

# --- summary; exit code = number of failures (0 = all good) ---
echo
if [ "$fails" -eq 0 ]; then echo "ALL PASS: you're ready."; else echo "$fails check(s) FAILED: send this whole output to the instructors."; fi
exit $fails
