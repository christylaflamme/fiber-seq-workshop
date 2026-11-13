#!/usr/bin/env bash
# Instructor only: build the UCSC track hub students open in step 12a, laid out like the lab's
# production (prod-reporter) hubs: the FIRE pipeline's GM12878 trackDb, finalized as production does
# (hidden tracks shown), plus production's CpG track layout, as one useOneFile hub.txt.
#   conda activate fiberseq-workshop
#   HUB_EMAIL=<contact> bash scripts/make_hub.sh <step outputs dir> <hub dir> [threads]
# Track data:
#   built here from the tutorial steps: fibers + decorators (05.fire.bam), FIRE/linker/nucleosome
#     coverage and % accessible (07.pileup.bed), with the FIRE pipeline's formulas (workflow/rules/track-hub.smk)
#   cut from the public genome-wide pileup: -log10 FDR (log_FDR column), for SNRPN and NAPA
#   the FIRE pipeline's own genome-wide GM12878 files, by URL: narrow + wide peaks, hap-differences,
#     unreliable-coverage regions
#   CpG (combined/H1/H2) and CpG differences (hap1 - hap2) from step 10 (10.cpg-by-hap.tsv: ft pileup --cpg
#     + our cpg_by_hap.py). Production uses pb-CpG-tools + methbat instead (Linux only), with p-value tiers
#     on the differences; we show the difference only.
#
# Output (<hub dir>):
#   hub.txt                    the hub: header + every track stanza in one file (useOneFile)
#   fire-description.html      the hub's description page
#   bb/  fire-fibers.bb, fire-fiber-decorators.bb   one row per fiber (bigBed)
#   bw/  coverage, % accessible, log_FDR and CpG signal tracks (bigWig)
# Templates it reads from tutorial/hub/: fire-trackDb.GM12878.txt, cpg-tracks.template, *.as, fire-description.html
#
# How UCSC tracks are made here: write a bedGraph (chrom, start, end, value) with awk, then convert it to
# bigWig with bedGraphToBigWig; feature BEDs are converted to bigBed with bedToBigBed. Both need chrom.sizes.

# stop at the first error, on unset variables, and on failures inside pipes
set -euo pipefail
IN=$(cd "${1:?usage: make_hub.sh <step outputs dir> <hub dir> [threads]}" && pwd)   # step outputs (e.g. backup/)
OUT=${2:?usage: make_hub.sh <step outputs dir> <hub dir> [threads]}                  # where to write the hub
T=${3:-4}                                                                             # threads
TPL=$(cd "$(dirname "$0")/../hub" && pwd)                                             # tutorial/hub/ templates
EMAIL=${HUB_EMAIL:-YOUR_EMAIL}                                                        # contact shown in UCSC
SAMPLE=GM12878
# the FIRE pipeline's public GM12878 results (genome-wide tracks are linked, not copied)
GW=https://s3.kopah.uw.edu/stergachis/public/FIRE/broadly-consented/GM12878
REGIONS="chr15:24850000-25050000 chr19:47480000-47525000"
# macOS only: bioconda's UCSC tools look for libpng in a Homebrew path; point them at the env's copy
if [ "$(uname)" = Darwin ] && [ -n "${CONDA_PREFIX:-}" ]; then export DYLD_FALLBACK_LIBRARY_PATH="$CONDA_PREFIX/lib"; fi

# fresh output folders + a scratch folder inside them (deleted at the end)
rm -rf "$OUT" && mkdir -p "$OUT/bb" "$OUT/bw"
OUT=$(cd "$OUT" && pwd)
TMP="$OUT/.tmp" && mkdir -p "$TMP"
# hg38 chromosome lengths, required by bedToBigBed / bedGraphToBigWig
curl -sSf -o "$TMP/chrom.sizes" https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes

# fibers + decorators (ft track-decorators, as in the FIRE pipeline)
# --bed12: one BED12 line per fiber (the fiber bar); --decorator: the FIREs, linkers, nucleosomes,
# m6A and 5mC drawn on top of it. Sorted, then converted to bigBed with the pipeline's .as schemas.
ft track-decorators -t $T --bed12 "$TMP/fibers.bed" --decorator "$TMP/decorators.bed" "$IN/05.fire.bam"
sort -k1,1 -k2,2n "$TMP/fibers.bed" > "$TMP/fibers.s.bed"
sort -k1,1 -k2,2n "$TMP/decorators.bed" > "$TMP/decorators.s.bed"
bedToBigBed -type=bed12+1 -as="$TPL/bed12_filter.as" "$TMP/fibers.s.bed" "$TMP/chrom.sizes" "$OUT/bb/fire-fibers.bb"
bedToBigBed -type=bed12+5 -as="$TPL/decoration.as" "$TMP/decorators.s.bed" "$TMP/chrom.sizes" "$OUT/bb/fire-fiber-decorators.bb"

# coverage tracks from the step 7 pileup (columns: coverage 4, fire 5, nuc 7, msp 8; _H1 9-13; _H2 14-18)
# NR>1 skips the header line. Formulas follow the FIRE pipeline:
#   FIRE coverage = fire; linker coverage = msp - fire (MSPs that are not FIREs); nucleosome coverage = nuc
#   % accessible = 100 * fire / coverage (only where coverage > 0), for all fibers, hap1 and hap2
P="$IN/07.pileup.bed"
awk 'NR>1 {print $1"\t"$2"\t"$3"\t"$5}'    "$P" > "$TMP/all.fire.bg"
awk 'NR>1 {print $1"\t"$2"\t"$3"\t"$8-$5}' "$P" > "$TMP/all.linker.bg"
awk 'NR>1 {print $1"\t"$2"\t"$3"\t"$7}'    "$P" > "$TMP/all.nucleosome.bg"
awk 'NR>1 && $4>0  {print $1"\t"$2"\t"$3"\t"100*$5/$4}'   "$P" > "$TMP/all.pa.bg"
awk 'NR>1 && $9>0  {print $1"\t"$2"\t"$3"\t"100*$10/$9}'  "$P" > "$TMP/hap1.pa.bg"
awk 'NR>1 && $14>0 {print $1"\t"$2"\t"$3"\t"100*$15/$14}' "$P" > "$TMP/hap2.pa.bg"
# bedGraph -> bigWig (file names match the bigDataUrl lines in the trackDb template)
for e in fire linker nucleosome; do bedGraphToBigWig "$TMP/all.$e.bg" "$TMP/chrom.sizes" "$OUT/bw/all.$e.coverage.bw"; done
for h in all hap1 hap2; do bedGraphToBigWig "$TMP/$h.pa.bg" "$TMP/chrom.sizes" "$OUT/bw/$h.percent.accessible.bw"; done

# -log10 FDR from the public genome-wide pileup (the public hub lacks log_FDR.bw)
# tabix pulls our two regions over the internet; the remote index it caches is deleted.
# awk finds the log_FDR column by its header name, then prints chrom/start/end/log_FDR as a bedGraph.
( cd "$TMP" && tabix -h "$GW/GM12878-fire-v0.1-pileup.bed.gz" $REGIONS > gw.pileup.bed && rm -f GM12878-fire-v0.1-pileup.bed.gz.tbi )
awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) if ($i=="log_FDR") c=i; next} {print $1"\t"$2"\t"$3"\t"$c}' "$TMP/gw.pileup.bed" \
  | sort -k1,1 -k2,2n > "$TMP/log_FDR.bg"
bedGraphToBigWig "$TMP/log_FDR.bg" "$TMP/chrom.sizes" "$OUT/bw/log_FDR.bw"

# trackDb: FIRE's GM12878 trackDb, finalized like production (hidden -> full), genome-wide tracks by URL
# sed edits, in order: show every track (visibility hide -> full); strip trailing spaces; point the four
# genome-wide bigBeds (peaks, wide peaks, hap-differences, unreliable regions) at the public pipeline files
# instead of local copies.
sed -E -e 's/visibility hide/visibility full/' -e 's/[[:space:]]+$//' \
  -e "s#bigDataUrl bb/(fire-peaks|fire-wide-peaks|hap_differences|unreliable-coverage-regions)\.bb#bigDataUrl $GW/trackHub-v0.1/bb/\1.bb#" \
  "$TPL/fire-trackDb.GM12878.txt" > "$TMP/trackDb.txt"
# CpG tracks from step 10 (10.cpg-by-hap.tsv columns: pct_all 6, pct_H1 9, pct_H2 12)
# one bedGraph per column (combined, hap1, hap2), skipping CpGs with no value ("nan")
C="$IN/10.cpg-by-hap.tsv"
for x in combined:6 hap1:9 hap2:12; do
  awk -F'\t' -v c=${x#*:} 'NR>1 && $c!="" && $c!="nan" {print $1"\t"$2"\t"$3"\t"$c}' "$C" > "$TMP/cpg.${x%:*}.bg"
done
# CpG differences: hap1 % - hap2 % (only CpGs with a value on both haplotypes)
awk -F'\t' 'NR>1 && $9!="" && $9!="nan" && $12!="" && $12!="nan" {print $1"\t"$2"\t"$3"\t"$9-$12}' "$C" > "$TMP/cpg.diffs_all.bg"
for x in combined hap1 hap2 diffs_all; do bedGraphToBigWig "$TMP/cpg.$x.bg" "$TMP/chrom.sizes" "$OUT/bw/cpg.$x.bw"; done
# append production's CpG track stanzas to the trackDb, with ${SAMPLE} filled in
printf '\n' >> "$TMP/trackDb.txt"
sed "s/\${SAMPLE}/$SAMPLE/g" "$TPL/cpg-tracks.template" >> "$TMP/trackDb.txt"

# one-file hub, as prod-reporter publishes them
# hub header (name, labels, useOneFile, contact, description page, genome), then all the track stanzas
{
  printf 'hub %s-FIRE-fiberseq\nshortLabel %s-FIRE-fiberseq\nlongLabel %s-FIRE-fiberseq\nuseOneFile on\nemail %s\ndescriptionUrl fire-description.html\n\ngenome hg38\n\n' \
    "$SAMPLE" "$SAMPLE" "$SAMPLE" "$EMAIL"
  cat "$TMP/trackDb.txt"
} > "$OUT/hub.txt"
cp "$TPL/fire-description.html" "$OUT/"
rm -rf "$TMP"

# every local bigDataUrl must exist
# (a track pointing at a missing file shows as an error in UCSC; URLs to the public pipeline are not checked)
awk '/bigDataUrl/ {print $2}' "$OUT/hub.txt" | grep -v '^https://' | while read -r f; do
  [ -s "$OUT/$f" ] || { echo "MISSING: $f" >&2; exit 1; }
done
# list what was written
find "$OUT" -type f | sort
