# *SNRPN* track hub (tutorial step 12a)

Built by `scripts/make_hub.sh`, laid out like the Stergachis lab's production hubs:

- `fire-trackDb.GM12878.txt`: the FIRE pipeline's own trackDb for GM12878, copied unchanged from the public
  hub (`https://s3.kopah.uw.edu/stergachis/public/FIRE/broadly-consented/GM12878/trackHub-v0.1/trackDb.txt`).
  `make_hub.sh` finalizes it as production does (hidden tracks shown) and writes one `useOneFile` hub.txt.
- `cpg-tracks.template`: production's CpG and CpG-differences track stanzas (difference track only; production
  also shows p-value tiers from methbat).
- `fire-description.html`: the FIRE pipeline's hub description page (from the same public hub).
- `bed12_filter.as`, `decoration.as`: from the [FIRE pipeline](https://github.com/fiberseq/FIRE)
  (`workflow/templates`, MIT license).

Track data: fibers, coverage and % accessible from the tutorial's own steps (FIRE's formulas); CpG from step 10
(production computes CpG with pb-CpG-tools instead); FDR from the
public genome-wide pileup; peaks, wide peaks, hap-differences and unreliable-coverage regions are the FIRE
pipeline's genome-wide GM12878 files, linked by URL.
