# Extended demo: Fiber-seq in the UCSC Genome Browser

**Data:** public Fiber-seq of the GM12878 lymphoblastoid cell line (a B-cell line from a female donor), processed with the FIRE pipeline. It's shown as a UCSC track hub on **hg38**. Nothing to install: click a link.

**Time:** ~20 minutes, 5 stops. Each stop has a link that opens the hub at that locus with the single-fiber tracks expanded.

> All coordinates below are **hg38**. Accessibility numbers are the maxima of the hub's own `% accessible` tracks over each element. We read them from the public bigWig/bigBed files on 2026-10-27; they are not new analysis.

---

## Before you start: the tracks

| Track | What it shows |
|---|---|
| **FIRE fibers (H1 / H2 / UNK)** | One row per sequenced molecule (fiber), split by haplotype (H1, H2, or unphased). **Red** = FIRE (a regulatory element on that fiber), **purple** = MSP (accessible patch that isn't a FIRE), **gray** = nucleosome |
| **% accessible** | % of fibers that have a FIRE overlapping each base. **Black** = all fibers, **red** = haplotype 1, **blue** = haplotype 2 |
| **FIRE peaks** | Peaks called where the FIRE score is a local maximum with FDR < 5% |
| **Wide FIRE peaks** | Peaks plus all regions below the FDR threshold, merged within 147 bp |
| **FIRE hap-differences** | Per-peak test of hap1 vs hap2 accessibility. ⚠ **Colors are reversed from % accessible here:** **blue** = more accessible on hap1, **red** = more accessible on hap2. Bright = passes the pipeline's significance threshold, dark = doesn't |
| **FIRE coverage** | Number of fibers with a FIRE / MSP / nucleosome at each base |
| **Unreliable coverage regions** | Excluded from FDR calculation (depth > 5 SD from the median) |

Track definitions are from the hub's own description page (`fire-description.html`). The hap-differences colors come from the FIRE pipeline source ([`workflow/scripts/hap-diffs.R`](https://github.com/fiberseq/FIRE/blob/main/workflow/scripts/hap-diffs.R)). The item name is `p-value_hap1%_hap2%`; click a feature to see it.

**Haplotype labels are arbitrary.** "H1" and "H2" come from phasing. They don't mean maternal and paternal unless the phasing used parental data.

---

## Stop 1: Orientation, a housekeeping promoter (*GAPDH*) — [open][gapdh]

`chr12:6,532,500-6,536,500`

- Walk through each track, top to bottom.
- Zoom in on the H1 fibers: each row is **one DNA molecule** from **one cell**. Find the regularly spaced gray blocks (nucleosomes, ~147 bp) and the red patch at the promoter.
- *Ask the room:* The promoter is open on ~78% of fibers. Why not 100%?
  - *Answer:* Each fiber is a snapshot from one cell at one moment. Even "always on" elements aren't open on every molecule. Bulk ATAC/DNase can't give you this number.

## Stop 2: On vs off, cell-type specificity — [*CD79A*][cd79a] vs [*ALB*][alb]

`chr19:41,872,000-41,884,000` and `chr4:73,398,000-73,410,000`

- *CD79A* (a B-cell receptor component) has an open promoter on ~70% of fibers.
- *ALB* (albumin, a liver gene) is ~1%: a solid wall of gray nucleosomes.
- Optional: turn on **ENCODE Regulation → DNase** (GM12878) in UCSC to compare Fiber-seq with bulk accessibility.

## Stop 3: Imprinting, one allele open and one closed — [*KCNQ1OT1* / KvDMR][kvdmr]

`chr11:2,695,000-2,705,000`

- The *KCNQ1OT1* promoter sits in an imprinting control region (KvDMR1). It is active only on the paternally inherited allele.
- `% accessible`: hap1 ~2% and hap2 ~88%. The hap-differences track shows a **bright red** (hap2-selective) feature here.
- Look at the fibers: in H2 nearly every fiber has a red FIRE, and in H1 almost none do.
- *Ask:* What would bulk ATAC show here? (~46%. It looks like a "medium" peak, and the two-allele story is invisible.)

## Stop 4: X inactivation — [*XIST*][xist], [*HPRT1*][hprt1], [*KDM6A*][kdm6a]

GM12878 is female: one X is active (Xa), one is inactive (Xi).

| Locus (hg38) | hap1 | hap2 | Reading |
|---|---|---|---|
| *XIST* promoter, `chrX:73,848,000-73,856,000` | ~11% | ~66% | *XIST* is expressed from the **Xi** → here hap2 behaves like the Xi |
| *HPRT1* promoter, `chrX:134,456,000-134,466,000` | ~81% | ~10% | silenced on the Xi → hap1 behaves like the Xa |
| *KDM6A* promoter, `chrX:44,868,000-44,878,000` | ~93% | ~87% | a known **escape gene**: open on both X's |

- *Ask:* Why are XIST and HPRT1 open on *opposite* haplotypes?
- In this sample the pattern lines up across >60 Mb of chrX (hap2 = Xi-like at both loci). That's an observation from these three loci, not a claim the pipeline makes.

## Stop 5: Explore

Type any gene into the UCSC search box. Things to try:
- A gene you work on: is its promoter open in B-cells?
- Zoom out to ~200 kb and compare FIRE peaks with ENCODE cCREs.
- Right-click a fiber track → configure → turn on **m6A** or **5mC** decorations to see the raw marks.

---

**Start-page link** (hub with default track settings): [open][start]

[start]: https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&hubClear=https://s3.kopah.uw.edu/stergachis/public/FIRE/broadly-consented/GM12878/trackHub-v0.1/hub.txt&udcTimeout=1&position=chr12:6532500-6536500
[gapdh]: https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&hubClear=https://s3.kopah.uw.edu/stergachis/public/FIRE/broadly-consented/GM12878/trackHub-v0.1/hub.txt&udcTimeout=1&position=chr12:6532500-6536500&GM12878-H1-FIRE-fibers=pack&GM12878-H2-FIRE-fibers=pack
[cd79a]: https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&hubClear=https://s3.kopah.uw.edu/stergachis/public/FIRE/broadly-consented/GM12878/trackHub-v0.1/hub.txt&udcTimeout=1&position=chr19:41872000-41884000&GM12878-H1-FIRE-fibers=pack&GM12878-H2-FIRE-fibers=pack
[alb]: https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&hubClear=https://s3.kopah.uw.edu/stergachis/public/FIRE/broadly-consented/GM12878/trackHub-v0.1/hub.txt&udcTimeout=1&position=chr4:73398000-73410000&GM12878-H1-FIRE-fibers=pack&GM12878-H2-FIRE-fibers=pack
[kvdmr]: https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&hubClear=https://s3.kopah.uw.edu/stergachis/public/FIRE/broadly-consented/GM12878/trackHub-v0.1/hub.txt&udcTimeout=1&position=chr11:2695000-2705000&GM12878-H1-FIRE-fibers=pack&GM12878-H2-FIRE-fibers=pack
[xist]: https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&hubClear=https://s3.kopah.uw.edu/stergachis/public/FIRE/broadly-consented/GM12878/trackHub-v0.1/hub.txt&udcTimeout=1&position=chrX:73848000-73856000&GM12878-H1-FIRE-fibers=pack&GM12878-H2-FIRE-fibers=pack
[hprt1]: https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&hubClear=https://s3.kopah.uw.edu/stergachis/public/FIRE/broadly-consented/GM12878/trackHub-v0.1/hub.txt&udcTimeout=1&position=chrX:134456000-134466000&GM12878-H1-FIRE-fibers=pack&GM12878-H2-FIRE-fibers=pack
[kdm6a]: https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&hubClear=https://s3.kopah.uw.edu/stergachis/public/FIRE/broadly-consented/GM12878/trackHub-v0.1/hub.txt&udcTimeout=1&position=chrX:44868000-44878000&GM12878-H1-FIRE-fibers=pack&GM12878-H2-FIRE-fibers=pack
