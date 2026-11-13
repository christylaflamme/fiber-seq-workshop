# An introductory workshop for analyzing Fiber-seq data

A 1-hour introduction to analyzing **Fiber-seq** data: single-molecule, long-read chromatin architecture. Workshop was originally written for the Cold Spring Harbor Laboratory (CSHL) Computational Genomics course (2026) and assumes no prior experience with Fiber-seq or the command line.

## 👉 Start here: [the hands-on tutorial](tutorial/README.md)

**[`tutorial/README.md`](tutorial/README.md) is the workshop:** [0. Setup](tutorial/README.md#0-setup-before-the-workshop) (do before the session if possible);
then steps 1–12 from single fibers to haplotype-selective chromatin and TF footprints.

*Tip: the `.md` files here can also be viewed with a Markdown reader, e.g. the **Markdown Reader** extension for Chrome (or read them on GitHub).*

## What you'll learn
1. What Fiber-seq measures and how (m6A stencilling + long-read sequencing)
2. How a single read gets annotated for nucleosomes, accessible patches (MSPs) and regulatory elements (FIREs)
3. How individual fibers pileup to become genome tracks: peaks and **% accessible**
4. How long reads can resolve two distinct haplotypes in the human genome, thereby allowing us to ascertain haplotype-selective accessibility (e.g. imprinting, X inactivation)
5. What else the same data can do: CpG methylation, footprinting, and related methods

## Contents
- **[`tutorial/`](tutorial/README.md): the hands-on workshop (start here)**
- [`slides/`](slides/): lecture slides
- [`qc/GM12878_report/GM12878.overview.html`](qc/GM12878_report/GM12878.overview.html): precomputed fiberseq-qc report for GM12878
- [`handout/cheatsheet.md`](handout/cheatsheet.md): one-page glossary and how to read the tracks
- [Resources](tutorial/README.md#resources) (end of the tutorial): papers, software and data

## Credits
Demo data: public GM12878 Fiber-seq from the Stergachis lab, processed with the [FIRE pipeline](https://github.com/fiberseq/FIRE), created by Mitchell Vollger. [FiberHMM](https://github.com/fiberseq/FiberHMM) was created by Tommy Tullius. See the [tutorial's Resources](tutorial/README.md#resources) for citations.

Thank you to Andrew Stergachis, Mitchell Vollger and Tommy Tullius for their input and expertise.

## UCSC extended demo *(optional)*
👉 **[GM12878 Fiber-seq in the UCSC Genome Browser](https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&hubClear=https://s3.kopah.uw.edu/stergachis/public/FIRE/broadly-consented/GM12878/trackHub-v0.1/hub.txt&udcTimeout=1&position=chr12:6532500-6536500&GM12878-H1-FIRE-fibers=pack&GM12878-H2-FIRE-fibers=pack)** (hg38)

Then follow the stops in [`UCSC_demo/UCSC_demo_guide.md`](UCSC_demo/UCSC_demo_guide.md).

## FiberBrowser tutorial *(optional)*
👉 **[FiberBrowser tutorial](https://github.com/fiberseq/FiberBrowser/blob/main/docs/tutorial.md)**: a guided tour of single-molecule footprint analysis in FiberBrowser.

## License
[CC BY 4.0](LICENSE)

---

*Materials developed in collaboration with [Claude Code](https://claude.com/claude-code) (Anthropic).*
