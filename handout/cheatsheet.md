# Fiber-seq cheat sheet

## Glossary
| Term | Meaning |
|---|---|
| **Hia5** | A non-specific N6-adenine methyltransferase that adds m6A to accessible adenines. |
| **m6A** | N6-methyladenine. Nearly absent from human DNA, so any m6A on a Fiber-seq read was put there by Hia5 (accessible DNA). |
| **5mC / CpG methylation** | Endogenous epigenetic mark that can be read from the same molecule with long-read sequencing. |
| **Fiber** | One Hia5-treated and sequenced DNA molecule (usually 10–20+ kb). |
| **Chromatin accessibility** | Open/accessible patch of DNA. |
| **Nucleosome** | Histone octamer in which DNA is wound to form basis of chromatin. ~147 bp stretch of inaccessibility on a fiber. |
| **Methyltransferase-sensitive patch (MSP)** | Any accessible stretch on a fiber. |
| **Linker regions** | Short accessible patches (MSPs) between nucleosomes. |
| **FIRE** | MSP corresponding to a Fiber-seq Inferred Regulatory Element. |
| **Chromatin actuation** | FIRE coverage / total coverage (% of fibers with accessibility). |
| **FIRE peak** | Genomic position where FIREs pile up beyond chance (FDR < 5%). |
| **Footprint** | Short stretch of DNA protected from m6A modification (e.g. TF, Pol, nucleosome). |
| **Haplotype (H1/H2)** | Reads are assigned to haplotype 1 or haplotype 2 using heterozygous variants on the same read. Labels are arbitrary unless the phasing used parental data. |
| **Haplotype-selective element** | A peak that is significantly more accessible on one haplotype versus the other (imprinting, X inactivation, cis-acting variants). |

## Reading a fiber track (FIRE hub colors)
- **Red** = FIRE · **purple** = MSP (accessible, not a FIRE) · **gray** = nucleosome
- `% accessible` lines: **black** = all, **red** = hap1, **blue** = hap2
- `hap-differences`: **blue** = hap1 more accessible, **red** = hap2 more accessible (⚠ the opposite of the % accessible colors). Bright = significant.

## Fiber-seq vs. ATAC-seq
| | ATAC-seq | Fiber-seq |
|---|---|---|
| How accessible DNA is marked | Tn5 cuts and inserts adapters | Hia5 methylates A's (m6A); DNA stays intact |
| One read is… | A short fragment between two cuts | One chromatin fiber, many kb (median 16.5 kb in the tutorial data) |
| Signal | Read counts: relative, population average | % of fibers open: an absolute fraction |
| Nucleosomes | Inferred from fragment sizes | Called on every molecule |
| Haplotypes | Only reads that cover a het SNP | Almost every read (817 H1 / 825 H2 of 1,644 in the tutorial data) |
| Also on the read | — | 5mC (CpG methylation) on the same molecule |

ATAC-seq: Buenrostro et al., *Nat Methods* 2013 · Fiber-seq: Stergachis et al., *Science* 2020

## Typical analysis path (for later)
1. Sequence on PacBio (or ONT) to get a BAM with m6A + 5mC calls (`MM`/`ML`). Current Revio runs call m6A on the instrument (older runs: keep kinetics, then `ft predict-m6a`)
2. `ft add-nucleosomes`: nucleosomes and MSPs
3. Align and phase (e.g. pbmm2 + HiPhase), or align to the sample's own assembly
4. FIRE pipeline: FIREs, peaks, % accessible, haplotype tests, UCSC track hub
5. Look at it (UCSC, IGV), then ask your question

Details: [fiberseq.github.io](https://fiberseq.github.io/)
