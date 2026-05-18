# Proposed diagnostic probe set — V3-V4 16S family-specific probes

Derived from 7,646 ASVs (QIIME2 / SILVA, V3-V4 amplicon 341F–805R).
Targets the four top Random Forest predictors distinguishing
AB_lower_contact (Muzzle + Pastern + Udder + Ventral abdomen)
from C_elevated (Dorsum + Forehead + Neck + Pectoral area);
OOB accuracy 97.9%, AUC = 0.999.

## Probe sequences

| ID       | Target family            | Sequence (5′→3′)         | Len | GC% | Tm (°C) |
|----------|--------------------------|--------------------------|-----|-----|---------|
| DEIN-V3  | Deinococcaceae           | CAGCCGCGGTAATACGGAGG     | 20  | 65  | 52.8    |
| INTR-V3  | Intrasporangiaceae       | TATTGGGCGTAAAGAGCTTG     | 20  | 45  | 44.6    |
| HUNG-V3  | Hungateiclostridiaceae   | TTACTGGGTGTAAAGGGCGTGT   | 22  | 50  | 49.7    |
| CHRI-V3  | Christensenellaceae      | AGGAAGCCCCGGCTAACTACGT   | 22  | 59  | 53.4    |

Tm: Marmur–Schildkraut–Doty formula, 50 mM NaCl (Primer3 standard conditions).

## Specificity (within-study evaluation)

Fraction of ASVs from each family matching the probe (exact substring or reverse complement):

| Probe    | Deinococcaceae | Intrasporangiaceae | Hungateiclostridiaceae | Christensenellaceae | All 7646 ASVs |
|----------|:-:|:-:|:-:|:-:|:-:|
| DEIN-V3  | **0.97** | 0.00 | 0.00 | 0.00 | 6.4% |
| INTR-V3  | 0.00 | **1.00** | 0.00 | 0.00 | 2.1% |
| HUNG-V3  | 0.00 | 0.00 | **0.82** | 0.00 | 4.2% |
| CHRI-V3  | 0.00 | 0.00 | 0.00 | **0.90** | 8.3% |

All four probes have zero cross-reactivity among the four target families.
Background hit rates (all ASVs) are low (2–8%), reflecting probe specificity
beyond the four-family panel.

## Biological context

| Probe    | Group enrichment  | log2FC (AB/C) | Source       |
|----------|-------------------|---------------|--------------|
| DEIN-V3  | **C_elevated**    | −4.30 ***     | env-shared   |
| INTR-V3  | **AB_lower**      | +1.56 ***     | animal-only  |
| HUNG-V3  | **AB_lower**      | +1.96 ***     | animal-only  |
| CHRI-V3  | **AB_lower**      | +2.40 ***     | env-shared   |

BH-adjusted p-values from Wilcoxon AB vs C (n=48 per group).

## Proposed applications

1. **TaqMan qPCR** — probes with fluorescent label + quencher flanked by
   universal 16S primers (e.g., 341F / 805R). Tm gap to universal primers
   ~10–15°C; LNA modifications (2–4 substitutions) recommended to raise
   probe Tm to 60–65°C.

2. **Hybridisation capture** — biotinylated 20–22-mers for selective enrichment
   of target-family amplicons before short-read sequencing (MyBaits / xGen style).

3. **Rapid field diagnostic** — lateral-flow or bead-array with all four probes
   simultaneously; a AB-pattern (INTR+, HUNG+, CHRI+, DEIN−) vs
   C-pattern (DEIN+, others−) gives a single contact-zone call.

## Caveats and required validation

- Within-study specificity only; full validation requires BLAST against
  the complete SILVA 138 NR99 database (16S) and NCBI nt.
- Sensitivity for Hungateiclostridiaceae (82%) and Christensenellaceae (90%)
  reflects within-family sequence diversity; probe pools (2–3 oligos per family)
  would improve coverage.
- Primer3 / BLAST-based primer design recommended before synthesis.
- Experimental validation (in vitro PCR, Sanger confirmation) required
  before clinical use.

## Derivation notes

- Sequences sourced from QIIME2 DADA2 ASVs (rep-seqs.qza), V3-V4 region.
- Alignment: MAFFT v7.505 (--auto mode).
- Conserved-window search: 20–22 nt sliding window, ≥85% column conservation,
  38–62% GC (before family-discriminability filter).
- Cross-family specificity computed as exact substring match (or RC) against
  all ASVs from each of the four target families.
- Script: probe_design steps embedded in picrust2_analysis.R session.
