# Skin microbiome of healthy Norik Muránský draught mares

**Body-site environmental contact determines skin microbiome composition in healthy Norik Muránský mares, with identification of a four-family diagnostic probe panel**

*Manuscript in preparation*

---

## Study overview

We characterised the baseline skin microbiome of 12 healthy adult Norik Muránský mares
(Dobšiná, Slovakia; age 5–15 years, mean 9.6 ± 3.4 SD) at eight anatomical sites plus
one environmental control, using 16S rRNA V3–V4 amplicon sequencing (QIIME2/DADA2,
SILVA 138) cross-validated with Kraken2 (confidence > 33%, min-hit-groups = 5).

**Key findings**

| Finding | Result |
|---|---|
| Primary community driver | Body-site topology (GCtS vs EtS), not age |
| Family richness GCtS vs EtS | Wilcoxon W = 1435, p = 0.038; mean 37.5 vs 33.4 families |
| Families differing GCtS vs EtS | 25 / 46 (BH-corrected Wilcoxon) |
| RF OOB accuracy (binary GCtS vs EtS) | 97.9% ± 0.73% SD (1000 trees, 10 seeds) |
| AUC (bootstrap n = 2000) | 0.999 (95% CI 0.995–1.000) |
| Youden's J | 0.979 (95% CI 0.939–1.000) |
| Horse identity (ICC) | < 0.01 — negligible batch effect |
| PICRUSt2 pathways significant | 262 / 540 (BH-corrected Wilcoxon) |

**Topology grouping**

| Group | Body sites | Ecological rationale |
|---|---|---|
| **GCtS** (ground-contact) | Pastern, Muzzle, Ventral abdomen, Udder | Direct soil / grass / faecal contact |
| **EtS** (elevated) | Dorsum, Forehead, Neck, Pectoral area | UV exposure, sweat, harness contact |

---

## Repository structure

```
.
├── README.md
├── .gitignore
│
├── data/
│   ├── metadata.csv                     # sample metadata (horse ID, site, age, sex)
│   ├── pca_family_contributions.csv     # normalised family relative abundances
│   └── taxonomy_export/
│       └── taxonomy.tsv                 # ASV → taxonomy (SILVA 138)
│
├── scripts/
│   ├── 01_richness_glm.R                # GLM richness ~ topology + age
│   ├── 02_anova_richness_topology.R     # ANOVA full model (with age)
│   ├── 03_anova_richness_no_age.R       # ANOVA without age; Tukey CLD; richness plots
│   ├── 04_wilcoxon_young_old.R          # Udder and Dorsum richness: young vs old
│   ├── 05_family_abundance_tests.R      # Kruskal-Wallis + Wilcoxon per family, BH
│   ├── 06_family_abundance_figure.R     # 4-panel abundance figure
│   ├── 07_rf_topology_predictors.R      # RF multiclass; MDA importance
│   ├── 08_rf_validation.R               # Multinomial, ROC, PCA, Kraken2 overlap
│   ├── 09_rf_diagnostics.R              # Convergence, stability, calibration, Youden
│   ├── 10_roc_smoothed.R                # Smoothed ROC with bootstrap CI
│   ├── 11_analysis_ACgroups.R           # Binary RF GCtS vs EtS; abundance; PCA
│   ├── 12_picrust2_install_run.sh       # PICRUSt2 installation + pipeline
│   ├── 13_picrust2_analysis.R           # MetaCyc pathway analysis
│   ├── 14_lmm_sensitivity.R             # ICC, LMM, within-horse permutation, noise injection
│   ├── 15_alternative_panel_validation.R# Panel A vs Panel B ROC, Youden, LOHO, noise
│   ├── 16_probe_population_validation.R # Per-probe ROC, permutation, LOHO cross-validation
│   └── 17_hypothesis_figure.py          # Conceptual figure: EtS dysbiosis hypothesis
│
├── probe_design/
│   └── diagnostic_probe_set.md          # Probe sequences and design notes
│
└── results/
    ├── figures/                         # Publication-ready figures (PDF + PNG)
    └── tables/                          # Result tables (CSV)
```

> **Raw reads** are deposited at NCBI SRA under accession **[TBD upon acceptance]**
> and are excluded from this repository.

---

## How to reproduce

### Prerequisites

**R packages:**

```r
install.packages(c("tidyverse", "randomForest", "pROC", "nnet",
                   "patchwork", "multcomp", "MASS", "lme4", "lmerTest"))
```

**Python packages (script 17):**

```bash
pip install matplotlib
```

**PICRUSt2** (conda, ~2 GB):

```bash
bash scripts/12_picrust2_install_run.sh
```

### Run order

```bash
Rscript scripts/01_richness_glm.R
Rscript scripts/02_anova_richness_topology.R
Rscript scripts/03_anova_richness_no_age.R
Rscript scripts/04_wilcoxon_young_old.R
Rscript scripts/05_family_abundance_tests.R
Rscript scripts/06_family_abundance_figure.R
Rscript scripts/07_rf_topology_predictors.R
Rscript scripts/08_rf_validation.R
Rscript scripts/09_rf_diagnostics.R
Rscript scripts/10_roc_smoothed.R
Rscript scripts/11_analysis_ACgroups.R
Rscript scripts/14_lmm_sensitivity.R

# After PICRUSt2 pipeline completes (~30–90 min):
bash scripts/12_picrust2_install_run.sh
Rscript scripts/13_picrust2_analysis.R

# Probe panel validation:
Rscript scripts/15_alternative_panel_validation.R
Rscript scripts/16_probe_population_validation.R

# Discussion figure:
python3 scripts/17_hypothesis_figure.py
```

All R scripts use `set.seed(42)` for reproducibility.

---

## Diagnostic probe panels

**Panel A (original — top RF predictors):**

| ID | Family | Sequence 5′→3′ | Tm (°C) | GC% |
|---|---|---|---|---|
| DEIN-V3 | *Deinococcaceae* | `CAGCCGCGGTAATACGGAGG` | 52.8 | 65 |
| INTR-V3 | *Intrasporangiaceae* | `TATTGGGCGTAAAGAGCTTG` | 44.6 | 45 |
| HUNG-V3 | *Hungateiclostridiaceae* | `TTACTGGGTGTAAAGGGCGTGT` | 49.7 | 50 |
| CHRI-V3 | *Christensenellaceae* | `AGGAAGCCCCGGCTAACTACGT` | 53.4 | 59 |

**Panel B (generalised — geographically widespread families):**
Deinococcaceae + Micrococcaceae + Lachnospiraceae + UCG-010
(OOB 97.9%, AUC 0.998; equivalent performance to Panel A with potential
cross-population applicability)

Classification rule: DEIN-V3 positive + gut-family probes negative → **EtS**.
Reverse pattern → **GCtS**.

---

## Data availability

| Resource | Location |
|---|---|
| Raw reads | NCBI SRA [accession TBD] |
| Processed data | `data/` directory in this repository |
| PICRUSt2 outputs | Regenerate with `scripts/12_picrust2_install_run.sh` |

---

## Citation

> Kováčová V et al. (2026) Body-site environmental contact determines skin microbiome
> composition in healthy Norik Muránský draught mares. *[journal, in preparation]*.

---

## License

Code: MIT License.
Data: CC BY 4.0 (upon publication).
