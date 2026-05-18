# Skin microbiome of healthy cold-blooded draft horses

**Body-site contact zone determines skin microbiome composition in healthy Norik Muránsky mares, with identification of a four-family diagnostic marker panel**

*Manuscript in preparation — Frontiers in Veterinary Science*

---

## Study overview

We characterised the baseline skin microbiome of 12 healthy, adult Norik Muránsky mares
(Dobšiná, Slovakia) at eight anatomical sites plus one environmental control,
using 16S rRNA V3-V4 amplicon sequencing (QIIME2/DADA2, SILVA 138) cross-validated
with Kraken2 (confidence > 33%, min-hit-groups = 5).

**Key findings**

| Finding | Result |
|---|---|
| Primary community driver | Body-site contact zone (AB vs C), not age |
| Family richness by site | ANOVA F(7,88) = 10.69, p < 0.001; Dorsum highest |
| Age effect | Not significant (drop1 F-test, p = 0.46) |
| AB vs C richness | Wilcoxon W = 1435, p = 0.038 |
| Families differing AB vs C | 25 / 46 (BH-corrected Wilcoxon) |
| RF OOB accuracy (binary AB vs C) | 97.9% (1000 trees, 10-seed SD = 0.73%) |
| AUC (binary RF, bootstrap n = 2000) | 0.999 (95% CI: 0.995–1.000) |
| Youden's J | 0.979 (95% CI: 0.939–1.000) |
| Minimal predictor set | 4 families → 93.8% accuracy (multinomial logistic) |

**Contact-zone grouping**

| Group | Body sites |
|---|---|
| **GCtS** | Left front pastern, Muzzle, Ventral abdomen, Udder |
| **EtS** | Dorsum, Forehead, Neck, Pectoral area |

---

## Repository structure

```
.
├── README.md
├── .gitignore
│
├── data/
│   ├── metadata.csv                     # sample metadata (horse ID, site, age, sex)
│   ├── family_table.tsv                 # raw family-level count table (QIIME2 export)
│   ├── pca_family_contributions.csv     # normalised family relative abundances (input for all R scripts)
│   └── taxonomy_export/
│       └── taxonomy.tsv                 # ASV → taxonomy assignments (SILVA 138)
│
├── qiime2_artifacts/
│   ├── table.qza                        # DADA2 feature table (ASV level)
│   ├── rep-seqs.qza                     # representative ASV sequences
│   └── taxonomy.qza                     # taxonomy assignments
│
├── scripts/
│   ├── 01_richness_glm.R                # GLM richness ~ topology + age; IQR outlier removal
│   ├── 02_anova_richness_topology.R     # ANOVA with age (full model — kept for transparency)
│   ├── 03_anova_richness_no_age.R       # ANOVA without age; Tukey CLD; richness plots
│   ├── 04_wilcoxon_young_old.R          # Udder and Dorsum richness: young vs old
│   ├── 05_family_abundance_tests.R      # Kruskal-Wallis + Wilcoxon per family, BH-corrected
│   ├── 06_family_abundance_figure.R     # 4-panel abundance figure
│   ├── 07_rf_topology_predictors.R      # RF multiclass A/B/C; MDA importance
│   ├── 08_rf_validation.R               # Overdispersion, multinomial, ROC, PCA, Kraken2 overlap
│   ├── 09_rf_diagnostics.R              # Convergence, seed stability, calibration, Youden's J
│   ├── 10_roc_smoothed.R                # Smoothed ROC with bootstrap CI (binormal, pROC)
│   ├── 11_analysis_ACgroups.R           # Binary RF AB vs C; all _onlyACgroups outputs
│   ├── 12_picrust2_install_run.sh       # PICRUSt2 installation + pipeline run
│   └── 13_picrust2_analysis.R           # Functional prediction analysis (post-PICRUSt2)
│
├── probe_design/
│   └── diagnostic_probe_set.md          # Four family-specific 16S V3-V4 probe sequences
│
└── results/
    ├── figures/                         # publication-ready figures (PDF + PNG)
    └── tables/                          # result tables (CSV)
```

> **Raw reads** (demux-paired.qza, ~1.4 GB) are deposited at NCBI SRA under
> accession **[TBD upon submission]** and are excluded from this repository.

---

## How to reproduce

### Prerequisites

**R packages** (install once):

```r
install.packages(c("tidyverse", "randomForest", "pROC", "nnet",
                   "patchwork", "multcomp", "MASS"))
```

**PICRUSt2** (conda, ~2 GB, first run only):

```bash
bash scripts/12_picrust2_install_run.sh
```

### Run order

Scripts are numbered to reflect execution order.
Each script reads from `data/` and writes to `results/`.

```bash
cd /path/to/repo

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

# After PICRUSt2 pipeline completes (~30–90 min):
bash scripts/12_picrust2_install_run.sh
Rscript scripts/13_picrust2_analysis.R
```

All scripts use `set.seed(42)` for reproducibility.

---

## Diagnostic probe set

Four taxon-specific 20–22-mer probes targeting the V3 sub-region of the 16S
V3-V4 amplicon (341F–805R) identify the binary contact-zone classification
with zero cross-reactivity among the four target families (within-study evaluation):

| ID | Family | Sequence 5'→3' | Tm (°C) | GC% | Sensitivity |
|---|---|---|---|---|---|
| DEIN-V3 | *Deinococcaceae* | `CAGCCGCGGTAATACGGAGG` | 52.8 | 65 | 97% |
| INTR-V3 | *Intrasporangiaceae* | `TATTGGGCGTAAAGAGCTTG` | 44.6 | 45 | 100% |
| HUNG-V3 | *Hungateiclostridiaceae* | `TTACTGGGTGTAAAGGGCGTGT` | 49.7 | 50 | 82% |
| CHRI-V3 | *Christensenellaceae* | `AGGAAGCCCCGGCTAACTACGT` | 53.4 | 59 | 90% |

DEIN-V3 positive / INTR-V3 HUNG-V3 CHRI-V3 negative → **EtS** site.
Reverse pattern → **GCtS** site.

Full derivation details: [`probe_design/diagnostic_probe_set.md`](probe_design/diagnostic_probe_set.md)

---

## Statistical methods summary

Family richness was modelled with a quasi-Poisson GLM per topology site
(IQR-based outlier exclusion per site), with topology and age as predictors.
Type II marginal F-tests (`drop1`) were used throughout; age was non-significant
(p = 0.46–0.56) and excluded from the simplified model. Tukey HSD post-hoc
comparisons used compact letter display (`multcomp::cld`). Family abundance
differences were tested with Kruskal-Wallis (topology, 4-level including
environment) and Wilcoxon (AB vs C, age binary), with Benjamini–Hochberg FDR
correction applied across all tests simultaneously. Random Forest classifiers
(1000 trees) used out-of-bag error for internal validation; stability was
confirmed across 10 random seeds (OOB SD = 0.73%). ROC curves used binormal
smoothing (`pROC::smooth`); AUC confidence intervals were obtained by bootstrap
resampling (n = 2000). Youden's J and its bootstrap CI were computed by
manual resampling to avoid boundary failures in `pROC::ci.coords`.
Functional predictions used PICRUSt2 v2.6.3 (stratified MetaCyc pathways,
EC numbers, KO terms). Probe sequences were designed by MAFFT alignment
(v7.505) of within-family ASVs and a sliding-window conserved-region search.

---

## Comparison with published studies

| | This study | Strompfová 2024 | O'Shaughnessy-Hunter 2021 | Styková 2025 |
|---|---|---|---|---|
| Breed | Norik Muránsky | Shetland pony | Mixed | Norik |
| Type | Cold-blooded | Cold-blooded | Mixed | Cold-blooded |
| n horses | 12 | 6 | 12 | 18 |
| Sites | 8 + env | 5 | 4 | 1 (pastern) |
| 16S region | **V3-V4** | V3-V4 | V4 | V4 |
| Primary driver | Contact zone | Individual | Season | Disease |
| Age effect | None (5–17 yr) | — | — | — |
| *Deinococcaceae* | ✓ (EtS) | ✓ (back) | ✗ | ✗ |
| Probe set proposed | ✓ | ✗ | ✗ | ✗ |

*Deinococcaceae* on dorsal/elevated sites is confirmed across two independent
cold-blooded equine studies; V4-only protocols appear to miss this taxon.

---

## Data availability

- Raw reads: NCBI SRA **[accession TBD]**
- Processed data: this repository (`data/` directory)
- PICRUSt2 outputs: regenerate using `scripts/12_picrust2_install_run.sh`

---

## Citation

> Kováčová V et al. (2026) Body-site contact zone determines skin microbiome
> composition in healthy cold-blooded draft horses, with identification of a
> four-family diagnostic marker panel. *Frontiers in Veterinary Science* [in preparation].

---

## License

Code: MIT License.
Data: CC BY 4.0 (upon publication).
