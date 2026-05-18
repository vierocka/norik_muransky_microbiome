#!/usr/bin/env Rscript
# 14_lmm_sensitivity.R — Sensitivity & robustness analysis
#
# Three checks beyond the main analysis:
#   1. ICC  — how much richness/abundance variance is between horses?
#   2. LMM  — linear mixed model with horse as random effect
#             (replaces quasi-Poisson GLM and Wilcoxon for formal test)
#   3. Permutation — within-horse label permutation to validate H0
#   4. Noise injection — 1 / 5 / 10 % multiplicative noise on RF stability

library(lme4)
library(lmerTest)   # Satterthwaite df for lmer p-values
library(tidyverse)
library(randomForest)

set.seed(42)
N_PERM  <- 999    # permutation replicates
N_SEEDS <- 10     # seeds for noise RF

# =============================================================================
# Data
# =============================================================================
contrib <- read_csv("data/pca_family_contributions.csv", show_col_types = FALSE)
meta    <- read_csv("data/metadata.csv", show_col_types = FALSE)

fam_cols <- setdiff(colnames(contrib), c("sample_id", "topology"))

dat <- contrib %>%
  left_join(select(meta, sample_id, age), by = "sample_id") %>%
  filter(topology != "Environment") %>%
  mutate(
    horse_id = sub("-.*", "", sample_id) %>% factor(),
    richness  = rowSums(across(all_of(fam_cols), ~ . > 0)),
    topo3 = case_when(
      topology %in% c("Left front pastern", "Muzzle")              ~ "GCtS_A",
      topology %in% c("Ventral abdomen",    "Udder")               ~ "GCtS_B",
      topology %in% c("Dorsum", "Forehead", "Neck", "Pectoral area") ~ "EtS"
    ) %>% factor(levels = c("GCtS_A", "GCtS_B", "EtS")),
    group2 = ifelse(topo3 == "EtS", "EtS", "GCtS") %>%
             factor(levels = c("GCtS", "EtS"))
  )

cat("=== Dataset ===\n")
cat("Horses:", n_distinct(dat$horse_id), " | Samples:", nrow(dat), "\n")
cat("Sites per horse:", nrow(dat) / n_distinct(dat$horse_id), "\n\n")

# =============================================================================
# 1. ICC — intraclass correlation for horse-level variance
# =============================================================================
cat("=== 1. ICC: between-horse variance ===\n")

icc_rich <- lmer(richness ~ 1 + (1|horse_id), data = dat, REML = TRUE)
vc <- as.data.frame(VarCorr(icc_rich))
icc_val  <- vc$vcov[1] / sum(vc$vcov)
cat(sprintf("Richness ICC (horse):  %.3f  (between=%.2f, within=%.2f)\n",
            icc_val, vc$vcov[1], vc$vcov[2]))

# ICC per family (mean across all 46)
eps <- 1e-5
icc_fam <- map_dbl(fam_cols, function(f) {
  d <- dat %>% mutate(lab = log(.data[[f]] + eps))
  m <- tryCatch(lmer(lab ~ 1 + (1|horse_id), data = d, REML = TRUE), error = function(e) NULL)
  if (is.null(m)) return(NA_real_)
  vc2 <- as.data.frame(VarCorr(m))
  vc2$vcov[1] / sum(vc2$vcov)
})
cat(sprintf("Family abundance ICC:  mean=%.3f  median=%.3f  range [%.3f, %.3f]\n\n",
            mean(icc_fam, na.rm=TRUE), median(icc_fam, na.rm=TRUE),
            min(icc_fam, na.rm=TRUE), max(icc_fam, na.rm=TRUE)))

icc_tbl <- tibble(family = fam_cols, icc = icc_fam) %>%
  bind_rows(tibble(family = "RICHNESS", icc = icc_val)) %>%
  arrange(desc(icc))
write.csv(icc_tbl, "lmm_icc.csv", row.names = FALSE)
cat("Saved lmm_icc.csv\n")

# =============================================================================
# 2a. LMM richness — glmer(Poisson) + (1|horse_id)
# =============================================================================
cat("\n=== 2a. LMM richness ~ topo3 + age + (1|horse_id) ===\n")
cat("(Compare to script 03 quasi-Poisson GLM)\n\n")

fit_rich <- glmer(richness ~ topo3 + age + (1|horse_id),
                  data   = dat,
                  family = poisson,
                  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

cat("Fixed effects (log-scale):\n")
coef_rich <- coef(summary(fit_rich))
print(round(coef_rich, 4))

cat("\nType II likelihood-ratio tests (drop1):\n")
drop_rich <- drop1(fit_rich, test = "Chisq")
print(drop_rich)

# IRR (incidence rate ratios) for topology
irr <- exp(fixef(fit_rich))
cat("\nIRR (exp(coef)) — topo3 relative to GCtS_A:\n")
print(round(irr, 3))

lmm_rich_out <- tibble(
  term     = names(fixef(fit_rich)),
  estimate = fixef(fit_rich),
  IRR      = exp(fixef(fit_rich)),
  z_value  = coef_rich[, "z value"],
  p_value  = coef_rich[, "Pr(>|z|)"]
)
write.csv(lmm_rich_out, "lmm_richness_glmm.csv", row.names = FALSE)
cat("Saved lmm_richness_glmm.csv\n")

# =============================================================================
# 2b. LMM family abundance ~ group2 + (1|horse_id) for all 46 families
# =============================================================================
cat("\n=== 2b. LMM per-family: GCtS vs EtS + (1|horse_id) ===\n")
cat("(Compare to script 11 Wilcoxon)\n\n")

run_lmm_fam <- function(f, data) {
  d <- data %>% mutate(log_ab = log(.data[[f]] + eps))
  fit <- tryCatch(
    suppressMessages(
      lmer(log_ab ~ group2 + (1|horse_id), data = d,
           control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(tibble(family=f, beta_EtS=NA, SE=NA, t=NA, df=NA, p_lmm=NA, converged=FALSE))
  ct <- coef(summary(fit))
  if (!"group2EtS" %in% rownames(ct)) return(tibble(family=f, beta_EtS=NA, SE=NA, t=NA, df=NA, p_lmm=NA, converged=FALSE))
  tibble(
    family    = f,
    beta_EtS  = ct["group2EtS", "Estimate"],
    SE        = ct["group2EtS", "Std. Error"],
    t         = ct["group2EtS", "t value"],
    df        = ct["group2EtS", "df"],
    p_lmm     = ct["group2EtS", "Pr(>|t|)"],
    converged = !any(grepl("Model failed", fit@optinfo$conv$lme4$messages %||% ""))
  )
}

lmm_fam <- map_dfr(fam_cols, run_lmm_fam, data = dat) %>%
  mutate(p_adj_lmm = p.adjust(p_lmm, method = "BH"),
         sig_lmm   = p_adj_lmm < 0.05)

# Load original Wilcoxon for comparison
wilc <- read_csv("famabund_wilcoxon_ACvsC_onlyACgroups.csv", show_col_types = FALSE) %>%
  select(family, log2FC, p_adj_wilcox = p_adj)

comparison <- lmm_fam %>%
  left_join(wilc, by = "family") %>%
  mutate(sig_wilcox = p_adj_wilcox < 0.05,
         agree      = sig_lmm == sig_wilcox) %>%
  arrange(p_adj_lmm)

write.csv(comparison, "lmm_vs_wilcoxon_comparison.csv", row.names = FALSE)

cat(sprintf("Significant BH<0.05 — LMM: %d  |  Wilcoxon: %d\n",
            sum(comparison$sig_lmm, na.rm=TRUE),
            sum(comparison$sig_wilcox == "***" | comparison$sig_wilcox == "**" | comparison$sig_wilcox == "*", na.rm=TRUE)))
cat(sprintf("Agreement: %d / %d families\n", sum(comparison$agree, na.rm=TRUE), nrow(comparison)))

discordant <- comparison %>% filter(!agree, !is.na(agree))
if (nrow(discordant) > 0) {
  cat("\nDiscordant families:\n")
  print(discordant %>% select(family, p_adj_lmm, p_adj_wilcox, sig_lmm))
} else {
  cat("No discordant families — LMM and Wilcoxon fully agree.\n")
}

cat("\nTop 15 by LMM:\n")
print(comparison %>% select(family, beta_EtS, p_adj_lmm, log2FC, p_adj_wilcox) %>% head(15), n=15)
cat("Saved lmm_vs_wilcoxon_comparison.csv\n")

# =============================================================================
# 3. Within-horse permutation test
# Permute topology labels WITHIN each horse to preserve blocking structure.
# Tests richness difference (EtS vs GCtS) and each significant family.
# =============================================================================
cat("\n=== 3. Within-horse permutation (n =", N_PERM, "replicates) ===\n")

# Observed group means
obs_rich_diff <- mean(dat$richness[dat$group2 == "GCtS"]) -
                 mean(dat$richness[dat$group2 == "EtS"])

sig_fams <- comparison %>% filter(sig_lmm) %>% pull(family)

perm_stat <- function(data, fam = NULL) {
  perm_dat <- data %>%
    group_by(horse_id) %>%
    mutate(group2_perm = sample(group2)) %>%
    ungroup()
  if (is.null(fam)) {
    mean(perm_dat$richness[perm_dat$group2_perm == "GCtS"]) -
    mean(perm_dat$richness[perm_dat$group2_perm == "EtS"])
  } else {
    mean(perm_dat[[fam]][perm_dat$group2_perm == "GCtS"]) -
    mean(perm_dat[[fam]][perm_dat$group2_perm == "EtS"])
  }
}

set.seed(42)
null_rich <- replicate(N_PERM, perm_stat(dat))
p_perm_rich <- mean(abs(null_rich) >= abs(obs_rich_diff))
cat(sprintf("Richness: obs diff=%.2f, permutation p=%.4f\n", obs_rich_diff, p_perm_rich))

perm_fam_out <- map_dfr(sig_fams, function(f) {
  obs <- mean(dat[[f]][dat$group2 == "GCtS"]) - mean(dat[[f]][dat$group2 == "EtS"])
  null <- replicate(N_PERM, perm_stat(dat, f))
  tibble(family = f, obs_diff = obs, p_perm = mean(abs(null) >= abs(obs)))
})

perm_fam_out <- perm_fam_out %>%
  mutate(p_adj_perm = p.adjust(p_perm, method = "BH")) %>%
  arrange(p_perm)

write.csv(perm_fam_out, "permutation_results.csv", row.names = FALSE)

cat(sprintf("Families significant by permutation (BH<0.05): %d / %d\n",
            sum(perm_fam_out$p_adj_perm < 0.05), nrow(perm_fam_out)))
cat("Saved permutation_results.csv\n")
cat("\nPermutation results (significant families):\n")
print(perm_fam_out %>% filter(p_adj_perm < 0.05) %>% select(family, obs_diff, p_perm, p_adj_perm), n=30)

# =============================================================================
# 4. Noise injection — RF OOB stability under added noise
# Add multiplicative Gaussian noise at 1%, 5%, 10% level
# =============================================================================
cat("\n=== 4. Noise injection: RF stability ===\n")

X_clean <- dat %>% select(all_of(fam_cols)) %>% as.matrix()
y       <- dat$group2

add_noise <- function(X, pct) {
  noise <- matrix(rnorm(prod(dim(X)), mean = 0, sd = pct), nrow = nrow(X))
  pmax(X + X * noise, 0)   # keep non-negative
}

noise_levels <- c(0, 0.01, 0.05, 0.10)

noise_results <- map_dfr(noise_levels, function(pct) {
  oob_vec <- map_dbl(seq_len(N_SEEDS), function(s) {
    set.seed(s)
    X_noisy <- if (pct == 0) X_clean else add_noise(X_clean, pct)
    rf <- randomForest(x = X_noisy, y = y, ntree = 500)
    1 - rf$err.rate[500, "OOB"]  # accuracy
  })
  # Top predictors
  set.seed(42)
  X_noisy <- if (pct == 0) X_clean else add_noise(X_clean, pct)
  rf_rep  <- randomForest(x = X_noisy, y = y, ntree = 1000, importance = TRUE)
  top5    <- names(sort(importance(rf_rep)[,"MeanDecreaseAccuracy"], dec=TRUE))[1:5]
  tibble(
    noise_pct    = pct * 100,
    mean_acc     = mean(oob_vec) * 100,
    sd_acc       = sd(oob_vec) * 100,
    top5_predictors = paste(top5, collapse = " | ")
  )
})

write.csv(noise_results, "noise_injection_rf.csv", row.names = FALSE)

cat("\nRF accuracy (OOB) under added noise:\n")
print(noise_results %>% select(noise_pct, mean_acc, sd_acc, top5_predictors))
cat("Saved noise_injection_rf.csv\n")

# =============================================================================
# Summary
# =============================================================================
cat("\n=== SUMMARY ===\n")
cat(sprintf("Horse ICC (richness):         %.3f\n", icc_val))
cat(sprintf("Horse ICC (family abundance): mean %.3f, median %.3f\n",
            mean(icc_fam, na.rm=TRUE), median(icc_fam, na.rm=TRUE)))
cat(sprintf("LMM significant families:    %d / 46\n", sum(lmm_fam$sig_lmm, na.rm=TRUE)))
cat(sprintf("LMM-Wilcoxon agreement:      %d / %d\n", sum(comparison$agree, na.rm=TRUE), nrow(comparison)))
cat(sprintf("Permutation richness p:       %.4f\n", p_perm_rich))
cat(sprintf("RF accuracy (0%% noise):       %.1f%% ± %.1f%%\n",
            noise_results$mean_acc[1], noise_results$sd_acc[1]))
cat(sprintf("RF accuracy (10%% noise):      %.1f%% ± %.1f%%\n",
            noise_results$mean_acc[4], noise_results$sd_acc[4]))
