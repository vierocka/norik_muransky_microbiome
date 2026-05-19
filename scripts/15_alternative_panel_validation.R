#!/usr/bin/env Rscript
# 15_alternative_panel_validation.R — Alternative generalised 4-family panel
#
# Panel A (original):   Deinococcaceae + Intrasporangiaceae + Hungateiclostridiaceae + Christensenellaceae
# Panel B (generalised): Deinococcaceae + Micrococcaceae    + Lachnospiraceae        + UCG-010
#
# Analyses for each panel (and head-to-head comparison):
#   1. Random Forest binary classifier  (GCtS vs EtS, OOB, N_SEEDS replicates)
#   2. ROC curve + AUC
#   3. Sensitivity / Specificity / Youden's J at Youden-optimal threshold (bootstrap CI)
#   4. Within-horse permutation test  (N_PERM, restricted permutation)
#   5. Noise injection stability       (1 / 2 / 5 / 10 % multiplicative Gaussian)

library(tidyverse)
library(randomForest)
library(pROC)

set.seed(42)
N_SEEDS  <- 50    # RF replicates for mean ± SD OOB
N_PERM   <- 999   # permutation replicates
N_BOOT   <- 2000  # bootstrap replicates for CI on ROC metrics
NOISE_LEVELS <- c(0.01, 0.02, 0.05, 0.10)  # 1, 2, 5, 10 %

PANEL_A <- c("Deinococcaceae", "Intrasporangiaceae", "Hungateiclostridiaceae", "Christensenellaceae")
PANEL_B <- c("Deinococcaceae", "Micrococcaceae",     "Lachnospiraceae",        "UCG-010")

OUT_TABLES  <- "results/tables"
OUT_FIGURES <- "results/figures"

# =============================================================================
# Data
# =============================================================================
contrib <- read_csv("data/pca_family_contributions.csv", show_col_types = FALSE)
meta    <- read_csv("data/metadata.csv",                 show_col_types = FALSE)

fam_cols <- setdiff(colnames(contrib), c("sample_id", "topology"))

dat <- contrib %>%
  left_join(select(meta, sample_id, age), by = "sample_id") %>%
  filter(topology != "Environment") %>%
  mutate(
    horse_id = sub("-.*", "", sample_id) %>% factor(),
    topo3 = case_when(
      topology %in% c("Left front pastern", "Muzzle")                    ~ "GCtS_A",
      topology %in% c("Ventral abdomen",    "Udder")                      ~ "GCtS_B",
      topology %in% c("Dorsum", "Forehead", "Neck", "Pectoral area")     ~ "EtS"
    ) %>% factor(levels = c("GCtS_A", "GCtS_B", "EtS")),
    group2 = ifelse(topo3 == "EtS", "EtS", "GCtS") %>%
             factor(levels = c("GCtS", "EtS"))
  )

cat("=== Dataset ===\n")
cat("Horses:", n_distinct(dat$horse_id), " | Samples:", nrow(dat), "\n")
cat("GCtS:", sum(dat$group2 == "GCtS"), " | EtS:", sum(dat$group2 == "EtS"), "\n\n")

# =============================================================================
# Helpers
# =============================================================================

add_noise <- function(X, pct) {
  noise <- matrix(rnorm(prod(dim(X)), mean = 0, sd = pct), nrow = nrow(X))
  pmax(X + X * noise, 0)
}

# Run N_SEEDS RF replicates; return mean OOB accuracy ± SD
rf_oob_replicated <- function(X, y, n_seeds = N_SEEDS, ntree = 500) {
  acc_vec <- map_dbl(seq_len(n_seeds), function(s) {
    set.seed(s)
    rf <- randomForest(x = X, y = y, ntree = ntree)
    1 - rf$err.rate[ntree, "OOB"]
  })
  list(mean = mean(acc_vec) * 100, sd = sd(acc_vec) * 100, vec = acc_vec * 100)
}

# Build a single RF (seed 42, ntree 1000) and return OOB vote probabilities
rf_oob_probs <- function(X, y, ntree = 1000) {
  set.seed(42)
  rf <- randomForest(x = X, y = y, ntree = ntree,
                     keep.inbag = TRUE)
  # OOB vote matrix: rows = samples, cols = classes
  votes <- rf$votes          # proportion votes for each class
  list(rf = rf, votes = votes)
}

# ROC + Youden metrics with bootstrap CIs
roc_youden_boot <- function(labels, probs_pos, n_boot = N_BOOT, pos_class = "EtS") {
  roc_obj <- roc(labels, probs_pos, levels = c("GCtS", "EtS"),
                 direction = "<", quiet = TRUE)
  auc_obs  <- as.numeric(auc(roc_obj))

  # Youden optimal threshold
  youden_idx <- which.max(roc_obj$sensitivities + roc_obj$specificities - 1)
  thr_obs    <- roc_obj$thresholds[youden_idx]
  sens_obs   <- roc_obj$sensitivities[youden_idx]
  spec_obs   <- roc_obj$specificities[youden_idx]
  j_obs      <- sens_obs + spec_obs - 1

  # Bootstrap
  n <- length(labels)
  boot_mat <- replicate(n_boot, {
    idx <- sample(n, n, replace = TRUE)
    lb  <- labels[idx]
    pr  <- probs_pos[idx]
    if (length(unique(lb)) < 2) return(rep(NA_real_, 4))
    r   <- tryCatch(roc(lb, pr, levels = c("GCtS", "EtS"),
                        direction = "<", quiet = TRUE), error = function(e) NULL)
    if (is.null(r)) return(rep(NA_real_, 4))
    yi  <- which.max(r$sensitivities + r$specificities - 1)
    c(auc  = as.numeric(auc(r)),
      sens = r$sensitivities[yi],
      spec = r$specificities[yi],
      j    = r$sensitivities[yi] + r$specificities[yi] - 1)
  })

  ci95 <- function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE)

  list(
    auc      = auc_obs,
    auc_ci   = ci95(boot_mat["auc", ]),
    sens     = sens_obs,
    sens_ci  = ci95(boot_mat["sens", ]),
    spec     = spec_obs,
    spec_ci  = ci95(boot_mat["spec", ]),
    j        = j_obs,
    j_ci     = ci95(boot_mat["j", ]),
    thr      = thr_obs,
    roc_obj  = roc_obj
  )
}

# Within-horse permutation p-value for OOB accuracy
perm_oob_acc <- function(X, y, horse_id, n_perm = N_PERM, ntree = 500) {
  set.seed(42)
  rf_obs  <- randomForest(x = X, y = y, ntree = ntree)
  obs_acc <- 1 - rf_obs$err.rate[ntree, "OOB"]

  horses <- levels(horse_id)
  null_acc <- replicate(n_perm, {
    perm_y <- y
    for (h in horses) {
      idx <- which(horse_id == h)
      perm_y[idx] <- sample(y[idx])
    }
    if (length(unique(perm_y)) < 2) return(NA_real_)
    rf_p <- randomForest(x = X, y = perm_y, ntree = ntree)
    1 - rf_p$err.rate[ntree, "OOB"]
  })

  null_acc <- null_acc[!is.na(null_acc)]
  p_val    <- mean(null_acc >= obs_acc)
  list(obs_acc = obs_acc * 100, p_perm = p_val, null_dist = null_acc * 100)
}

# Noise stability across NOISE_LEVELS + baseline (0)
noise_stability <- function(X_clean, y, noise_levels = NOISE_LEVELS, n_seeds = N_SEEDS) {
  all_levels <- c(0, noise_levels)
  map_dfr(all_levels, function(pct) {
    acc_vec <- map_dbl(seq_len(n_seeds), function(s) {
      set.seed(s)
      X_n <- if (pct == 0) X_clean else add_noise(X_clean, pct)
      rf  <- randomForest(x = X_n, y = y, ntree = 500)
      1 - rf$err.rate[500, "OOB"]
    })
    # top predictors at seed 42
    set.seed(42)
    X_n <- if (pct == 0) X_clean else add_noise(X_clean, pct)
    rf_rep <- randomForest(x = X_n, y = y, ntree = 1000, importance = TRUE)
    top_k <- min(ncol(X_n), 4L)
    top4   <- names(sort(importance(rf_rep)[, "MeanDecreaseAccuracy"],
                         decreasing = TRUE))[seq_len(top_k)]
    tibble(
      noise_pct = pct * 100,
      mean_acc  = mean(acc_vec) * 100,
      sd_acc    = sd(acc_vec)   * 100,
      top_predictors = paste(top4, collapse = " | ")
    )
  })
}

# =============================================================================
# Run both panels
# =============================================================================

panels <- list(
  "Panel_A_original"    = PANEL_A,
  "Panel_B_generalised" = PANEL_B
)

all_roc_data   <- list()
all_metrics    <- list()
all_noise      <- list()
all_perm       <- list()

for (pname in names(panels)) {
  fams <- panels[[pname]]
  cat(sprintf("\n=== %s ===\n", pname))
  cat(sprintf("Families: %s\n\n", paste(fams, collapse = ", ")))

  X <- dat %>% select(all_of(fams)) %>% as.matrix()
  y <- dat$group2
  h <- dat$horse_id

  # --- 1. OOB accuracy (replicated) ---
  oob <- rf_oob_replicated(X, y)
  cat(sprintf("OOB accuracy: %.1f%% ± %.1f%%\n", oob$mean, oob$sd))

  # --- 2+3. ROC + Youden ---
  rv <- rf_oob_probs(X, y)
  probs_ets <- rv$votes[, "EtS"]
  metrics   <- roc_youden_boot(y, probs_ets)

  cat(sprintf("AUC:      %.4f  [%.4f, %.4f]\n",
              metrics$auc, metrics$auc_ci[1], metrics$auc_ci[2]))
  cat(sprintf("Youden J: %.4f  [%.4f, %.4f]\n",
              metrics$j, metrics$j_ci[1], metrics$j_ci[2]))
  cat(sprintf("Sens:     %.4f  [%.4f, %.4f]\n",
              metrics$sens, metrics$sens_ci[1], metrics$sens_ci[2]))
  cat(sprintf("Spec:     %.4f  [%.4f, %.4f]\n",
              metrics$spec, metrics$spec_ci[1], metrics$spec_ci[2]))
  cat(sprintf("Optimal threshold: %.4f\n", metrics$thr))

  roc_df <- data.frame(
    panel       = pname,
    fpr         = 1 - metrics$roc_obj$specificities,
    tpr         = metrics$roc_obj$sensitivities,
    threshold   = metrics$roc_obj$thresholds
  )
  all_roc_data[[pname]] <- roc_df

  all_metrics[[pname]] <- tibble(
    panel     = pname,
    families  = paste(fams, collapse = "; "),
    mean_oob_acc = oob$mean,
    sd_oob_acc   = oob$sd,
    auc       = metrics$auc,
    auc_lo    = metrics$auc_ci[1],
    auc_hi    = metrics$auc_ci[2],
    youden_j  = metrics$j,
    j_lo      = metrics$j_ci[1],
    j_hi      = metrics$j_ci[2],
    sens      = metrics$sens,
    sens_lo   = metrics$sens_ci[1],
    sens_hi   = metrics$sens_ci[2],
    spec      = metrics$spec,
    spec_lo   = metrics$spec_ci[1],
    spec_hi   = metrics$spec_ci[2],
    opt_thr   = metrics$thr
  )

  # --- 4. Within-horse permutation ---
  cat(sprintf("\nPermutation test (N=%d)...\n", N_PERM))
  perm_res <- perm_oob_acc(X, y, h)
  cat(sprintf("Observed OOB: %.1f%%, permutation p = %.4f\n",
              perm_res$obs_acc, perm_res$p_perm))

  all_perm[[pname]] <- tibble(
    panel       = pname,
    obs_acc     = perm_res$obs_acc,
    p_perm      = perm_res$p_perm
  )

  # --- 5. Noise stability ---
  cat("\nNoise injection stability...\n")
  noise_df <- noise_stability(X, y)
  noise_df$panel <- pname
  all_noise[[pname]] <- noise_df
  print(noise_df %>% select(noise_pct, mean_acc, sd_acc, top_predictors))
}

# =============================================================================
# Save tables
# =============================================================================

metrics_tbl <- bind_rows(all_metrics)
noise_tbl   <- bind_rows(all_noise)
perm_tbl    <- bind_rows(all_perm)
roc_tbl     <- bind_rows(all_roc_data)

write.csv(metrics_tbl, file.path(OUT_TABLES, "panel_validation_metrics.csv"),     row.names = FALSE)
write.csv(noise_tbl,   file.path(OUT_TABLES, "panel_noise_stability.csv"),        row.names = FALSE)
write.csv(perm_tbl,    file.path(OUT_TABLES, "panel_permutation_test.csv"),       row.names = FALSE)
write.csv(roc_tbl,     file.path(OUT_TABLES, "panel_roc_curves.csv"),             row.names = FALSE)

cat("\nSaved panel_validation_metrics.csv\n")
cat("Saved panel_noise_stability.csv\n")
cat("Saved panel_permutation_test.csv\n")
cat("Saved panel_roc_curves.csv\n")

# =============================================================================
# Figure 1 — ROC curves (both panels, one plot)
# =============================================================================

panel_cols <- c(
  "Panel_A_original"    = "#1f77b4",
  "Panel_B_generalised" = "#d62728"
)
panel_labels <- c(
  "Panel_A_original"    = "Panel A: DEIN + INTR + HUNG + CHRI",
  "Panel_B_generalised" = "Panel B: DEIN + MICC + LACHN + UCG-010"
)

auc_labels <- metrics_tbl %>%
  mutate(label = sprintf("%s\nAUC = %.3f [%.3f–%.3f]",
                         panel_labels[panel], auc, auc_lo, auc_hi))

p_roc <- ggplot(roc_tbl, aes(x = fpr, y = tpr, colour = panel)) +
  geom_line(linewidth = 1.1) +
  geom_abline(linetype = "dashed", colour = "grey60") +
  annotate("text", x = 0.55, y = 0.25,
           label = auc_labels$label[auc_labels$panel == "Panel_A_original"],
           colour = panel_cols["Panel_A_original"], size = 3.5, hjust = 0) +
  annotate("text", x = 0.55, y = 0.08,
           label = auc_labels$label[auc_labels$panel == "Panel_B_generalised"],
           colour = panel_cols["Panel_B_generalised"], size = 3.5, hjust = 0) +
  scale_colour_manual(values = panel_cols, labels = panel_labels) +
  coord_equal() +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom",
        legend.title    = element_blank()) +
  labs(title    = "ROC curves: 4-family panel classifiers (GCtS vs EtS)",
       subtitle = "OOB vote probabilities, N=12 horses",
       x = "1 – Specificity (FPR)", y = "Sensitivity (TPR)")

ggsave(file.path(OUT_FIGURES, "panel_roc_curves.pdf"), p_roc, width = 6, height = 6)
ggsave(file.path(OUT_FIGURES, "panel_roc_curves.png"), p_roc, width = 6, height = 6, dpi = 300)
cat("Saved panel_roc_curves.pdf/.png\n")

# =============================================================================
# Figure 2 — Noise stability bar chart
# =============================================================================

noise_plot_df <- noise_tbl %>%
  mutate(noise_label = paste0(noise_pct, "%"),
         noise_label = factor(noise_label, levels = paste0(c(0, 1, 2, 5, 10), "%")),
         panel_label = panel_labels[panel])

p_noise <- ggplot(noise_plot_df,
                  aes(x = noise_label, y = mean_acc, fill = panel,
                      group = panel)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
  geom_errorbar(aes(ymin = mean_acc - sd_acc, ymax = mean_acc + sd_acc),
                position = position_dodge(width = 0.7), width = 0.25) +
  scale_fill_manual(values = panel_cols, labels = panel_labels) +
  scale_y_continuous(limits = c(80, 100), oob = scales::oob_keep) +
  theme_classic(base_size = 12) +
  theme(legend.position  = "bottom",
        legend.title     = element_blank()) +
  labs(title    = "RF OOB accuracy under multiplicative noise",
       subtitle  = sprintf("Mean ± SD over %d random seeds, ntree=500", N_SEEDS),
       x = "Added noise level", y = "OOB accuracy (%)")

ggsave(file.path(OUT_FIGURES, "panel_noise_stability.pdf"), p_noise, width = 7, height = 5)
ggsave(file.path(OUT_FIGURES, "panel_noise_stability.png"), p_noise, width = 7, height = 5, dpi = 300)
cat("Saved panel_noise_stability.pdf/.png\n")

# =============================================================================
# Figure 3 — Head-to-head metrics comparison (forest-plot style)
# =============================================================================

metric_long <- metrics_tbl %>%
  transmute(
    panel,
    panel_label = panel_labels[panel],
    AUC           = auc,      AUC_lo      = auc_lo,   AUC_hi      = auc_hi,
    `Youden J`    = youden_j, `Youden J_lo` = j_lo,   `Youden J_hi` = j_hi,
    Sensitivity   = sens,     Sensitivity_lo = sens_lo, Sensitivity_hi = sens_hi,
    Specificity   = spec,     Specificity_lo = spec_lo, Specificity_hi = spec_hi
  ) %>%
  pivot_longer(cols = -c(panel, panel_label), names_to = "var", values_to = "val") %>%
  mutate(
    bound  = case_when(
      str_ends(var, "_lo") ~ "lo",
      str_ends(var, "_hi") ~ "hi",
      TRUE                 ~ "est"
    ),
    metric = str_remove(var, "_(lo|hi)$")
  ) %>%
  select(-var) %>%
  pivot_wider(id_cols = c(panel, panel_label, metric),
              names_from = bound, values_from = val) %>%
  mutate(metric = factor(metric,
                         levels = c("AUC", "Youden J", "Sensitivity", "Specificity")))

p_forest <- ggplot(metric_long,
                   aes(x = est, y = metric, colour = panel)) +
  geom_linerange(aes(xmin = lo, xmax = hi),
                 position = position_dodge(width = 0.5), linewidth = 0.7) +
  geom_point(aes(shape = panel), size = 3,
             position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = panel_cols, labels = panel_labels) +
  scale_shape_manual(values = c(16, 17), labels = panel_labels) +
  scale_x_continuous(limits = c(0.85, 1.005), breaks = seq(0.85, 1, 0.05)) +
  theme_classic(base_size = 12) +
  theme(legend.position  = "bottom",
        legend.title     = element_blank(),
        axis.title.y     = element_blank()) +
  labs(title    = "Panel classifier performance: head-to-head comparison",
       subtitle  = sprintf("Point estimates with 95%% bootstrap CI (B=%d)", N_BOOT),
       x = "Metric value")

ggsave(file.path(OUT_FIGURES, "panel_metrics_comparison.pdf"), p_forest, width = 7, height = 4)
ggsave(file.path(OUT_FIGURES, "panel_metrics_comparison.png"), p_forest, width = 7, height = 4, dpi = 300)
cat("Saved panel_metrics_comparison.pdf/.png\n")

# =============================================================================
# Summary
# =============================================================================
cat("\n=== SUMMARY ===\n")
for (i in seq_len(nrow(metrics_tbl))) {
  r <- metrics_tbl[i, ]
  cat(sprintf("\n%s\n", r$panel))
  cat(sprintf("  OOB accuracy:  %.1f%% ± %.1f%%\n", r$mean_oob_acc, r$sd_oob_acc))
  cat(sprintf("  AUC:           %.4f [%.4f, %.4f]\n", r$auc, r$auc_lo, r$auc_hi))
  cat(sprintf("  Youden J:      %.4f [%.4f, %.4f]\n", r$youden_j, r$j_lo, r$j_hi))
  cat(sprintf("  Sensitivity:   %.4f [%.4f, %.4f]\n", r$sens, r$sens_lo, r$sens_hi))
  cat(sprintf("  Specificity:   %.4f [%.4f, %.4f]\n", r$spec, r$spec_lo, r$spec_hi))

  pm <- perm_tbl %>% filter(panel == r$panel)
  cat(sprintf("  Permutation p: %.4f  (obs OOB = %.1f%%)\n", pm$p_perm, pm$obs_acc))
}
