#!/usr/bin/env Rscript
# 16_probe_population_validation.R — Population-level probe validation
#
# Tests each of the 8 candidate probe families as a standalone predictor and
# validates both 4-family panels at the population level.
#
# Analyses:
#   1. Per-probe ROC + AUC + Youden J / Sens / Spec  (bootstrap CI, B=2000)
#   2. Per-probe within-horse permutation p-value      (N=999)
#   3. Per-probe noise stability                       (1, 2, 5, 10 %)
#   4. Leave-one-horse-out (LOHO) cross-validation
#      — per probe (logistic regression)
#      — per panel (Random Forest)
#      Simulates how probes/panels generalise to an unseen horse population.

library(tidyverse)
library(pROC)
library(randomForest)

set.seed(42)
N_BOOT   <- 2000
N_PERM   <- 999
N_SEEDS  <- 20
NOISE_LEVELS <- c(0.01, 0.02, 0.05, 0.10)

PANEL_A <- c("Deinococcaceae", "Intrasporangiaceae", "Hungateiclostridiaceae", "Christensenellaceae")
PANEL_B <- c("Deinococcaceae", "Micrococcaceae",     "Lachnospiraceae",        "UCG-010")
ALL_PROBES <- union(PANEL_A, PANEL_B)   # 7 unique families

OUT_TABLES  <- "results/tables"
OUT_FIGURES <- "results/figures"

# =============================================================================
# Data
# =============================================================================
contrib <- read_csv("data/pca_family_contributions.csv", show_col_types = FALSE)
meta    <- read_csv("data/metadata.csv",                 show_col_types = FALSE)

dat <- contrib %>%
  left_join(select(meta, sample_id, age), by = "sample_id") %>%
  filter(topology != "Environment") %>%
  mutate(
    horse_id = sub("-.*", "", sample_id) %>% factor(),
    topo3 = case_when(
      topology %in% c("Left front pastern", "Muzzle")                  ~ "GCtS_A",
      topology %in% c("Ventral abdomen",    "Udder")                    ~ "GCtS_B",
      topology %in% c("Dorsum", "Forehead", "Neck", "Pectoral area")   ~ "EtS"
    ) %>% factor(levels = c("GCtS_A", "GCtS_B", "EtS")),
    group2 = ifelse(topo3 == "EtS", "EtS", "GCtS") %>%
             factor(levels = c("GCtS", "EtS"))
  )

cat("=== Dataset ===\n")
cat("Horses:", n_distinct(dat$horse_id),
    " | Samples:", nrow(dat),
    " | GCtS:", sum(dat$group2 == "GCtS"),
    " | EtS:", sum(dat$group2 == "EtS"), "\n\n")

horses <- levels(dat$horse_id)

# =============================================================================
# Helpers
# =============================================================================

roc_metrics_boot <- function(labels, score, n_boot = N_BOOT) {
  # direction = "auto" so pROC picks the better orientation for each family.
  # GCtS-enriched families (higher abundance → GCtS, not EtS) would otherwise
  # show AUC < 0.5 with direction = "<". Auto-detection gives the undirected
  # discriminatory power, equivalent to max(AUC, 1-AUC).
  roc_obj  <- roc(labels, score, levels = c("GCtS", "EtS"),
                  direction = "auto", quiet = TRUE)
  auc_obs  <- as.numeric(auc(roc_obj))
  yi       <- which.max(roc_obj$sensitivities + roc_obj$specificities - 1)
  thr_obs  <- roc_obj$thresholds[yi]
  sens_obs <- roc_obj$sensitivities[yi]
  spec_obs <- roc_obj$specificities[yi]
  j_obs    <- sens_obs + spec_obs - 1

  n <- length(labels)
  boot_mat <- replicate(n_boot, {
    idx <- sample(n, n, replace = TRUE)
    lb  <- labels[idx]; sc <- score[idx]
    if (length(unique(lb)) < 2) return(rep(NA_real_, 4))
    r <- tryCatch(roc(lb, sc, levels = c("GCtS", "EtS"),
                      direction = "auto", quiet = TRUE),
                  error = function(e) NULL)
    if (is.null(r)) return(rep(NA_real_, 4))
    yi2 <- which.max(r$sensitivities + r$specificities - 1)
    c(auc  = as.numeric(auc(r)),
      sens = r$sensitivities[yi2],
      spec = r$specificities[yi2],
      j    = r$sensitivities[yi2] + r$specificities[yi2] - 1)
  })

  ci95 <- function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE)
  list(roc_obj = roc_obj,
       auc  = auc_obs,  auc_ci  = ci95(boot_mat["auc",]),
       sens = sens_obs, sens_ci = ci95(boot_mat["sens",]),
       spec = spec_obs, spec_ci = ci95(boot_mat["spec",]),
       j    = j_obs,    j_ci    = ci95(boot_mat["j",]),
       thr  = thr_obs)
}

perm_p <- function(labels, score, horse_id, n_perm = N_PERM) {
  roc_obs <- roc(labels, score, levels = c("GCtS","EtS"),
                 direction = "auto", quiet = TRUE)
  obs_auc <- as.numeric(auc(roc_obs))
  null_auc <- replicate(n_perm, {
    perm_score <- score
    for (h in levels(horse_id)) {
      idx <- which(horse_id == h)
      perm_score[idx] <- sample(score[idx])
    }
    r <- tryCatch(roc(labels, perm_score, levels=c("GCtS","EtS"),
                      direction="auto", quiet=TRUE), error=function(e) NULL)
    if (is.null(r)) return(NA_real_)
    as.numeric(auc(r))
  })
  null_auc <- null_auc[!is.na(null_auc)]
  mean(null_auc >= obs_auc)
}

noise_probe_stability <- function(ab, labels, noise_levels = NOISE_LEVELS, n_seeds = N_SEEDS) {
  all_levels <- c(0, noise_levels)
  map_dfr(all_levels, function(pct) {
    auc_vec <- map_dbl(seq_len(n_seeds), function(s) {
      set.seed(s)
      sc <- if (pct == 0) ab else pmax(ab + ab * rnorm(length(ab), 0, pct), 0)
      r  <- tryCatch(roc(labels, sc, levels=c("GCtS","EtS"),
                         direction="auto", quiet=TRUE), error=function(e) NULL)
      if (is.null(r)) return(NA_real_)
      as.numeric(auc(r))
    })
    tibble(noise_pct = pct * 100,
           mean_auc  = mean(auc_vec, na.rm=TRUE),
           sd_auc    = sd(auc_vec,   na.rm=TRUE))
  })
}

add_noise <- function(X, pct) {
  noise <- matrix(rnorm(prod(dim(X)), 0, pct), nrow = nrow(X))
  pmax(X + X * noise, 0)
}

# =============================================================================
# 1. Per-probe ROC metrics + permutation + noise
# =============================================================================
cat("=== 1. Per-probe validation ===\n\n")

probe_roc_rows <- list()
probe_roc_curves <- list()
probe_perm_rows  <- list()
probe_noise_rows <- list()

for (fam in ALL_PROBES) {
  cat(sprintf("  %s ...\n", fam))
  ab     <- dat[[fam]]
  labels <- dat$group2
  horse  <- dat$horse_id

  # ROC + Youden
  m <- roc_metrics_boot(labels, ab)

  probe_roc_rows[[fam]] <- tibble(
    family  = fam,
    panel   = case_when(fam %in% PANEL_A & fam %in% PANEL_B ~ "Both",
                        fam %in% PANEL_A ~ "Panel A only",
                        TRUE             ~ "Panel B only"),
    auc     = m$auc,  auc_lo  = m$auc_ci[1],  auc_hi  = m$auc_ci[2],
    youden_j= m$j,    j_lo    = m$j_ci[1],    j_hi    = m$j_ci[2],
    sens    = m$sens, sens_lo = m$sens_ci[1],  sens_hi = m$sens_ci[2],
    spec    = m$spec, spec_lo = m$spec_ci[1],  spec_hi = m$spec_ci[2],
    opt_thr = m$thr
  )

  probe_roc_curves[[fam]] <- tibble(
    family    = fam,
    fpr       = 1 - m$roc_obj$specificities,
    tpr       = m$roc_obj$sensitivities,
    threshold = m$roc_obj$thresholds
  )

  # Permutation
  p_val <- perm_p(labels, ab, horse)
  probe_perm_rows[[fam]] <- tibble(family = fam, p_perm = p_val)

  # Noise
  nd <- noise_probe_stability(ab, labels)
  nd$family <- fam
  probe_noise_rows[[fam]] <- nd
}

probe_roc_tbl   <- bind_rows(probe_roc_rows) %>%
  left_join(bind_rows(probe_perm_rows), by = "family") %>%
  arrange(desc(auc))
probe_noise_tbl <- bind_rows(probe_noise_rows)
probe_roc_curve_tbl <- bind_rows(probe_roc_curves)

write.csv(probe_roc_tbl,       file.path(OUT_TABLES, "probe_roc_metrics.csv"),       row.names=FALSE)
write.csv(probe_noise_tbl,     file.path(OUT_TABLES, "probe_noise_stability.csv"),   row.names=FALSE)
write.csv(probe_roc_curve_tbl, file.path(OUT_TABLES, "probe_roc_curves.csv"),        row.names=FALSE)
cat("Saved probe_roc_metrics.csv | probe_noise_stability.csv | probe_roc_curves.csv\n\n")

cat("Per-probe summary (sorted by AUC):\n")
print(probe_roc_tbl %>%
  select(family, panel, auc, auc_lo, auc_hi, youden_j, sens, spec, p_perm) %>%
  mutate(across(where(is.numeric), ~round(., 3))), n = 10)

# =============================================================================
# 2. Leave-one-horse-out (LOHO) cross-validation
# =============================================================================
cat("\n=== 2. Leave-one-horse-out cross-validation ===\n\n")

# --- 2a. Per-probe: logistic regression, predict held-out horse ---
cat("  2a. Per-probe logistic regression LOHO ...\n")

loho_probe <- map_dfr(ALL_PROBES, function(fam) {
  preds <- map_dfr(horses, function(h) {
    train <- dat %>% filter(horse_id != h)
    test  <- dat %>% filter(horse_id == h)
    form  <- as.formula(paste0("group2 ~ `", fam, "`"))
    fit   <- tryCatch(
      glm(form, data = train, family = binomial),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)
    prob_ets <- predict(fit, newdata = test, type = "response")
    tibble(
      horse  = h,
      true   = test$group2,
      prob   = prob_ets,
      pred   = ifelse(prob_ets > 0.5, "EtS", "GCtS") %>%
               factor(levels = c("GCtS", "EtS"))
    )
  })
  acc <- mean(preds$pred == preds$true)
  roc_obj <- tryCatch(
    roc(preds$true, preds$prob, levels=c("GCtS","EtS"),
        direction="<", quiet=TRUE),
    error = function(e) NULL
  )
  auc_val <- if (is.null(roc_obj)) NA_real_ else as.numeric(auc(roc_obj))
  tibble(family   = fam,
         loho_acc = acc * 100,
         loho_auc = auc_val)
})

cat("  Per-probe LOHO results:\n")
print(loho_probe %>% arrange(desc(loho_auc)) %>%
  mutate(across(where(is.numeric), ~round(., 3))), n = 10)

# --- 2b. Per-panel: Random Forest LOHO ---
cat("\n  2b. Panel RF LOHO ...\n")

panels <- list(Panel_A = PANEL_A, Panel_B = PANEL_B)

loho_panel <- map_dfr(names(panels), function(pname) {
  fams <- panels[[pname]]
  preds <- map_dfr(horses, function(h) {
    train <- dat %>% filter(horse_id != h)
    test  <- dat %>% filter(horse_id == h)
    X_tr  <- train %>% select(all_of(fams)) %>% as.matrix()
    X_te  <- test  %>% select(all_of(fams)) %>% as.matrix()
    y_tr  <- train$group2
    set.seed(42)
    rf  <- randomForest(x = X_tr, y = y_tr, ntree = 500)
    votes <- predict(rf, newdata = X_te, type = "vote")
    tibble(
      horse = h,
      true  = test$group2,
      prob  = votes[, "EtS"],
      pred  = predict(rf, newdata = X_te)
    )
  })
  acc <- mean(preds$pred == preds$true)
  roc_obj <- tryCatch(
    roc(preds$true, preds$prob, levels=c("GCtS","EtS"),
        direction="<", quiet=TRUE),
    error = function(e) NULL
  )
  auc_val <- if (is.null(roc_obj)) NA_real_ else as.numeric(auc(roc_obj))
  tibble(panel    = pname,
         families = paste(fams, collapse="; "),
         loho_acc = acc * 100,
         loho_auc = auc_val)
})

cat("  Panel LOHO results:\n")
print(loho_panel %>% mutate(across(where(is.numeric), ~round(., 3))))

write.csv(loho_probe, file.path(OUT_TABLES, "probe_loho_cv.csv"),  row.names=FALSE)
write.csv(loho_panel, file.path(OUT_TABLES, "panel_loho_cv.csv"), row.names=FALSE)
cat("Saved probe_loho_cv.csv | panel_loho_cv.csv\n")

# =============================================================================
# Figures
# =============================================================================

fam_cols <- setNames(
  c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
    "#9467bd", "#8c564b", "#e377c2"),
  ALL_PROBES
)
panel_membership <- c(
  "Deinococcaceae"               = "Both panels",
  "Intrasporangiaceae"           = "Panel A only",
  "Hungateiclostridiaceae"       = "Panel A only",
  "Christensenellaceae"          = "Panel A only",
  "Micrococcaceae"               = "Panel B only",
  "Lachnospiraceae"              = "Panel B only",
  "UCG-010"                      = "Panel B only"
)

# Figure 1 — Per-probe ROC curves
p_probe_roc <- ggplot(probe_roc_curve_tbl,
                      aes(x = fpr, y = tpr, colour = family)) +
  geom_line(linewidth = 0.9) +
  geom_abline(linetype = "dashed", colour = "grey60") +
  scale_colour_manual(values = fam_cols,
                      name   = "Family") +
  coord_equal() +
  theme_classic(base_size = 11) +
  theme(legend.position = "right",
        legend.text     = element_text(size = 8),
        legend.key.size = unit(0.5, "cm")) +
  labs(title    = "ROC curves — individual probe families (GCtS vs EtS)",
       subtitle  = sprintf("Score = raw family abundance, B=%d bootstrap CIs", N_BOOT),
       x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)")

ggsave(file.path(OUT_FIGURES, "probe_roc_curves.pdf"), p_probe_roc, width = 7, height = 5)
ggsave(file.path(OUT_FIGURES, "probe_roc_curves.png"), p_probe_roc, width = 7, height = 5, dpi = 300)
cat("Saved probe_roc_curves.pdf/.png\n")

# Figure 2 — Per-probe forest plot: AUC + Youden J
forest_df <- probe_roc_tbl %>%
  mutate(family_label = paste0(family,
                               ifelse(p_perm < 0.05, "*", ""),
                               ifelse(p_perm < 0.001, "**", "")),
         family_label = reorder(family_label, auc),
         membership   = panel_membership[family])

p_forest <- ggplot(forest_df, aes(y = family_label, colour = membership)) +
  # AUC
  geom_linerange(aes(xmin = auc_lo, xmax = auc_hi, x = auc),
                 position = position_nudge(y = 0.15), linewidth = 0.8) +
  geom_point(aes(x = auc, shape = "AUC"),
             position = position_nudge(y = 0.15), size = 3) +
  # Youden J
  geom_linerange(aes(xmin = j_lo, xmax = j_hi, x = youden_j),
                 position = position_nudge(y = -0.15), linewidth = 0.8) +
  geom_point(aes(x = youden_j, shape = "Youden J"),
             position = position_nudge(y = -0.15), size = 3) +
  scale_colour_manual(values = c("Both panels"  = "#2ca02c",
                                 "Panel A only" = "#1f77b4",
                                 "Panel B only" = "#d62728"),
                      name = NULL) +
  scale_shape_manual(values = c("AUC" = 16, "Youden J" = 17), name = NULL) +
  scale_x_continuous(limits = c(0.3, 1.02), breaks = seq(0.3, 1, 0.1)) +
  geom_vline(xintercept = 0.5, linetype = "dotted", colour = "grey60") +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom",
        axis.title.y    = element_blank()) +
  labs(title    = "Per-probe discriminatory power (GCtS vs EtS)",
       subtitle  = "95% bootstrap CI | * p_perm<0.05  ** p_perm<0.001",
       x = "Metric value")

ggsave(file.path(OUT_FIGURES, "probe_forest_plot.pdf"), p_forest, width = 7, height = 5)
ggsave(file.path(OUT_FIGURES, "probe_forest_plot.png"), p_forest, width = 7, height = 5, dpi = 300)
cat("Saved probe_forest_plot.pdf/.png\n")

# Figure 3 — Noise stability per probe (line plot)
noise_plot_df <- probe_noise_tbl %>%
  mutate(noise_label = paste0(noise_pct, "%"),
         noise_label = factor(noise_label, levels = paste0(c(0,1,2,5,10),"%")))

p_noise_probe <- ggplot(noise_plot_df,
                        aes(x = noise_pct, y = mean_auc,
                            colour = family, group = family)) +
  geom_line(linewidth = 0.9) +
  geom_ribbon(aes(ymin = mean_auc - sd_auc, ymax = mean_auc + sd_auc,
                  fill = family), alpha = 0.15, colour = NA) +
  scale_colour_manual(values = fam_cols, name = "Family") +
  scale_fill_manual(  values = fam_cols, name = "Family") +
  scale_x_continuous(breaks = c(0, 1, 2, 5, 10),
                     labels = paste0(c(0, 1, 2, 5, 10), "%")) +
  scale_y_continuous(limits = c(0.0, 1.0), breaks = seq(0, 1, 0.2)) +
  theme_classic(base_size = 11) +
  theme(legend.position = "right",
        legend.text     = element_text(size = 8)) +
  labs(title    = "Per-probe AUC stability under multiplicative noise",
       subtitle  = sprintf("Mean +/- SD over %d seeds", N_SEEDS),
       x = "Added noise level", y = "AUC")

ggsave(file.path(OUT_FIGURES, "probe_noise_stability.pdf"), p_noise_probe, width = 7, height = 5)
ggsave(file.path(OUT_FIGURES, "probe_noise_stability.png"), p_noise_probe, width = 7, height = 5, dpi = 300)
cat("Saved probe_noise_stability.pdf/.png\n")

# Figure 4 — LOHO accuracy: probes vs panels
loho_probe_plot <- loho_probe %>%
  mutate(type = "Single probe",
         label = family) %>%
  select(label, type, loho_acc, loho_auc)

loho_panel_plot <- loho_panel %>%
  mutate(type  = "4-family panel",
         label = panel) %>%
  select(label, type, loho_acc, loho_auc)

loho_all <- bind_rows(loho_probe_plot, loho_panel_plot) %>%
  mutate(label = reorder(label, loho_auc))

p_loho <- ggplot(loho_all, aes(x = loho_auc, y = label, colour = type, shape = type)) +
  geom_point(size = 4) +
  geom_vline(xintercept = 0.5, linetype = "dotted", colour = "grey60") +
  scale_colour_manual(values = c("Single probe" = "#7f7f7f",
                                 "4-family panel" = "#d62728"),
                      name = NULL) +
  scale_shape_manual(values  = c("Single probe" = 16,
                                 "4-family panel" = 18),
                     name = NULL) +
  scale_x_continuous(limits = c(0.3, 1.02), breaks = seq(0.3, 1, 0.1)) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom",
        axis.title.y    = element_blank()) +
  labs(title    = "Leave-one-horse-out AUC: single probes vs. 4-family panels",
       subtitle  = "Logistic regression (probes) | Random Forest (panels) | N=12 horses",
       x = "LOHO AUC")

ggsave(file.path(OUT_FIGURES, "probe_loho_comparison.pdf"), p_loho, width = 7, height = 5)
ggsave(file.path(OUT_FIGURES, "probe_loho_comparison.png"), p_loho, width = 7, height = 5, dpi = 300)
cat("Saved probe_loho_comparison.pdf/.png\n")

# =============================================================================
# Summary
# =============================================================================
cat("\n=== SUMMARY ===\n\n")
cat("Per-probe AUC and permutation p (sorted by AUC):\n")
print(probe_roc_tbl %>%
  select(family, panel, auc, auc_lo, auc_hi, youden_j, p_perm) %>%
  mutate(across(where(is.numeric), ~round(., 3))), n = 10)

cat("\nLOHO cross-validation:\n")
cat("  Probes:\n")
print(loho_probe %>% arrange(desc(loho_auc)) %>%
  mutate(across(where(is.numeric), ~round(., 3))), n=10)
cat("  Panels:\n")
print(loho_panel %>% select(panel, loho_acc, loho_auc) %>%
  mutate(across(where(is.numeric), ~round(., 3))))
