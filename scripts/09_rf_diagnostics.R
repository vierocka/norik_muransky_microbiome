library(tidyverse)
library(randomForest)
library(pROC)
library(patchwork)

setwd("/home/veve/Dropbox/kone/qiime2")

contrib <- read.csv("pca_family_contributions.csv", check.names = FALSE)
fam_cols <- setdiff(colnames(contrib), c("sample_id", "topology"))
shorten_fam <- function(x) sub(
  "d__Bacteria;p__Patescibacteria;c__Saccharimonadia;o__Saccharimonadales;__",
  "Saccharimonadales", x)

dat <- contrib %>%
  filter(topology != "Environment") %>%
  mutate(topo3 = case_when(
    topology %in% c("Left front pastern", "Muzzle")   ~ "A_ground_contact",
    topology %in% c("Ventral abdomen",    "Udder")     ~ "B_near_ground",
    topology %in% c("Dorsum", "Forehead", "Neck",
                    "Pectoral area")                    ~ "C_elevated"
  ) %>% factor(levels = c("A_ground_contact", "B_near_ground", "C_elevated")))

X <- dat %>% select(all_of(fam_cols)) %>%
     mutate(across(everything(), as.numeric)) %>%
     rename_with(shorten_fam)

groups    <- levels(dat$topo3)
topo3_cols <- c("A_ground_contact" = "#2ca02c",
                "B_near_ground"    = "#17becf",
                "C_elevated"       = "#d62728")

# =========================================================
# 1. CONVERGENCE: OOB error vs number of trees
# Run once with 2000 trees and inspect the error trajectory
# =========================================================
cat("=== 1. RF convergence (2000 trees) ===\n")
set.seed(42)
rf2k <- randomForest(x = X, y = dat$topo3, ntree = 2000, importance = FALSE)

oob_df <- as.data.frame(rf2k$err.rate) %>%
  mutate(trees = row_number()) %>%
  pivot_longer(-trees, names_to = "class", values_to = "error")

# Check error at key checkpoints
for (n in c(100, 250, 500, 750, 1000, 1500, 2000)) {
  cat(sprintf("  ntrees=%4d  OOB=%.3f  A=%.3f  B=%.3f  C=%.3f\n", n,
              rf2k$err.rate[n, "OOB"],
              rf2k$err.rate[n, "A_ground_contact"],
              rf2k$err.rate[n, "B_near_ground"],
              rf2k$err.rate[n, "C_elevated"]))
}

p_conv <- ggplot(oob_df, aes(x = trees, y = error * 100, colour = class)) +
  geom_line(linewidth = 0.6, alpha = 0.85) +
  geom_vline(xintercept = 1000, linetype = "dashed", colour = "grey40",
             linewidth = 0.5) +
  annotate("text", x = 1020, y = 18, label = "n=1000", size = 3,
           hjust = 0, colour = "grey40") +
  scale_colour_manual(
    values = c("OOB" = "black", topo3_cols),
    labels = c("OOB" = "OOB (overall)", groups),
    name   = NULL) +
  theme_classic(base_size = 10) +
  theme(legend.position = c(0.75, 0.75),
        legend.text = element_text(size = 8)) +
  labs(title    = "RF convergence: OOB error vs number of trees",
       subtitle = "Dashed line = 1000 trees used in analysis",
       x = "Number of trees", y = "OOB error (%)")

# =========================================================
# 2. STABILITY: OOB error across 10 random seeds
# Tests whether results depend on random seed (variance across runs)
# =========================================================
cat("\n=== 2. RF stability across 10 random seeds ===\n")
seeds <- c(1, 7, 42, 99, 123, 256, 314, 500, 777, 999)
stab <- map_dfr(seeds, function(s) {
  set.seed(s)
  rf_s <- randomForest(x = X, y = dat$topo3, ntree = 1000, importance = FALSE)
  tibble(seed = s,
         OOB   = rf_s$err.rate[1000, "OOB"] * 100,
         A_err = rf_s$err.rate[1000, "A_ground_contact"] * 100,
         B_err = rf_s$err.rate[1000, "B_near_ground"]    * 100,
         C_err = rf_s$err.rate[1000, "C_elevated"]       * 100)
})
print(stab)
cat(sprintf("\nOOB mean=%.2f%%  SD=%.2f%%  range=[%.2f%%, %.2f%%]\n",
            mean(stab$OOB), sd(stab$OOB), min(stab$OOB), max(stab$OOB)))

p_stab <- stab %>%
  pivot_longer(-seed, names_to = "metric", values_to = "error") %>%
  mutate(metric = factor(metric, levels = c("OOB","A_err","B_err","C_err"))) %>%
  ggplot(aes(x = factor(seed), y = error, colour = metric, group = metric)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  scale_colour_manual(values = c("OOB" = "black", "A_err" = "#2ca02c",
                                  "B_err" = "#17becf", "C_err" = "#d62728"),
                      labels = c("OOB (overall)", groups), name = NULL) +
  theme_classic(base_size = 10) +
  theme(legend.position = c(0.75, 0.75), legend.text = element_text(size = 8)) +
  labs(title    = "RF stability: OOB error across 10 random seeds",
       subtitle = sprintf("OOB mean=%.2f%%, SD=%.2f%%", mean(stab$OOB), sd(stab$OOB)),
       x = "Random seed", y = "OOB error (%)")

# =========================================================
# 3. PROBABILITY CALIBRATION
# Checks if predicted probs are well-calibrated (no "overdispersion" of probs)
# Bins predicted probabilities and compares to observed fraction
# =========================================================
set.seed(42)
rf_main <- randomForest(x = X, y = dat$topo3, ntree = 1000, importance = FALSE)
probs   <- predict(rf_main, type = "prob")

cal_df <- map_dfr(groups, function(g) {
  obs  <- as.integer(dat$topo3 == g)
  pred <- probs[, g]
  bins <- cut(pred, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)
  tibble(group = g, prob_bin = bins, obs = obs, pred = pred) %>%
    group_by(group, prob_bin) %>%
    summarise(mean_pred = mean(pred), frac_obs = mean(obs), n = n(), .groups = "drop")
})

p_cal <- ggplot(cal_df, aes(x = mean_pred, y = frac_obs, colour = group)) +
  geom_abline(linetype = "dashed", colour = "grey60") +
  geom_point(aes(size = n), alpha = 0.8) +
  geom_line(alpha = 0.6) +
  scale_colour_manual(values = topo3_cols, name = NULL) +
  scale_size_continuous(range = c(1.5, 6), name = "n samples") +
  theme_classic(base_size = 10) +
  theme(legend.position = c(0.2, 0.75), legend.text = element_text(size = 8)) +
  labs(title    = "Probability calibration plot",
       subtitle = "Points near diagonal = well-calibrated; no overdispersion of probabilities",
       x = "Mean predicted probability", y = "Observed fraction")

# =========================================================
# 4. BOOTSTRAP CI STABILITY: compare 500 / 1000 / 2000 replicates
# =========================================================
cat("\n=== 4. Bootstrap CI stability (500 / 1000 / 2000 replicates) ===\n")
ci_stability <- map_dfr(c(500, 1000, 2000), function(nb) {
  map_dfr(groups, function(g) {
    set.seed(42)
    r  <- roc(ifelse(dat$topo3 == g, 1, 0), probs[, g], quiet = TRUE)
    ci <- ci.auc(r, method = "bootstrap", boot.n = nb,
                 conf.level = 0.95, progress = "none")
    tibble(boot_n = nb, group = g,
           auc    = round(as.numeric(r$auc), 4),
           ci_lo  = round(ci[1], 4),
           ci_hi  = round(ci[3], 4),
           width  = round(ci[3] - ci[1], 4))
  })
}, .progress = FALSE)
print(ci_stability)
cat("\n-> CI width stable across bootstrap replicates? ",
    ifelse(max(ci_stability$width) - min(ci_stability$width) < 0.01,
           "YES (width range < 0.01)", "CHECK"), "\n")

# =========================================================
# 5. YOUDEN'S INDEX + CI + FULL ROC METRICS
# Youden's J = sensitivity + specificity - 1 (maximised at optimal threshold)
# CI via bootstrap on coords at the optimal (Youden) threshold
# =========================================================
cat("\n=== 5. Youden's index and ROC metrics at optimal threshold ===\n")

youden_tbl <- map_dfr(groups, function(g) {
  resp <- ifelse(dat$topo3 == g, 1, 0)
  r    <- roc(resp, probs[, g], quiet = TRUE)

  # Optimal threshold (Youden)
  best <- coords(r, "best", ret = c("threshold","sensitivity","specificity","youden"),
                 best.method = "youden", transpose = FALSE)

  # PPV, NPV, F1 at optimal threshold
  pred_class <- ifelse(probs[, g] >= best$threshold, 1, 0)
  tp <- sum(pred_class == 1 & resp == 1)
  fp <- sum(pred_class == 1 & resp == 0)
  tn <- sum(pred_class == 0 & resp == 0)
  fn <- sum(pred_class == 0 & resp == 1)
  ppv <- ifelse((tp + fp) > 0, tp / (tp + fp), NA)
  npv <- ifelse((tn + fn) > 0, tn / (tn + fn), NA)
  f1  <- ifelse((2*tp + fp + fn) > 0, 2*tp / (2*tp + fp + fn), NA)

  # Youden's J = sensitivity + specificity - 1  (pROC returns sens+spec, so subtract 1)
  J <- best$sensitivity + best$specificity - 1

  # Bootstrap CI for Youden's J via manual resampling
  # (ci.coords fails at boundary when sens=1.0 for all bootstrap samples)
  set.seed(42)
  n    <- length(resp)
  J_boot <- replicate(2000, {
    idx  <- sample(n, n, replace = TRUE)
    r_b  <- tryCatch(
      roc(resp[idx], probs[idx, g], quiet = TRUE),
      error = function(e) NULL)
    if (is.null(r_b)) return(NA_real_)
    best_b <- coords(r_b, "best", ret = c("sensitivity","specificity"),
                     best.method = "youden", transpose = FALSE)
    best_b$sensitivity[1] + best_b$specificity[1] - 1
  })
  J_ci <- quantile(J_boot, c(0.025, 0.975), na.rm = TRUE)

  tibble(
    group       = g,
    threshold   = round(best$threshold, 3),
    sensitivity = round(best$sensitivity, 3),
    specificity = round(best$specificity, 3),
    Youden_J    = round(J, 3),
    youd_lo     = round(J_ci[1], 3),
    youd_hi     = round(J_ci[2], 3),
    PPV         = round(ppv, 3),
    NPV         = round(npv, 3),
    F1          = round(f1,  3)
  )
})

print(as.data.frame(youden_tbl))
write.csv(youden_tbl, "rf_roc_youden_metrics.csv", row.names = FALSE)
cat("Saved rf_roc_youden_metrics.csv\n")

# =========================================================
# 6. COMPOSE DIAGNOSTIC FIGURE
# =========================================================
fig <- (p_conv | p_stab) / (p_cal | plot_spacer()) +
  plot_annotation(
    title = "Random Forest diagnostics: convergence, stability, calibration",
    theme = theme(plot.title = element_text(size = 12, face = "bold"))
  )

ggsave("rf_diagnostics.png", fig, width = 13, height = 9, dpi = 300)
ggsave("rf_diagnostics.pdf", fig, width = 13, height = 9)
cat("Saved rf_diagnostics.png/.pdf\n")
