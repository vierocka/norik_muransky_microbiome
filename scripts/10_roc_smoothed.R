library(tidyverse)
library(randomForest)
library(pROC)

# NOTE: No independent test set available — ROC curves reflect training-set
# performance. OOB error from RF provides an approximately unbiased estimate
# of generalisation error (93.8% accuracy), but formal external validation
# requires a held-out cohort.

set.seed(42)
setwd("/home/veve/Dropbox/kone/qiime2")

contrib <- read.csv("pca_family_contributions.csv", check.names = FALSE)
imp_tbl <- read.csv("rf_family_importance.csv",     check.names = FALSE)

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

# -------------------------
# Refit RF and get predicted probabilities
# -------------------------
rf <- randomForest(x = X, y = dat$topo3, ntree = 1000,
                   importance = FALSE, keep.forest = TRUE)
probs  <- predict(rf, type = "prob")
groups <- levels(dat$topo3)

# -------------------------
# Compute raw ROC, smoothed ROC, and 95% CI for AUC
# Smoothing: binormal method (fits Gaussian to score distributions in each class)
# CI: bootstrap, 2000 replicates
# -------------------------
cat("Computing bootstrap CIs (2000 replicates) — may take ~30 s...\n")

roc_ls <- lapply(groups, function(g) {
  resp <- ifelse(dat$topo3 == g, 1, 0)
  list(
    raw    = roc(resp, probs[, g], quiet = TRUE),
    smooth = smooth(roc(resp, probs[, g], quiet = TRUE), method = "binormal"),
    ci     = ci.auc(roc(resp, probs[, g], quiet = TRUE),
                    method = "bootstrap", boot.n = 2000, conf.level = 0.95,
                    progress = "none")
  )
})
names(roc_ls) <- groups

# -------------------------
# Build data frame of smoothed curves
# -------------------------
roc_df <- map_dfr(groups, function(g) {
  s <- roc_ls[[g]]$smooth
  tibble(group       = g,
         specificity = s$specificities,
         sensitivity = s$sensitivities)
})

# Legend labels with AUC + 95% CI
legend_labels <- sapply(groups, function(g) {
  ci  <- roc_ls[[g]]$ci
  auc <- as.numeric(roc_ls[[g]]$raw$auc)
  sprintf("%s\nAUC = %.3f (95%% CI: %.3f-%.3f)", g, auc, ci[1], ci[3])
})

cat("\n=== AUC with 95% bootstrap CI ===\n")
for (g in groups) {
  ci <- roc_ls[[g]]$ci
  cat(sprintf("%-20s AUC = %.3f  95%% CI: %.3f - %.3f\n",
              g, as.numeric(roc_ls[[g]]$raw$auc), ci[1], ci[3]))
}

topo3_cols <- c("A_ground_contact" = "#2ca02c",
                "B_near_ground"    = "#17becf",
                "C_elevated"       = "#d62728")

p_roc <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity, colour = group)) +
  geom_line(linewidth = 1.1) +
  geom_abline(linetype = "dashed", colour = "grey60", linewidth = 0.6) +
  scale_colour_manual(values = topo3_cols, labels = legend_labels, name = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position  = c(0.64, 0.22),
        legend.text      = element_text(size = 8),
        legend.key.height = unit(1.1, "cm"),
        legend.background = element_rect(fill = "white", colour = "grey80",
                                         linewidth = 0.3),
        plot.title    = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8)) +
  labs(title    = "ROC curves: body-site contact zone (A/B/C) classification",
       subtitle = paste("Random Forest, 1000 trees | one-vs-rest | binormal smoothing (pROC)",
                        "| 95% CI: bootstrap (n=2000)", sep = "\n"),
       x = "1 - Specificity (False Positive Rate)",
       y = "Sensitivity (True Positive Rate)")

ggsave("roc_smoothed.png", p_roc, width = 6.5, height = 5.5, dpi = 300)
cat("Saved roc_smoothed.png\n")
