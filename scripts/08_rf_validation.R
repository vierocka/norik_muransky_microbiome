library(tidyverse)
library(randomForest)
library(pROC)
library(nnet)
library(patchwork)

set.seed(42)
setwd("/home/veve/Dropbox/kone/qiime2")

# -------------------------
# Load data
# -------------------------
contrib  <- read.csv("pca_family_contributions.csv", check.names = FALSE)
meta     <- read.csv("/home/veve/Dropbox/kone/kraken2_reports/metadata.csv",
                     check.names = FALSE) %>%
            mutate(across(everything(), trimws))
imp_tbl  <- read.csv("rf_family_importance.csv",     check.names = FALSE)

fam_cols <- setdiff(colnames(contrib), c("sample_id", "topology"))

shorten_fam <- function(x) sub(
  "d__Bacteria;p__Patescibacteria;c__Saccharimonadia;o__Saccharimonadales;__",
  "Saccharimonadales", x)

top8  <- imp_tbl$family[1:8]
top4  <- imp_tbl$family[1:4]

# Animal samples with topo3 grouping
dat <- contrib %>%
  filter(topology != "Environment") %>%
  mutate(topo3 = case_when(
    topology %in% c("Left front pastern", "Muzzle")   ~ "GCtS_A",
    topology %in% c("Ventral abdomen",    "Udder")     ~ "GCtS_B",
    topology %in% c("Dorsum", "Forehead", "Neck",
                    "Pectoral area")                    ~ "EtS"
  ) %>% factor(levels = c("GCtS_A", "GCtS_B", "EtS")))

X_all <- dat %>% select(all_of(fam_cols)) %>%
         mutate(across(everything(), as.numeric)) %>%
         rename_with(shorten_fam)

topo3_cols <- c("GCtS_A" = "#2ca02c",
                "GCtS_B"    = "#17becf",
                "EtS"       = "#d62728")
topology_cols <- c(
  "Dorsum"             = "#1f77b4", "Forehead"           = "#ff7f0e",
  "Left front pastern" = "#2ca02c", "Muzzle"             = "#9467bd",
  "Neck"               = "#8c564b", "Pectoral area"      = "#7f7f7f",
  "Udder"              = "#bcbd22", "Ventral abdomen"    = "#17becf"
)

# =========================================================
# 1. OVERDISPERSION: raw counts for top 8 families
# =========================================================
fam_raw   <- read.delim("family_table.tsv", check.names = FALSE,
                         comment.char = "", row.names = 1)
colnames(fam_raw) <- sub("_S\\d+$", "", colnames(fam_raw))
fam_mat   <- as.matrix(fam_raw); mode(fam_mat) <- "numeric"
fam_names <- sub("^.*f__", "", rownames(fam_mat))
fam_agg   <- rowsum(fam_mat, fam_names)
# Animal samples only
anim_ids  <- meta$sample_id[meta$topology != "Environment"]

cat("=== Overdispersion check: raw counts for top 8 families ===\n")
cat(sprintf("%-35s  %8s  %8s  %8s  %8s  %s\n",
            "Family", "Mean", "Variance", "Var/Mean", "Zeros%", "Verdict"))
disp_tbl <- map_dfr(top8, function(f) {
  # match family name (top8 may be shortened)
  raw_name <- names(fam_names)[fam_names == f]
  if (length(raw_name) == 0) raw_name <- names(fam_names)[grepl(f, fam_names)]
  counts <- if (length(raw_name) > 0) as.numeric(fam_agg[raw_name[1], anim_ids]) else rep(0, length(anim_ids))
  mn  <- mean(counts);  vr  <- var(counts)
  disp <- round(vr / mn, 1)
  cat(sprintf("%-35s  %8.1f  %8.1f  %8.1f  %7.0f%%  %s\n",
              f, mn, vr, disp,
              100 * mean(counts == 0),
              ifelse(disp > 10, "overdispersed -> NB", "moderate -> Poisson ok")))
  tibble(family = f, mean_count = round(mn,1), variance = round(vr,1),
         dispersion = disp, pct_zeros = round(100*mean(counts==0),0))
})
write.csv(disp_tbl, "top8_overdispersion.csv", row.names = FALSE)
cat("Saved top8_overdispersion.csv\n")

# =========================================================
# 2. MULTINOMIAL LOGISTIC REGRESSION: top 3 and top 4
# =========================================================
fit_multinom <- function(preds, label) {
  df <- X_all %>% select(all_of(preds)) %>%
        mutate(group = dat$topo3)
  f  <- as.formula(paste("group ~", paste(paste0("`", preds, "`"), collapse = " + ")))
  m  <- multinom(f, data = df, trace = FALSE)
  acc <- mean(predict(m) == df$group)
  cat(sprintf("\n=== Multinomial logistic (%s) | accuracy = %.1f%% ===\n",
              label, acc * 100))
  cat("Coefficients:\n"); print(round(coef(m), 3))
  cat("Confusion matrix:\n")
  print(table(predicted = predict(m), actual = df$group))
  list(model = m, data = df, acc = acc, preds = preds)
}

res3 <- fit_multinom(top4[1:3], "top 3")
res4 <- fit_multinom(top4,      "top 4")

# =========================================================
# 3. ROC CURVES (one-vs-rest, RF predicted probabilities)
# =========================================================
rf_full <- randomForest(x = X_all, y = dat$topo3, ntree = 1000,
                        importance = FALSE, keep.forest = TRUE)
probs   <- predict(rf_full, type = "prob")

groups  <- levels(dat$topo3)
roc_ls  <- lapply(groups, function(g) {
  roc(response  = ifelse(dat$topo3 == g, 1, 0),
      predictor = probs[, g], quiet = TRUE)
})
names(roc_ls) <- groups

cat("\n=== AUC (one-vs-rest, RF) ===\n")
auc_tbl <- tibble(
  group = groups,
  AUC   = sapply(roc_ls, \(r) round(as.numeric(auc(r)), 3))
)
print(auc_tbl)
write.csv(auc_tbl, "rf_roc_auc.csv", row.names = FALSE)

roc_df <- map_dfr(groups, function(g) {
  r <- roc_ls[[g]]
  tibble(group = g, specificity = r$specificities, sensitivity = r$sensitivities)
})

roc_labels <- auc_tbl %>%
  mutate(label = sprintf("%s\nAUC=%.3f", group, AUC))

p_roc <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity, colour = group)) +
  geom_line(linewidth = 1) +
  geom_abline(linetype = "dashed", colour = "grey60") +
  scale_colour_manual(values = topo3_cols,
                      labels = setNames(
                        sprintf("%s  AUC=%.3f", groups, auc_tbl$AUC), groups),
                      name = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = c(0.65, 0.2),
        legend.text     = element_text(size = 8)) +
  labs(title    = "ROC curves: contact zone classification (one-vs-rest)",
       subtitle = "Random Forest, 1000 trees | OOB accuracy = 93.8%",
       x = "1 - Specificity", y = "Sensitivity")

ggsave("rf_roc_curves.pdf", p_roc, width = 6, height = 5)
ggsave("rf_roc_curves.png", p_roc, width = 6, height = 5, dpi = 300)
cat("Saved rf_roc_curves.pdf/.png\n")

# =========================================================
# 4. PCA ON TOP FAMILIES
# =========================================================
pca_plot <- function(X_sub, colour_var, colour_vals, colour_label,
                     title_str, stem, shape_var = NULL) {
  pca  <- prcomp(X_sub, center = TRUE, scale. = FALSE)
  imp  <- summary(pca)$importance
  df   <- as.data.frame(pca$x[, 1:2]) %>%
          mutate(colour = colour_var,
                 shape  = if (!is.null(shape_var)) shape_var else "all")
  p <- ggplot(df, aes(x = PC1, y = PC2, colour = colour)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_colour_manual(values = colour_vals, name = colour_label) +
    theme_classic(base_size = 11) +
    labs(title    = title_str,
         subtitle = sprintf("%d families | PC1=%.1f%% | PC2=%.1f%%",
                            ncol(X_sub),
                            imp[2, 1] * 100, imp[2, 2] * 100),
         x = sprintf("PC1 (%.1f%%)", imp[2,1]*100),
         y = sprintf("PC2 (%.1f%%)", imp[2,2]*100))
  ggsave(paste0(stem, ".pdf"), p, width = 7, height = 5)
  ggsave(paste0(stem, ".png"), p, width = 7, height = 5, dpi = 300)
  cat("Saved", stem, ".pdf/.png\n")
  p
}

X4 <- X_all %>% select(all_of(top4))
X8 <- X_all %>% select(all_of(top8))

# Top 4 — by A/B/C group
p4_abc  <- pca_plot(X4, as.character(dat$topo3), topo3_cols,
                    "Contact zone",
                    "PCA on top 4 RF predictors — contact zone",
                    "pca_top4_topo3")

# Top 8 — by A/B/C group
p8_abc  <- pca_plot(X8, as.character(dat$topo3), topo3_cols,
                    "Contact zone",
                    "PCA on top 8 RF predictors — contact zone",
                    "pca_top8_topo3")

# Top 8 — by individual body site
p8_site <- pca_plot(X8, dat$topology, topology_cols,
                    "Body site",
                    "PCA on top 8 RF predictors — body site",
                    "pca_top8_site")

# =========================================================
# 5. KRAKEN2 OVERLAP
# =========================================================
kr_fams <- colnames(read.csv(
  "/home/veve/Dropbox/kone/kraken2_reports/pca_family_contributions.csv",
  check.names = FALSE, nrows = 1))
kr_fams <- setdiff(kr_fams, c("sample_id", "topology"))

cat("\n=== Kraken2 overlap: top 8 RF predictors ===\n")
overlap_tbl <- tibble(
  rank        = 1:8,
  family      = top8,
  source      = imp_tbl$source[1:8],
  in_kraken2  = top8 %in% kr_fams
)
print(overlap_tbl)
write.csv(overlap_tbl, "top8_kraken2_overlap.csv", row.names = FALSE)
cat("Saved top8_kraken2_overlap.csv\n")
