library(multcomp)
library(tidyverse)
library(randomForest)
library(pROC)
library(nnet)
library(patchwork)

# A+B merged into GCtS (ground-contact topology sites)
# vs EtS (elevated topology sites)
# Rationale: A vs B non-significant in Tukey (p=0.588); both overlap in PCA;
# biologically both share soil/grass environmental exposure.

set.seed(42)
setwd("/home/veve/Dropbox/kone/qiime2")

contrib <- read.csv("pca_family_contributions.csv", check.names = FALSE)
meta    <- read.csv("/home/veve/Dropbox/kone/kraken2_reports/metadata.csv",
                    check.names = FALSE) %>%
           mutate(across(everything(), trimws))

fam_cols <- setdiff(colnames(contrib), c("sample_id", "topology"))
shorten_fam <- function(x) sub(
  "d__Bacteria;p__Patescibacteria;c__Saccharimonadia;o__Saccharimonadales;__",
  "Saccharimonadales", x)

# -------------------------
# Build dataset: 2 groups
# GCtS = Muzzle, Pastern, Udder, Ventral abdomen
# EtS  = Dorsum, Forehead, Neck, Pectoral area
# -------------------------
dat <- contrib %>%
  filter(topology != "Environment") %>%
  mutate(
    richness = rowSums(across(all_of(fam_cols), as.numeric) > 0),
    group2   = case_when(
      topology %in% c("Left front pastern", "Muzzle",
                      "Ventral abdomen",    "Udder")   ~ "GCtS",
      topology %in% c("Dorsum", "Forehead", "Neck",
                      "Pectoral area")                  ~ "EtS"
    ) %>% factor(levels = c("GCtS", "EtS"))
  ) %>%
  left_join(meta %>% select(sample_id, age), by = "sample_id")

group2_cols <- c("GCtS" = "#2ca02c", "EtS" = "#d62728")
topology_cols <- c(
  "Dorsum"             = "#1f77b4", "Forehead"           = "#ff7f0e",
  "Left front pastern" = "#2ca02c", "Muzzle"             = "#9467bd",
  "Neck"               = "#8c564b", "Pectoral area"      = "#7f7f7f",
  "Udder"              = "#bcbd22", "Ventral abdomen"    = "#17becf"
)

cat("=== Group sizes ===\n")
print(dat %>% count(group2, topology) %>% arrange(group2, topology))

X <- dat %>% select(all_of(fam_cols)) %>%
     mutate(across(everything(), as.numeric)) %>%
     rename_with(shorten_fam)

# =========================================================
# 1. RICHNESS: ANOVA + Wilcoxon GCtS vs EtS
# =========================================================
cat("\n=== 1. Richness: GCtS vs EtS ===\n")
wt_rich <- wilcox.test(richness ~ group2, data = dat, exact = FALSE)
m_rich  <- lm(richness ~ group2, data = dat)
cat(sprintf("Wilcoxon W=%.0f, p=%.4f\n", wt_rich$statistic, wt_rich$p.value))
cat(sprintf("Linear model F=%.3f, p=%.4f\n",
            summary(m_rich)$fstatistic[1],
            pf(summary(m_rich)$fstatistic[1],
               summary(m_rich)$fstatistic[2],
               summary(m_rich)$fstatistic[3], lower.tail=FALSE)))
rich_sum <- dat %>% group_by(group2) %>%
  summarise(n=n(), mean=round(mean(richness),1), sd=round(sd(richness),1),
            median=median(richness), .groups="drop")
print(rich_sum)
write.csv(rich_sum, "richness_ACgroups_onlyACgroups.csv", row.names=FALSE)

p_rich <- ggplot(dat, aes(x=group2, y=richness, fill=group2)) +
  geom_boxplot(alpha=0.8, outlier.shape=NA, width=0.5) +
  geom_jitter(aes(colour=group2), width=0.12, size=2, alpha=0.65) +
  scale_fill_manual(values=group2_cols)  +
  scale_colour_manual(values=group2_cols) +
  theme_classic(base_size=12) +
  theme(legend.position="none",
        plot.title=element_text(size=10,face="bold"),
        plot.subtitle=element_text(size=8)) +
  labs(title="Family richness: ground-contact (GCtS) vs elevated (EtS)",
       subtitle=sprintf("Wilcoxon W=%.0f, p=%.4f", wt_rich$statistic, wt_rich$p.value),
       x=NULL, y="Families present")
ggsave("richness_boxplot_onlyACgroups.pdf", p_rich, width=5, height=5)
ggsave("richness_boxplot_onlyACgroups.png", p_rich, width=5, height=5, dpi=300)
cat("Saved richness_boxplot_onlyACgroups.pdf/.png\n")

# =========================================================
# 2. FAMILY ABUNDANCE: Wilcoxon per family, GCtS vs EtS (BH)
# =========================================================
cat("\n=== 2. Family abundance: Wilcoxon GCtS vs EtS (BH) ===\n")
famab <- map_dfr(colnames(X), function(f) {
  ac  <- X[[f]][dat$group2 == "GCtS"]
  ce  <- X[[f]][dat$group2 == "EtS"]
  wt  <- wilcox.test(ac, ce, exact=FALSE)
  tibble(family    = f,
         med_AC    = round(median(ac)*100, 4),
         med_C     = round(median(ce)*100, 4),
         log2FC    = round(log2((mean(ac)+1e-6)/(mean(ce)+1e-6)), 3),
         W         = round(wt$statistic, 1),
         p         = wt$p.value)
}) %>%
  mutate(p_adj = p.adjust(p, method="BH"),
         sig   = case_when(p_adj<0.001~"***", p_adj<0.01~"**",
                           p_adj<0.05~"*",   p_adj<0.1~".", TRUE~"ns")) %>%
  arrange(p_adj)
cat(sprintf("%d / %d families significant (BH p<0.05)\n",
            sum(famab$p_adj<0.05), nrow(famab)))
print(famab %>% filter(p_adj<0.05), n=Inf)
write.csv(famab, "famabund_wilcoxon_ACvsC_onlyACgroups.csv", row.names=FALSE)
cat("Saved famabund_wilcoxon_ACvsC_onlyACgroups.csv\n")

# =========================================================
# 3. RANDOM FOREST: binary GCtS vs EtS
# =========================================================
cat("\n=== 3. Random Forest: GCtS vs EtS (1000 trees) ===\n")
rf <- randomForest(x=X, y=dat$group2, ntree=1000, importance=TRUE)
cat("OOB confusion matrix:\n"); print(rf$confusion)
cat(sprintf("OOB error: %.2f%%\n", rf$err.rate[1000,"OOB"]*100))

# Convergence check
cat("\nConvergence:\n")
for (n in c(100,250,500,750,1000)) {
  cat(sprintf("  ntrees=%4d  OOB=%.3f  GCtS=%.3f  EtS=%.3f\n", n,
              rf$err.rate[n,"OOB"],
              rf$err.rate[n,"GCtS"],
              rf$err.rate[n,"EtS"]))
}

# Stability across seeds
stab <- map_dfr(c(1,7,42,99,123,256,314,500,777,999), function(s) {
  set.seed(s)
  rf_s <- randomForest(x=X, y=dat$group2, ntree=1000, importance=FALSE)
  tibble(seed=s, OOB=rf_s$err.rate[1000,"OOB"]*100,
         GCtS_err=rf_s$err.rate[1000,"GCtS"]*100,
         EtS_err =rf_s$err.rate[1000,"EtS"]*100)
})
cat(sprintf("\nStability (10 seeds): OOB mean=%.2f%%  SD=%.2f%%\n",
            mean(stab$OOB), sd(stab$OOB)))

# Variable importance
imp <- importance(rf, type=1) %>% as.data.frame() %>%
  rownames_to_column("family") %>%
  rename(MDA=MeanDecreaseAccuracy) %>%
  arrange(desc(MDA))

env_row     <- contrib %>% filter(topology=="Environment") %>%
               select(all_of(fam_cols)) %>% mutate(across(everything(),as.numeric))
env_shared  <- shorten_fam(fam_cols[as.numeric(env_row[1,])>0])
imp <- imp %>% mutate(source=ifelse(family %in% env_shared,"env-shared","animal-only"))

cat("\nTop 10 predictors:\n")
print(imp %>% slice_head(n=10))
write.csv(imp, "rf_importance_onlyACgroups.csv", row.names=FALSE)
cat("Saved rf_importance_onlyACgroups.csv\n")

top8 <- imp$family[1:8]
top4 <- imp$family[1:4]

# =========================================================
# 4. ROC CURVE: binary (single clean curve)
# =========================================================
probs  <- predict(rf, type="prob")[,"GCtS"]
roc_obj <- roc(ifelse(dat$group2=="GCtS",1,0), probs, quiet=TRUE)

cat(sprintf("\nAUC = %.3f\n", as.numeric(auc(roc_obj))))

# Bootstrap CI for AUC
set.seed(42)
ci_auc <- ci.auc(roc_obj, method="bootstrap", boot.n=2000,
                 conf.level=0.95, progress="none")
cat(sprintf("95%% CI (bootstrap, n=2000): %.3f - %.3f\n", ci_auc[1], ci_auc[3]))

# Bootstrap CI stability
ci_stab <- map_dfr(c(500,1000,2000), function(nb) {
  set.seed(42)
  ci <- ci.auc(roc_obj, method="bootstrap", boot.n=nb,
               conf.level=0.95, progress="none")
  tibble(boot_n=nb, auc=round(as.numeric(auc(roc_obj)),4),
         ci_lo=round(ci[1],4), ci_hi=round(ci[3],4),
         width=round(ci[3]-ci[1],4))
})
cat("\nBootstrap CI stability:\n"); print(ci_stab)

# Youden's J
best  <- coords(roc_obj, "best", ret=c("threshold","sensitivity","specificity","youden"),
                best.method="youden", transpose=FALSE)
J     <- best$sensitivity + best$specificity - 1

set.seed(42)
n_obs <- length(dat$group2)
resp_bin <- ifelse(dat$group2=="GCtS",1,0)
J_boot <- replicate(2000, {
  idx <- sample(n_obs, n_obs, replace=TRUE)
  r_b <- tryCatch(roc(resp_bin[idx], probs[idx], quiet=TRUE), error=function(e) NULL)
  if (is.null(r_b)) return(NA_real_)
  b   <- coords(r_b,"best",ret=c("sensitivity","specificity"),
                best.method="youden",transpose=FALSE)
  b$sensitivity[1]+b$specificity[1]-1
})
J_ci <- quantile(J_boot, c(0.025,0.975), na.rm=TRUE)

pred_class <- ifelse(probs >= best$threshold, 1, 0)
tp <- sum(pred_class==1 & resp_bin==1); fp <- sum(pred_class==1 & resp_bin==0)
tn <- sum(pred_class==0 & resp_bin==0); fn <- sum(pred_class==0 & resp_bin==1)

roc_metrics <- tibble(
  AUC         = round(as.numeric(auc(roc_obj)),3),
  AUC_lo      = round(ci_auc[1],3),
  AUC_hi      = round(ci_auc[3],3),
  threshold   = round(best$threshold,3),
  sensitivity = round(best$sensitivity,3),
  specificity = round(best$specificity,3),
  Youden_J    = round(J,3),
  Youden_lo   = round(J_ci[1],3),
  Youden_hi   = round(J_ci[2],3),
  PPV         = round(tp/(tp+fp),3),
  NPV         = round(tn/(tn+fn),3),
  F1          = round(2*tp/(2*tp+fp+fn),3)
)
cat("\n=== ROC metrics at optimal threshold ===\n")
print(as.data.frame(roc_metrics))
write.csv(roc_metrics, "roc_metrics_onlyACgroups.csv", row.names=FALSE)
cat("Saved roc_metrics_onlyACgroups.csv\n")

# Smoothed ROC plot (binormal; falls back to raw if near-perfect separation)
roc_smooth <- tryCatch(
  smooth(roc_obj, method = "binormal"),
  error = function(e) {
    cat("Binormal smoothing not applicable (near-perfect separation); using raw ROC.\n")
    roc_obj
  }
)
roc_df <- tibble(spec=roc_smooth$specificities, sens=roc_smooth$sensitivities)

p_roc <- ggplot(roc_df, aes(x=1-spec, y=sens)) +
  geom_line(colour="#2ca02c", linewidth=1.2) +
  geom_abline(linetype="dashed", colour="grey60") +
  annotate("point", x=1-best$specificity, y=best$sensitivity,
           colour="black", size=3, shape=18) +
  annotate("text", x=1-best$specificity+0.03, y=best$sensitivity-0.03,
           label=sprintf("Youden J=%.3f\n(95%% CI: %.3f-%.3f)",
                         J, J_ci[1], J_ci[2]),
           size=3, hjust=0) +
  theme_classic(base_size=11) +
  theme(plot.title=element_text(size=10,face="bold"),
        plot.subtitle=element_text(size=8)) +
  labs(title="ROC: ground-contact (GCtS) vs elevated (EtS) - binary RF",
       subtitle=sprintf("AUC=%.3f (95%% CI: %.3f-%.3f) | raw ROC | bootstrap n=2000",
                        as.numeric(auc(roc_obj)), ci_auc[1], ci_auc[3]),
       x="1 - Specificity", y="Sensitivity")
ggsave("roc_smoothed_onlyACgroups.png", p_roc, width=6, height=5, dpi=300)
ggsave("roc_smoothed_onlyACgroups.pdf", p_roc, width=6, height=5)
cat("Saved roc_smoothed_onlyACgroups.pdf/.png\n")

# =========================================================
# 5. PCA ON TOP PREDICTORS
# =========================================================
pca_plot2 <- function(X_sub, colour_var, colour_vals, colour_label, title_str, stem) {
  pca <- prcomp(X_sub, center=TRUE, scale.=FALSE)
  imp <- summary(pca)$importance
  df  <- as.data.frame(pca$x[,1:2]) %>% mutate(colour=colour_var)
  p <- ggplot(df, aes(x=PC1, y=PC2, colour=colour)) +
    geom_point(size=3, alpha=0.8) +
    scale_colour_manual(values=colour_vals, name=colour_label) +
    theme_classic(base_size=11) +
    labs(title    = title_str,
         subtitle = sprintf("%d families | PC1=%.1f%% | PC2=%.1f%%",
                            ncol(X_sub), imp[2,1]*100, imp[2,2]*100),
         x=sprintf("PC1 (%.1f%%)",imp[2,1]*100),
         y=sprintf("PC2 (%.1f%%)",imp[2,2]*100))
  ggsave(paste0(stem,".pdf"), p, width=7, height=5)
  ggsave(paste0(stem,".png"), p, width=7, height=5, dpi=300)
  cat("Saved", stem, ".pdf/.png\n"); p
}

X4 <- X %>% select(all_of(top4))
X8 <- X %>% select(all_of(top8))

p4_ac  <- pca_plot2(X4, as.character(dat$group2), group2_cols,
                    "Contact group",
                    "PCA top 4 predictors: GCtS vs EtS groups",
                    "pca_top4_ACgroups_onlyACgroups")

p8_ac  <- pca_plot2(X8, as.character(dat$group2), group2_cols,
                    "Contact group",
                    "PCA top 8 predictors: GCtS vs EtS groups",
                    "pca_top8_ACgroups_onlyACgroups")

p8_site <- pca_plot2(X8, dat$topology, topology_cols,
                     "Body site",
                     "PCA top 8 predictors: individual body sites",
                     "pca_top8_sites_onlyACgroups")

# =========================================================
# 6. IMPORTANCE PLOT
# =========================================================
imp_plot <- imp %>% slice_head(n=20) %>%
  mutate(family=factor(family, levels=rev(family)))

src_cols <- c("animal-only"="#d62728", "env-shared"="#7f7f7f")

p_imp <- ggplot(imp_plot, aes(x=MDA, y=family, fill=source)) +
  geom_col(width=0.7, alpha=0.85) +
  scale_fill_manual(values=src_cols,
                    labels=c("animal-only"="Absent from environment",
                             "env-shared" ="Present in environment"),
                    name=NULL) +
  theme_classic(base_size=11) +
  theme(legend.position=c(0.68,0.12), legend.text=element_text(size=9),
        plot.title=element_text(size=10,face="bold"),
        plot.subtitle=element_text(size=8)) +
  labs(title=sprintf("Top 20 family predictors: GCtS vs EtS (binary RF)"),
       subtitle=sprintf("OOB accuracy=%.1f%% | seed stability SD=%.2f%%",
                        (1-rf$err.rate[1000,"OOB"])*100, sd(stab$OOB)),
       x="Mean Decrease Accuracy", y=NULL)
ggsave("rf_importance_onlyACgroups.pdf", p_imp, width=8, height=6)
ggsave("rf_importance_onlyACgroups.png", p_imp, width=8, height=6, dpi=300)
cat("Saved rf_importance_onlyACgroups.pdf/.png\n")

# =========================================================
# 7. DIAGNOSTICS FIGURE (convergence + stability + calibration)
# =========================================================
oob_df <- as.data.frame(rf$err.rate) %>%
  mutate(trees=row_number()) %>%
  pivot_longer(-trees, names_to="class", values_to="error")

p_conv <- ggplot(oob_df, aes(x=trees, y=error*100, colour=class)) +
  geom_line(linewidth=0.6, alpha=0.85) +
  scale_colour_manual(values=c("OOB"="black","GCtS"="#2ca02c",
                                "EtS"="#d62728"),
                      labels=c("OOB"="OOB (overall)","GCtS"="GCtS","EtS"="EtS"),
                      name=NULL) +
  theme_classic(base_size=10) +
  theme(legend.position=c(0.75,0.75), legend.text=element_text(size=8)) +
  labs(title="RF convergence: OOB error vs trees",
       x="Number of trees", y="OOB error (%)")

p_stab <- stab %>%
  pivot_longer(-seed, names_to="metric", values_to="error") %>%
  mutate(metric=factor(metric,levels=c("OOB","GCtS_err","EtS_err"))) %>%
  ggplot(aes(x=factor(seed), y=error, colour=metric, group=metric)) +
  geom_line(linewidth=0.7) + geom_point(size=2) +
  scale_colour_manual(values=c("OOB"="black","GCtS_err"="#2ca02c","EtS_err"="#d62728"),
                      labels=c("OOB (overall)","GCtS","EtS"), name=NULL) +
  theme_classic(base_size=10) +
  theme(legend.position=c(0.75,0.75), legend.text=element_text(size=8)) +
  labs(title=sprintf("RF stability: 10 seeds | SD=%.2f%%", sd(stab$OOB)),
       x="Random seed", y="OOB error (%)")

# Calibration
cal_df <- tibble(obs=resp_bin, pred=probs) %>%
  mutate(bin=cut(pred, breaks=seq(0,1,by=0.1), include.lowest=TRUE)) %>%
  group_by(bin) %>%
  summarise(mean_pred=mean(pred), frac_obs=mean(obs), n=n(), .groups="drop")

p_cal <- ggplot(cal_df, aes(x=mean_pred, y=frac_obs)) +
  geom_abline(linetype="dashed", colour="grey60") +
  geom_point(aes(size=n), colour="#2ca02c", alpha=0.8) +
  geom_line(colour="#2ca02c", alpha=0.6) +
  scale_size_continuous(range=c(2,7), name="n") +
  theme_classic(base_size=10) +
  theme(legend.position=c(0.15,0.75)) +
  labs(title="Probability calibration (GCtS vs EtS)",
       x="Mean predicted probability", y="Observed fraction GCtS")

fig_diag <- (p_conv | p_stab | p_cal) +
  plot_annotation(title="RF diagnostics: binary GCtS vs EtS classification",
                  theme=theme(plot.title=element_text(size=11,face="bold")))
ggsave("rf_diagnostics_onlyACgroups.pdf", fig_diag, width=14, height=5)
ggsave("rf_diagnostics_onlyACgroups.png", fig_diag, width=14, height=5, dpi=300)
cat("Saved rf_diagnostics_onlyACgroups.pdf/.png\n")
