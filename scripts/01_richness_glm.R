library(tidyverse)
library(patchwork)
library(MASS, exclude = "select")   # glm.nb without masking dplyr::select

setwd("/home/veve/Dropbox/kone/qiime2")

# -------------------------
# Data — all 12 horses; outliers identified per tissue by IQR rule
# -------------------------
contrib <- read.csv("pca_family_contributions.csv", check.names = FALSE)
meta    <- read.csv("/home/veve/Dropbox/kone/kraken2_reports/metadata.csv",
                    check.names = FALSE) %>%
           mutate(across(everything(), trimws), age = as.numeric(age))

fam_cols <- setdiff(colnames(contrib), c("sample_id", "topology"))

dat <- contrib %>%
  filter(topology != "Environment") %>%
  mutate(richness = rowSums(across(all_of(fam_cols), as.numeric) > 0)) %>%
  select(sample_id, topology, richness) %>%
  left_join(meta %>% select(sample_id, age), by = "sample_id") %>%
  mutate(horse     = sub("-.*", "", sample_id),
         age_group = ifelse(age <= 7, "young (<=7)", "old (>=8)") %>%
                     factor(levels = c("young (<=7)", "old (>=8)"))) %>%
  # Flag IQR outliers within each topology
  group_by(topology) %>%
  mutate(Q1      = quantile(richness, 0.25),
         Q3      = quantile(richness, 0.75),
         IQR_val = Q3 - Q1,
         outlier = richness < Q1 - 1.5 * IQR_val |
                   richness > Q3 + 1.5 * IQR_val) %>%
  ungroup()

focus_sites <- c("Dorsum", "Forehead", "Left front pastern", "Muzzle",
                 "Neck", "Pectoral area", "Udder", "Ventral abdomen")

# -------------------------
# Per-site GLM: richness ~ age (continuous)
# Step 1: Poisson. Step 2: check dispersion (residual deviance / df).
# Step 3: if dispersion > 1.5 use quasi-Poisson; if > 3 also try NB.
# Step 4: p-value from drop1(test = "F") for quasi; LRT for NB/Poisson.
# -------------------------
fit_site <- function(site_name) {
  d <- dat %>% filter(topology == site_name, !outlier)

  # Poisson baseline
  m_pois <- glm(richness ~ age, data = d, family = poisson)
  disp   <- m_pois$deviance / m_pois$df.residual

  if (disp > 1.5) {
    # Quasi-Poisson
    m_qp  <- glm(richness ~ age, data = d, family = quasipoisson)
    d1    <- drop1(m_qp, test = "F")
    pval  <- d1["age", "Pr(>F)"]
    model <- "quasi-Poisson"
    fit   <- m_qp

    # NB for comparison if very overdispersed
    nb_note <- ""
    if (disp > 3) {
      tryCatch({
        m_nb   <- glm.nb(richness ~ age, data = d)
        nb_p   <- drop1(m_nb, test = "Chisq")["age", "Pr(>Chi)"]
        nb_note <- sprintf(" | NB p=%.3f", nb_p)
      }, error = function(e) {})
    }
    pval_note <- sprintf("%.4f%s", pval, nb_note)
  } else {
    d1    <- drop1(m_pois, test = "LRT")
    pval  <- d1["age", "Pr(>Chi)"]
    model <- "Poisson"
    fit   <- m_pois
    pval_note <- sprintf("%.4f", pval)
  }

  slope  <- coef(fit)["age"]
  se     <- sqrt(vcov(fit)["age", "age"])

  tibble(
    site        = site_name,
    n           = nrow(d),
    mean_rich   = round(mean(d$richness), 1),
    dispersion  = round(disp, 2),
    model       = model,
    slope_log   = round(slope, 4),
    SE          = round(se, 4),
    p_raw       = pval,
    pval_label  = pval_note
  )
}

glm_res <- map_dfr(focus_sites, fit_site) %>%
  mutate(p_adj = p.adjust(p_raw, method = "BH"),
         sig   = case_when(p_adj < 0.001 ~ "***",
                           p_adj < 0.01  ~ "**",
                           p_adj < 0.05  ~ "*",
                           TRUE          ~ "ns"))

cat("=== GLM richness ~ age | focus sites (per-tissue IQR outliers excluded) ===\n")
print(glm_res)
write.csv(glm_res, "richness_age_glm_results.csv", row.names = FALSE)
cat("Saved richness_age_glm_results.csv\n")

# Overdispersion details per site
cat("\n=== Dispersion check (residual deviance / df) ===\n")
for (s in focus_sites) {
  d  <- dat %>% filter(topology == s, !outlier)
  mp <- glm(richness ~ age, data = d, family = poisson)
  cat(sprintf("  %-20s  deviance=%.2f  df=%d  dispersion=%.2f  -> %s\n",
              s, mp$deviance, mp$df.residual,
              mp$deviance / mp$df.residual,
              ifelse(mp$deviance / mp$df.residual > 1.5, "quasi-Poisson", "Poisson")))
}

# -------------------------
# Colour palettes
# -------------------------
age_cols <- c("young (<=7)" = "#4575b4", "old (>=8)" = "#d73027")

topology_cols <- c(
  "Dorsum"             = "#1f77b4", "Forehead"           = "#ff7f0e",
  "Left front pastern" = "#2ca02c", "Muzzle"             = "#9467bd",
  "Neck"               = "#8c564b", "Pectoral area"      = "#7f7f7f",
  "Udder"              = "#bcbd22", "Ventral abdomen"    = "#17becf"
)

# -------------------------
# Plot 1: age_group_richness_tissue (all 8 sites, horse 8 excluded)
# Same style as clustering_age_groups.R
# -------------------------
dat_plot <- dat %>% filter(!outlier)

n_young <- sum(dat_plot$age_group == "young (<=7)" & !duplicated(dat_plot$horse))
n_old   <- sum(dat_plot$age_group == "old (>=8)"   & !duplicated(dat_plot$horse))

# Per-topology n for subtitle
topo_n <- dat_plot %>% count(topology) %>%
  mutate(label = paste0(topology, " (n=", n, ")")) %>%
  pull(label) %>% paste(collapse = ", ")

p_tissue <- dat_plot %>%
  ggplot(aes(x = topology, y = richness, fill = age_group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA,
               position = position_dodge(width = 0.75)) +
  geom_jitter(aes(colour = age_group),
              position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
              size = 1.5, alpha = 0.7) +
  scale_fill_manual(values   = age_cols) +
  scale_colour_manual(values = age_cols) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  labs(title    = "Richness per tissue by age group (per-tissue IQR outliers excluded)",
       subtitle = sprintf("Young <=7 yr | Old >=8 yr | outliers removed per topology"),
       x = NULL, y = "Families present",
       fill = "Age group", colour = "Age group")

ggsave("age_group_richness_tissue_outliers_removed.pdf", p_tissue, width = 10, height = 5)
ggsave("age_group_richness_tissue_outliers_removed.png", p_tissue, width = 10, height = 5, dpi = 300)
cat("Saved age_group_richness_tissue_outliers_removed.pdf/.png\n")

# -------------------------
# Plot 2: richness ~ continuous age per focus site with GLM fit + sig label
# -------------------------
# Build label data for significance annotation
sig_labels <- glm_res %>%
  select(site, model, dispersion, p_adj, sig) %>%
  left_join(
    dat %>% filter(topology %in% focus_sites) %>%
      group_by(topology) %>%
      summarise(x_pos = max(age) * 0.85,
                y_pos = max(richness) * 1.02, .groups = "drop"),
    by = c("site" = "topology")
  ) %>%
  mutate(label = sprintf("%s\ndisp=%.2f\np_adj=%.3f %s",
                         model, dispersion, p_adj, sig))

p_glm <- dat %>%
  filter(topology %in% focus_sites, !outlier) %>%
  mutate(topology = factor(topology, levels = focus_sites)) %>%
  ggplot(aes(x = age, y = richness)) +
  geom_point(aes(colour = age_group), size = 2.5, alpha = 0.8) +
  geom_smooth(method = "glm", method.args = list(family = "poisson"),
              se = TRUE, colour = "grey30", linewidth = 0.8) +
  geom_text(data = sig_labels %>% mutate(topology = factor(site, levels = focus_sites)),
            aes(x = x_pos, y = y_pos, label = label),
            size = 2.8, hjust = 1, vjust = 1, inherit.aes = FALSE) +
  scale_colour_manual(values = age_cols) +
  facet_wrap(~topology, scales = "free_y", nrow = 2) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(title    = "Richness vs continuous age | all body sites",
       subtitle = "GLM Poisson / quasi-Poisson fit | per-tissue IQR outliers excluded | BH-adjusted p",
       x = "Age (years)", y = "Families present", colour = "Age group")

ggsave("richness_age_glm_plot_per_tissue.pdf", p_glm, width = 14, height = 8)
ggsave("richness_age_glm_plot_per_tissue.png", p_glm, width = 14, height = 8, dpi = 300)
cat("Saved richness_age_glm_plot_per_tissue.pdf/.png\n")

# -------------------------
# Udder: Wilcoxon young vs old (per-tissue outliers excluded)
# -------------------------
udder <- dat %>% filter(topology == "Udder", !outlier)

cat("\n=== Udder richness: young (<=7) vs old (>=8) ===\n")
cat(sprintf("n young = %d | n old = %d\n",
            sum(udder$age_group == "young (<=7)"),
            sum(udder$age_group == "old (>=8)")))
cat(sprintf("Median young = %.1f | Median old = %.1f\n",
            median(udder$richness[udder$age_group == "young (<=7)"]),
            median(udder$richness[udder$age_group == "old (>=8)"])))

wt_udder <- wilcox.test(richness ~ age_group, data = udder, exact = FALSE)
cat(sprintf("Wilcoxon W = %.0f | p = %.4f\n", wt_udder$statistic, wt_udder$p.value))

# Boxplot with jitter
p_udder <- udder %>%
  ggplot(aes(x = age_group, y = richness, fill = age_group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, width = 0.5) +
  geom_jitter(aes(colour = age_group), width = 0.1, size = 3, alpha = 0.8) +
  scale_fill_manual(values   = age_cols) +
  scale_colour_manual(values = age_cols) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        plot.title    = element_text(size = 7),
        plot.subtitle = element_text(size = 6)) +
  labs(title    = "Udder richness: young vs old (per-tissue IQR outliers excluded)",
       subtitle = sprintf("Wilcoxon W=%.0f, p=%.4f | quasi-Poisson slope=-0.033, p_raw=0.173",
                          wt_udder$statistic, wt_udder$p.value),
       x = "Age group", y = "Families present")

ggsave("udder_richness_young_vs_old.pdf", p_udder, width = 5, height = 5)
ggsave("udder_richness_young_vs_old.png", p_udder, width = 5, height = 5, dpi = 300)
cat("Saved udder_richness_young_vs_old.pdf/.png\n")

# -------------------------
# Dorsum: Wilcoxon young vs old (per-tissue outliers excluded)
# -------------------------
dorsum <- dat %>% filter(topology == "Dorsum", !outlier)

cat("\n=== Dorsum richness: young (<=7) vs old (>=8) ===\n")
cat(sprintf("n young = %d | n old = %d\n",
            sum(dorsum$age_group == "young (<=7)"),
            sum(dorsum$age_group == "old (>=8)")))
cat(sprintf("Median young = %.1f | Median old = %.1f\n",
            median(dorsum$richness[dorsum$age_group == "young (<=7)"]),
            median(dorsum$richness[dorsum$age_group == "old (>=8)"])))

wt_dorsum <- wilcox.test(richness ~ age_group, data = dorsum, exact = FALSE)
cat(sprintf("Wilcoxon W = %.0f | p = %.4f\n", wt_dorsum$statistic, wt_dorsum$p.value))

p_dorsum <- dorsum %>%
  ggplot(aes(x = age_group, y = richness, fill = age_group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, width = 0.5) +
  geom_jitter(aes(colour = age_group), width = 0.1, size = 3, alpha = 0.8) +
  scale_fill_manual(values   = age_cols) +
  scale_colour_manual(values = age_cols) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        plot.title    = element_text(size = 7),
        plot.subtitle = element_text(size = 6)) +
  labs(title    = "Dorsum richness: young vs old (per-tissue IQR outliers excluded)",
       subtitle = sprintf("Wilcoxon W=%.0f, p=%.4f | Poisson slope=+0.009, p_raw=0.965",
                          wt_dorsum$statistic, wt_dorsum$p.value),
       x = "Age group", y = "Families present")

ggsave("dorsum_richness_young_vs_old.pdf", p_dorsum, width = 5, height = 5)
ggsave("dorsum_richness_young_vs_old.png", p_dorsum, width = 5, height = 5, dpi = 300)
cat("Saved dorsum_richness_young_vs_old.pdf/.png\n")
