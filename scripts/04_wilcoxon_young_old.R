library(tidyverse)

setwd("/home/veve/Dropbox/kone/qiime2")

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
                     factor(levels = c("young (<=7)", "old (>=8)")))

age_cols <- c("young (<=7)" = "#4575b4", "old (>=8)" = "#d73027")

# -------------------------
# Helper: print stats + plot for one tissue, all samples
# -------------------------
test_tissue <- function(tissue_name, slope_label) {
  d <- dat %>% filter(topology == tissue_name)

  cat(sprintf("\n=== %s richness: young (<=7) vs old (>=8) | ALL samples ===\n", tissue_name))
  cat(sprintf("n young = %d | n old = %d\n",
              sum(d$age_group == "young (<=7)"),
              sum(d$age_group == "old (>=8)")))
  cat(sprintf("Median young = %.1f | Median old = %.1f\n",
              median(d$richness[d$age_group == "young (<=7)"]),
              median(d$richness[d$age_group == "old (>=8)"])))

  wt <- wilcox.test(richness ~ age_group, data = d, exact = FALSE)
  cat(sprintf("Wilcoxon W = %.0f | p = %.4f\n", wt$statistic, wt$p.value))

  stem <- tolower(gsub(" ", "_", tissue_name))

  p <- d %>%
    ggplot(aes(x = age_group, y = richness, fill = age_group)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA, width = 0.5) +
    geom_jitter(aes(colour = age_group), width = 0.1, size = 3, alpha = 0.8) +
    scale_fill_manual(values   = age_cols) +
    scale_colour_manual(values = age_cols) +
    theme_classic(base_size = 13) +
    theme(legend.position = "none",
          plot.title    = element_text(size = 7),
          plot.subtitle = element_text(size = 6)) +
    labs(title    = sprintf("%s richness: young vs old (all samples)", tissue_name),
         subtitle = sprintf("Wilcoxon W=%.0f, p=%.4f | %s", wt$statistic, wt$p.value, slope_label),
         x = "Age group", y = "Families present")

  ggsave(sprintf("%s_richness_young_vs_old_allsamples.pdf", stem), p, width = 5, height = 5)
  ggsave(sprintf("%s_richness_young_vs_old_allsamples.png", stem), p, width = 5, height = 5, dpi = 300)
  cat(sprintf("Saved %s_richness_young_vs_old_allsamples.pdf/.png\n", stem))
}

test_tissue("Udder",  "quasi-Poisson slope=-0.033, p_raw=0.173")
test_tissue("Dorsum", "Poisson slope=+0.009, p_raw=0.965")

# -------------------------
# Dorsum: age cut-off at 10 (young <10 vs old >=10), all samples
# -------------------------
dorsum10 <- dat %>%
  filter(topology == "Dorsum") %>%
  mutate(age_group10 = ifelse(age < 10, "young (<10)", "old (>=10)") %>%
                       factor(levels = c("young (<10)", "old (>=10)")))

cat("\n=== Dorsum richness: young (<10) vs old (>=10) | ALL samples ===\n")
cat(sprintf("n young = %d | n old = %d\n",
            sum(dorsum10$age_group10 == "young (<10)"),
            sum(dorsum10$age_group10 == "old (>=10)")))
cat(sprintf("Median young = %.1f | Median old = %.1f\n",
            median(dorsum10$richness[dorsum10$age_group10 == "young (<10)"]),
            median(dorsum10$richness[dorsum10$age_group10 == "old (>=10)"])))
cat("Ages young:", sort(dorsum10$age[dorsum10$age_group10 == "young (<10)"]), "\n")
cat("Ages old  :", sort(dorsum10$age[dorsum10$age_group10 == "old (>=10)"]),  "\n")

wt10 <- wilcox.test(richness ~ age_group10, data = dorsum10, exact = FALSE)
cat(sprintf("Wilcoxon W = %.0f | p = %.4f\n", wt10$statistic, wt10$p.value))

age_cols10 <- c("young (<10)" = "#4575b4", "old (>=10)" = "#d73027")

p_d10 <- dorsum10 %>%
  ggplot(aes(x = age_group10, y = richness, fill = age_group10)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, width = 0.5) +
  geom_jitter(aes(colour = age_group10), width = 0.1, size = 3, alpha = 0.8) +
  scale_fill_manual(values   = age_cols10) +
  scale_colour_manual(values = age_cols10) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        plot.title    = element_text(size = 7),
        plot.subtitle = element_text(size = 6)) +
  labs(title    = "Dorsum richness: young vs old (cut-off age 10, all samples)",
       subtitle = sprintf("Wilcoxon W=%.0f, p=%.4f", wt10$statistic, wt10$p.value),
       x = "Age group", y = "Families present")

ggsave("dorsum_richness_cutoff10_allsamples.pdf", p_d10, width = 5, height = 5)
ggsave("dorsum_richness_cutoff10_allsamples.png", p_d10, width = 5, height = 5, dpi = 300)
cat("Saved dorsum_richness_cutoff10_allsamples.pdf/.png\n")
