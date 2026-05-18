library(tidyverse)

setwd("/home/veve/Dropbox/kone/qiime2")

contrib <- read.csv("pca_family_contributions.csv", check.names = FALSE)
meta    <- read.csv("/home/veve/Dropbox/kone/kraken2_reports/metadata.csv",
                    check.names = FALSE) %>%
           mutate(across(everything(), trimws), age = as.numeric(age))

fam_cols <- setdiff(colnames(contrib), c("sample_id", "topology"))

# Long format: one row per sample × family (environment included)
dat_long <- contrib %>%
  pivot_longer(all_of(fam_cols), names_to = "family", values_to = "abundance") %>%
  mutate(abundance = as.numeric(abundance)) %>%
  left_join(meta %>% select(sample_id, age), by = "sample_id") %>%
  mutate(
    age2  = ifelse(age <= 7, "young (<=7)", "old (>=8)") %>%
            factor(levels = c("young (<=7)", "old (>=8)")),
    age3  = case_when(age  < 7  ~ "young (<7)",
                      age <= 13 ~ "mid (7-13)",
                      TRUE      ~ "old (>13)") %>%
            factor(levels = c("young (<7)", "mid (7-13)", "old (>13)")),
    topo4 = case_when(
      topology %in% c("Left front pastern", "Muzzle")     ~ "A_ground_contact",
      topology %in% c("Ventral abdomen",    "Udder")       ~ "B_near_ground",
      topology %in% c("Dorsum", "Forehead", "Neck",
                      "Pectoral area")                      ~ "C_elevated",
      topology == "Environment"                             ~ "D_environment"
    ) %>% factor(levels = c("A_ground_contact", "B_near_ground",
                             "C_elevated", "D_environment")),
    topology = factor(topology)
  )

# Age-based tests use only animal samples
dat_animal <- dat_long %>% filter(topology != "Environment")

families <- unique(dat_long$family)

# -------------------------
# Helper functions
# -------------------------
kw_test <- function(df, group_var) {
  map_dfr(families, function(f) {
    d  <- df %>% filter(family == f)
    kt <- kruskal.test(d$abundance ~ d[[group_var]])
    tibble(family = f, chi2 = round(kt$statistic, 3),
           df_kw  = kt$parameter, p = kt$p.value)
  }) %>%
    mutate(p_adj = p.adjust(p, method = "BH"),
           sig   = case_when(p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
                             p_adj < 0.05  ~ "*",   p_adj < 0.1  ~ ".",
                             TRUE ~ "ns")) %>%
    arrange(p_adj)
}

spearman_test <- function(df) {
  map_dfr(families, function(f) {
    d  <- df %>% filter(family == f)
    ct <- cor.test(d$abundance, d$age, method = "spearman", exact = FALSE)
    tibble(family = f, rho = round(ct$estimate, 3), p = ct$p.value)
  }) %>%
    mutate(p_adj = p.adjust(p, method = "BH"),
           sig   = case_when(p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
                             p_adj < 0.05  ~ "*",   p_adj < 0.1  ~ ".",
                             TRUE ~ "ns")) %>%
    arrange(p_adj)
}

# -------------------------
# 1. Topology: 8 animal sites + environment (9 levels)
# -------------------------
cat("=== 1. Kruskal-Wallis: abundance ~ topology (8 sites + environment) ===\n")
kw_topo8 <- kw_test(dat_long, "topology")
print(kw_topo8, n = Inf)
write.csv(kw_topo8, "famabund_kw_topology8.csv", row.names = FALSE)
cat("Saved famabund_kw_topology8.csv\n")

sig_topo8 <- kw_topo8 %>% filter(p_adj < 0.05) %>% pull(family)
if (length(sig_topo8) > 0) {
  cat(sprintf("\nPost-hoc pairwise Wilcoxon for %d significant families:\n",
              length(sig_topo8)))
  pw_topo8 <- map_dfr(sig_topo8, function(f) {
    d  <- dat_long %>% filter(family == f)
    pw <- pairwise.wilcox.test(d$abundance, d$topology,
                               p.adjust.method = "BH", exact = FALSE)
    as.data.frame(pw$p.value) %>%
      rownames_to_column("site1") %>%
      pivot_longer(-site1, names_to = "site2", values_to = "p_adj") %>%
      filter(!is.na(p_adj)) %>%
      mutate(family = f, sig = case_when(p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
                                         p_adj < 0.05  ~ "*",   TRUE ~ "ns"))
  })
  write.csv(pw_topo8, "famabund_posthoc_topology8.csv", row.names = FALSE)
  cat("Saved famabund_posthoc_topology8.csv\n")
}

# -------------------------
# 2. Topology: 4 groups (A/B/C + environment)
# -------------------------
cat("\n=== 2. Kruskal-Wallis: abundance ~ topo4 (A/B/C + environment) ===\n")
kw_topo4 <- kw_test(dat_long, "topo4")
print(kw_topo4, n = Inf)
write.csv(kw_topo4, "famabund_kw_topo4.csv", row.names = FALSE)
cat("Saved famabund_kw_topo4.csv\n")

sig_topo4 <- kw_topo4 %>% filter(p_adj < 0.05) %>% pull(family)
if (length(sig_topo4) > 0) {
  cat(sprintf("\nPost-hoc pairwise Wilcoxon for %d significant families:\n",
              length(sig_topo4)))
  pw_topo4 <- map_dfr(sig_topo4, function(f) {
    d  <- dat_long %>% filter(family == f)
    pw <- pairwise.wilcox.test(d$abundance, d$topo4,
                               p.adjust.method = "BH", exact = FALSE)
    as.data.frame(pw$p.value) %>%
      rownames_to_column("group1") %>%
      pivot_longer(-group1, names_to = "group2", values_to = "p_adj") %>%
      filter(!is.na(p_adj)) %>%
      mutate(family = f, sig = case_when(p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
                                         p_adj < 0.05  ~ "*",   TRUE ~ "ns"))
  })
  write.csv(pw_topo4, "famabund_posthoc_topo4.csv", row.names = FALSE)
  cat("Saved famabund_posthoc_topo4.csv\n")
}

# -------------------------
# 3. Age: 2 groups (animal samples only)
# -------------------------
cat("\n=== 3. Wilcoxon: abundance ~ age2 (young <=7 vs old >=8) ===\n")
kw_age2 <- map_dfr(families, function(f) {
  d     <- dat_animal %>% filter(family == f)
  young <- d$abundance[d$age2 == "young (<=7)"]
  old   <- d$abundance[d$age2 == "old (>=8)"]
  wt    <- wilcox.test(young, old, exact = FALSE)
  tibble(family    = f,
         W         = round(wt$statistic, 1),
         med_young = round(median(young), 4),
         med_old   = round(median(old),   4),
         p         = wt$p.value)
}) %>%
  mutate(p_adj = p.adjust(p, method = "BH"),
         sig   = case_when(p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
                           p_adj < 0.05  ~ "*",   p_adj < 0.1  ~ ".",
                           TRUE ~ "ns")) %>%
  arrange(p_adj)
print(kw_age2, n = Inf)
write.csv(kw_age2, "famabund_wilcoxon_age2.csv", row.names = FALSE)
cat("Saved famabund_wilcoxon_age2.csv\n")

# -------------------------
# 4. Age: 3 groups (animal samples only)
# -------------------------
cat("\n=== 4. Kruskal-Wallis: abundance ~ age3 (3 groups) ===\n")
kw_age3 <- kw_test(dat_animal, "age3")
print(kw_age3, n = Inf)
write.csv(kw_age3, "famabund_kw_age3.csv", row.names = FALSE)
cat("Saved famabund_kw_age3.csv\n")

# -------------------------
# 5. Age: continuous Spearman (animal samples only)
# -------------------------
cat("\n=== 5. Spearman: abundance ~ age (continuous) ===\n")
sp_age <- spearman_test(dat_animal)
print(sp_age, n = Inf)
write.csv(sp_age, "famabund_spearman_age.csv", row.names = FALSE)
cat("Saved famabund_spearman_age.csv\n")

# -------------------------
# Summary: combined significance across all tests
# -------------------------
summary_tbl <- tibble(family = families) %>%
  left_join(kw_topo8 %>% select(family, p_adj) %>% rename(p_topo8 = p_adj), by = "family") %>%
  left_join(kw_topo4 %>% select(family, p_adj) %>% rename(p_topo4 = p_adj), by = "family") %>%
  left_join(kw_age2  %>% select(family, p_adj) %>% rename(p_age2  = p_adj), by = "family") %>%
  left_join(kw_age3  %>% select(family, p_adj) %>% rename(p_age3  = p_adj), by = "family") %>%
  left_join(sp_age   %>% select(family, p_adj) %>% rename(p_spear = p_adj), by = "family") %>%
  mutate(across(starts_with("p_"), \(x) round(x, 4))) %>%
  arrange(p_topo8)

cat("\n=== Summary: BH-adjusted p across all tests ===\n")
cat("Columns: p_topo8=9 levels (8 sites+env) | p_topo4=A/B/C+env | p_age2=young/old | p_age3=3 groups | p_spear=Spearman\n")
print(summary_tbl, n = Inf)
write.csv(summary_tbl, "famabund_summary_all_tests.csv", row.names = FALSE)
cat("Saved famabund_summary_all_tests.csv\n")
