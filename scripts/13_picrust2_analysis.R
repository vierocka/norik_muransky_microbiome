library(tidyverse)
library(patchwork)

# NOTE: Run picrust2_install_run.sh first to generate picrust2_output/

set.seed(42)
setwd("/home/veve/Dropbox/kone/qiime2")

# =========================================================
# CONFIGURATION — families to examine
# =========================================================

# Gut-associated families enriched in AB_lower_contact
gut_families <- c(
  "Lachnospiraceae", "Ruminococcaceae", "Eggerthellaceae",
  "Christensenellaceae", "Hungateiclostridiaceae", "Anaerovoracaceae",
  "Peptostreptococcales-Tissierellales", "Oscillospiraceae",
  "[Eubacterium]_coprostanoligenes_group"
)

# C_elevated-specific (elevated harness sites; Deinococcaceae top differentiator)
c_specific_families <- c("Deinococcaceae", "Micrococcaceae")

# Top RF predictors (binary AC vs C, ranked by MDA)
rf_top_families <- c(
  "Deinococcaceae", "Intrasporangiaceae", "Hungateiclostridiaceae",
  "Christensenellaceae", "uncultured", "Eggerthellaceae",
  "UCG-010", "Micrococcaceae", "Peptostreptococcales-Tissierellales", "Dietziaceae"
)

# Focus families for figures (union of above)
focus_families <- unique(c(gut_families, c_specific_families, rf_top_families))

# =========================================================
# LOAD METADATA
# =========================================================
meta <- read.csv("/home/veve/Dropbox/kone/kraken2_reports/metadata.csv",
                 check.names = FALSE) %>%
        mutate(across(everything(), trimws)) %>%
        filter(topology != "Environment") %>%
        mutate(group2 = case_when(
          topology %in% c("Left front pastern", "Muzzle",
                          "Ventral abdomen",    "Udder")  ~ "AB_lower_contact",
          topology %in% c("Dorsum", "Forehead", "Neck",
                          "Pectoral area")                 ~ "C_elevated"
        ) %>% factor(levels = c("AB_lower_contact", "C_elevated")))

# PICRUSt2 BIOM sample IDs have _SXX suffix; strip to match metadata
strip_suffix <- function(x) sub("_S[0-9]+$", "", x)

group2_cols <- c("AB_lower_contact" = "#2ca02c", "C_elevated" = "#d62728")

# =========================================================
# LOAD TAXONOMY (ASV → family)
# taxonomy_export/taxonomy.tsv: Feature ID \t Taxon \t Confidence
# =========================================================
tax <- read.delim("taxonomy_export/taxonomy.tsv", check.names = FALSE) %>%
  rename(asv = "Feature ID", taxon = "Taxon") %>%
  mutate(
    family = str_extract(taxon, "f__[^;]+") %>%
             sub("f__", "", .) %>%
             str_trim(),
    family = ifelse(is.na(family) | family == "", "unclassified", family)
  ) %>%
  select(asv, family)

cat(sprintf("Taxonomy: %d ASVs, %d unique families\n",
            nrow(tax), n_distinct(tax$family)))

# =========================================================
# HELPER: load unstratified table and return tidy data frame
# rows = pathway/EC, cols = sample IDs
# =========================================================
load_unstrat <- function(path) {
  tbl <- read.delim(path, check.names = FALSE, comment.char = "")
  # First column is pathway/function ID
  id_col <- colnames(tbl)[1]
  tbl %>%
    pivot_longer(-all_of(id_col), names_to = "sample_raw", values_to = "abundance") %>%
    rename(feature = all_of(id_col)) %>%
    mutate(sample_id = strip_suffix(sample_raw)) %>%
    select(feature, sample_id, abundance)
}

# =========================================================
# HELPER: load stratified table
# Format: pathway \t sample|asv \t abundance
# =========================================================
load_strat <- function(path) {
  tbl <- read.delim(path, check.names = FALSE, comment.char = "")
  id_col  <- colnames(tbl)[1]
  tbl %>%
    pivot_longer(-all_of(id_col), names_to = "sample_taxon", values_to = "abundance") %>%
    rename(feature = all_of(id_col)) %>%
    mutate(
      sample_raw = sub("\\|.*", "", sample_taxon),
      asv        = sub(".*\\|", "", sample_taxon),
      sample_id  = strip_suffix(sample_raw)
    ) %>%
    select(feature, sample_id, asv, abundance) %>%
    left_join(tax, by = "asv")
}

# =========================================================
# 1. LOAD PATHWAY ABUNDANCES (unstratified)
# =========================================================
cat("\nLoading pathway abundances (unstratified)...\n")
path_unstrat_file <- "picrust2_output/pathways_out/path_abun_unstrat.tsv"
if (!file.exists(path_unstrat_file)) stop("Run picrust2_install_run.sh first.")

pw_unstrat <- load_unstrat(path_unstrat_file) %>%
  inner_join(meta %>% select(sample_id, group2), by = "sample_id")

cat(sprintf("Loaded: %d pathways x %d samples\n",
            n_distinct(pw_unstrat$feature), n_distinct(pw_unstrat$sample_id)))

# =========================================================
# 2. WILCOXON per pathway: AB vs C, BH-corrected
# =========================================================
cat("Testing pathways AB vs C (Wilcoxon + BH)...\n")

pw_test <- pw_unstrat %>%
  group_by(feature) %>%
  filter(sum(abundance > 0) >= 4) %>%   # present in at least 4 samples
  summarise(
    mean_AB = mean(abundance[group2 == "AB_lower_contact"]),
    mean_C  = mean(abundance[group2 == "C_elevated"]),
    log2FC  = log2((mean_AB + 1e-6) / (mean_C + 1e-6)),
    W       = wilcox.test(abundance[group2 == "AB_lower_contact"],
                          abundance[group2 == "C_elevated"],
                          exact = FALSE)$statistic,
    p       = wilcox.test(abundance[group2 == "AB_lower_contact"],
                          abundance[group2 == "C_elevated"],
                          exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p, method = "BH"),
    sig   = case_when(p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
                      p_adj < 0.05  ~ "*",   p_adj < 0.1  ~ ".", TRUE ~ "ns"),
    direction = ifelse(log2FC > 0, "AB_higher", "C_higher")
  ) %>%
  arrange(p_adj)

n_sig <- sum(pw_test$p_adj < 0.05, na.rm = TRUE)
cat(sprintf("%d / %d pathways significant (BH p<0.05)\n", n_sig, nrow(pw_test)))
cat(sprintf("  AB-higher: %d | C-higher: %d\n",
            sum(pw_test$p_adj < 0.05 & pw_test$direction == "AB_higher", na.rm = TRUE),
            sum(pw_test$p_adj < 0.05 & pw_test$direction == "C_higher",  na.rm = TRUE)))

write.csv(pw_test, "picrust2_pathway_ABvsC.csv", row.names = FALSE)
cat("Saved picrust2_pathway_ABvsC.csv\n")

# =========================================================
# 3. STRATIFIED: per-family pathway contributions
# Sum contributions by family per sample, then test AB vs C
# =========================================================
cat("\nLoading stratified pathway output (may be slow for large files)...\n")
path_strat_file <- "picrust2_output/pathways_out/path_abun_strat.tsv"

pw_strat_raw <- load_strat(path_strat_file)

# Sum by family per sample per pathway
pw_family <- pw_strat_raw %>%
  group_by(feature, sample_id, family) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  inner_join(meta %>% select(sample_id, group2), by = "sample_id")

cat(sprintf("Stratified: %d pathway x family x sample combinations\n", nrow(pw_family)))

# =========================================================
# 4. FOR EACH FOCUS FAMILY: top contributing pathways
# =========================================================
cat("\n=== Top MetaCyc pathways per focus family ===\n")

family_pathway_tbl <- map_dfr(focus_families, function(fam) {
  df <- pw_family %>%
    filter(family == fam) %>%
    group_by(feature, group2) %>%
    summarise(mean_abund = mean(abundance), .groups = "drop") %>%
    group_by(feature) %>%
    summarise(
      mean_AB    = mean_abund[group2 == "AB_lower_contact"],
      mean_C     = mean_abund[group2 == "C_elevated"],
      total_mean = mean(mean_abund),
      .groups = "drop"
    ) %>%
    filter(total_mean > 0) %>%
    arrange(desc(total_mean)) %>%
    slice_head(n = 10) %>%
    mutate(query_family = fam)
  df
})

write.csv(family_pathway_tbl, "picrust2_top_pathways_per_family.csv", row.names = FALSE)
cat("Saved picrust2_top_pathways_per_family.csv\n")

# =========================================================
# 5. WILCOXON per pathway per family: which pathways differ
#    between AB and C within each focus family's contribution
# =========================================================
cat("Testing per-family pathway contributions AB vs C...\n")

family_pw_test <- map_dfr(focus_families, function(fam) {
  df <- pw_family %>%
    filter(family == fam) %>%
    group_by(feature) %>%
    filter(sum(abundance > 0) >= 4) %>%
    summarise(
      log2FC = log2((mean(abundance[group2 == "AB_lower_contact"]) + 1e-9) /
                    (mean(abundance[group2 == "C_elevated"]) + 1e-9)),
      p      = tryCatch(
        wilcox.test(abundance[group2 == "AB_lower_contact"],
                    abundance[group2 == "C_elevated"],
                    exact = FALSE)$p.value,
        error = function(e) NA_real_),
      .groups = "drop"
    ) %>%
    mutate(query_family = fam)
})

family_pw_test <- family_pw_test %>%
  mutate(p_adj = p.adjust(p, method = "BH"),
         sig   = case_when(p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
                           p_adj < 0.05  ~ "*",   TRUE ~ "ns")) %>%
  arrange(query_family, p_adj)

write.csv(family_pw_test, "picrust2_family_pathway_test.csv", row.names = FALSE)
cat("Saved picrust2_family_pathway_test.csv\n")

# =========================================================
# 6. FIGURE A: heatmap of top pathways per group of families
# Shows mean contribution per sample group
# =========================================================

# Top 20 significant pathways overall (AB vs C)
top_pw <- pw_test %>%
  filter(p_adj < 0.05) %>%
  arrange(p_adj) %>%
  slice_head(n = 20) %>%
  pull(feature)

if (length(top_pw) > 0) {
  heat_df <- pw_unstrat %>%
    filter(feature %in% top_pw) %>%
    group_by(feature, group2) %>%
    summarise(mean_abund = mean(abundance), .groups = "drop") %>%
    mutate(log_abund = log10(mean_abund + 1),
           feature   = factor(feature, levels = rev(top_pw)))

  p_heat <- ggplot(heat_df, aes(x = group2, y = feature, fill = log_abund)) +
    geom_tile(colour = "white") +
    scale_fill_gradient2(low = "#2c7bb6", mid = "#ffffbf", high = "#d7191c",
                         midpoint = median(heat_df$log_abund),
                         name = "log10(mean\nabundance+1)") +
    theme_classic(base_size = 9) +
    theme(axis.text.x  = element_text(angle = 20, hjust = 1),
          axis.text.y  = element_text(size = 7),
          plot.title   = element_text(size = 10, face = "bold"),
          plot.subtitle = element_text(size = 8)) +
    labs(title    = "Top 20 MetaCyc pathways: AB_lower_contact vs C_elevated",
         subtitle = sprintf("Wilcoxon + BH correction | %d sig. pathways total", n_sig),
         x = NULL, y = NULL)

  ggsave("picrust2_pathway_heatmap.pdf", p_heat, width = 6.5, height = 7)
  ggsave("picrust2_pathway_heatmap.png", p_heat, width = 6.5, height = 7, dpi = 300)
  cat("Saved picrust2_pathway_heatmap.pdf/.png\n")
}

# =========================================================
# 7. FIGURE B: gut families vs Deinococcaceae — top pathways
# Lollipop of log2FC (AB vs C) for their top contributed pathways
# =========================================================

plot_family_lollipop <- function(fam, col_val, n_top = 12) {
  df <- family_pw_test %>%
    filter(query_family == fam, !is.na(p)) %>%
    inner_join(
      pw_family %>%
        filter(family == fam) %>%
        group_by(feature) %>%
        summarise(total_mean = mean(abundance), .groups = "drop"),
      by = "feature"
    ) %>%
    arrange(desc(total_mean)) %>%
    slice_head(n = n_top) %>%
    mutate(feature = factor(feature, levels = rev(feature)))

  if (nrow(df) == 0) return(NULL)

  ggplot(df, aes(x = log2FC, y = feature)) +
    geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed") +
    geom_segment(aes(x = 0, xend = log2FC, yend = feature),
                 colour = col_val, linewidth = 0.8) +
    geom_point(aes(colour = sig, size = -log10(p + 1e-10))) +
    scale_colour_manual(values = c("***" = "#d62728", "**" = "#ff7f0e",
                                   "*" = "#1f77b4", "ns" = "grey60"),
                        name = "BH adj.p") +
    scale_size_continuous(range = c(2, 6), name = "-log10(p)") +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(size = 9, face = "bold")) +
    labs(title = fam,
         x = "log2FC (AB / C)", y = NULL)
}

# Panel: gut families
gut_plots <- map(gut_families[gut_families %in% unique(pw_family$family)],
                 plot_family_lollipop, col_val = "#2ca02c")
gut_plots <- compact(gut_plots)

if (length(gut_plots) >= 2) {
  n_col <- min(3, length(gut_plots))
  fig_gut <- wrap_plots(gut_plots, ncol = n_col) +
    plot_annotation(
      title = "MetaCyc pathway contributions: gut-associated families (AB enriched)",
      subtitle = "log2FC > 0 = higher in AB_lower_contact | log2FC < 0 = higher in C_elevated",
      theme = theme(plot.title    = element_text(size = 11, face = "bold"),
                    plot.subtitle = element_text(size = 8))
    )
  ggsave("picrust2_gut_families_pathways.pdf", fig_gut,
         width = 5 * n_col, height = 4 * ceiling(length(gut_plots) / n_col))
  ggsave("picrust2_gut_families_pathways.png", fig_gut,
         width = 5 * n_col, height = 4 * ceiling(length(gut_plots) / n_col), dpi = 300)
  cat("Saved picrust2_gut_families_pathways.pdf/.png\n")
}

# Panel: C-specific (Deinococcaceae + Micrococcaceae)
c_plots <- map(c_specific_families[c_specific_families %in% unique(pw_family$family)],
               plot_family_lollipop, col_val = "#d62728")
c_plots <- compact(c_plots)
if (length(c_plots) >= 1) {
  fig_c <- wrap_plots(c_plots, ncol = min(2, length(c_plots))) +
    plot_annotation(
      title = "MetaCyc pathway contributions: C_elevated-specific families",
      theme = theme(plot.title = element_text(size = 11, face = "bold"))
    )
  ggsave("picrust2_C_specific_pathways.pdf", fig_c, width = 10, height = 5)
  ggsave("picrust2_C_specific_pathways.png", fig_c, width = 10, height = 5, dpi = 300)
  cat("Saved picrust2_C_specific_pathways.pdf/.png\n")
}

# Panel: top RF predictors
rf_plots <- map(rf_top_families[rf_top_families %in% unique(pw_family$family)],
                plot_family_lollipop, col_val = "#9467bd")
rf_plots <- compact(rf_plots)
if (length(rf_plots) >= 2) {
  fig_rf <- wrap_plots(rf_plots, ncol = min(3, length(rf_plots))) +
    plot_annotation(
      title = "MetaCyc pathway contributions: top RF predictor families",
      theme = theme(plot.title = element_text(size = 11, face = "bold"))
    )
  ggsave("picrust2_RF_top_pathways.pdf", fig_rf,
         width = 5 * min(3, length(rf_plots)),
         height = 4 * ceiling(length(rf_plots) / 3))
  ggsave("picrust2_RF_top_pathways.png", fig_rf,
         width = 5 * min(3, length(rf_plots)),
         height = 4 * ceiling(length(rf_plots) / 3), dpi = 300)
  cat("Saved picrust2_RF_top_pathways.pdf/.png\n")
}

# =========================================================
# 8. LOAD KO STRATIFIED: check specific functional categories
# KO annotations: look for UV repair (K03723 recA etc.),
# SCFA fermentation (K00929 butanoate kinase etc.)
# =========================================================
cat("\n=== KO-level analysis: UV resistance and SCFA markers ===\n")

ko_strat_file <- "picrust2_output/KO_metagenome_out/pred_metagenome_strat.tsv"
ko_strat <- load_strat(ko_strat_file)

# Key KO markers
uv_kos  <- c("K03701", "K03723", "K06925", "K10563")  # recA, uvrA/B, etc.
scfa_kos <- c(
  "K00929",  # butyrate kinase (butyrate)
  "K01034",  # acetate CoA transferase (butyrate)
  "K00656",  # formate acetyltransferase
  "K00634",  # phosphate acetyltransferase (acetate)
  "K00925",  # acetate kinase
  "K01026",  # propionate CoA transferase (propionate)
  "K02030"   # propionate kinase
)

ko_summary <- function(kos, label, fams) {
  df <- ko_strat %>%
    filter(feature %in% kos, family %in% fams) %>%
    group_by(feature, sample_id, family) %>%
    summarise(abundance = sum(abundance), .groups = "drop") %>%
    inner_join(meta %>% select(sample_id, group2), by = "sample_id") %>%
    group_by(feature, family, group2) %>%
    summarise(mean_abund = mean(abundance), .groups = "drop")
  cat(sprintf("\n--- %s (%d KOs) ---\n", label, length(kos)))
  if (nrow(df) > 0) print(as.data.frame(df))
  else cat("  None of these KOs found in these families.\n")
  df
}

uv_df  <- ko_summary(uv_kos,  "UV/DNA-repair KOs", c_specific_families)
scfa_df <- ko_summary(scfa_kos, "SCFA fermentation KOs", gut_families)

# Save combined KO summary
ko_combined <- bind_rows(
  uv_df  %>% mutate(category = "UV/DNA-repair"),
  scfa_df %>% mutate(category = "SCFA-fermentation")
)
write.csv(ko_combined, "picrust2_KO_functional_summary.csv", row.names = FALSE)
cat("\nSaved picrust2_KO_functional_summary.csv\n")

# =========================================================
# 9. ANIMAL-ONLY FAMILIES: do they have distinct functional profiles?
# =========================================================
cat("\n=== Animal-only families: functional enrichment vs env-shared ===\n")

# Load env/animal-only classification from RF results
imp_tbl <- read.csv("rf_importance_onlyACgroups.csv", check.names = FALSE)
anim_only <- imp_tbl$family[imp_tbl$source == "animal-only"]
env_shared <- imp_tbl$family[imp_tbl$source == "env-shared"]

cat(sprintf("Animal-only: %d families | Env-shared: %d families (from RF top 46)\n",
            length(anim_only), length(env_shared)))

# Mean pathway contribution: animal-only families vs env-shared families
pw_anim_vs_env <- pw_family %>%
  mutate(fam_class = case_when(
    family %in% anim_only  ~ "animal-only",
    family %in% env_shared ~ "env-shared",
    TRUE                   ~ "other"
  )) %>%
  filter(fam_class != "other") %>%
  group_by(feature, sample_id, fam_class) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  inner_join(meta %>% select(sample_id, group2), by = "sample_id") %>%
  group_by(feature, fam_class) %>%
  filter(sum(abundance > 0) >= 4) %>%
  summarise(
    mean_abund = mean(abundance),
    p = tryCatch(
      wilcox.test(abundance[fam_class == "animal-only"],
                  abundance[fam_class == "env-shared"],
                  exact = FALSE)$p.value,
      error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  distinct(feature, .keep_all = TRUE) %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  arrange(p_adj)

n_anim_sig <- sum(pw_anim_vs_env$p_adj < 0.05, na.rm = TRUE)
cat(sprintf("%d pathways differ between animal-only vs env-shared families (BH p<0.05)\n",
            n_anim_sig))
write.csv(pw_anim_vs_env, "picrust2_animal_vs_envshared_pathways.csv", row.names = FALSE)
cat("Saved picrust2_animal_vs_envshared_pathways.csv\n")

# =========================================================
# 10. PRINT FINAL SUMMARY
# =========================================================
cat("\n========== PICRUST2 ANALYSIS SUMMARY ==========\n")
cat(sprintf("Pathways tested (AB vs C): %d | Significant: %d\n",
            nrow(pw_test), n_sig))
cat(sprintf("AB_lower_contact-enriched: %d | C_elevated-enriched: %d\n",
            sum(pw_test$p_adj < 0.05 & pw_test$direction == "AB_higher", na.rm = TRUE),
            sum(pw_test$p_adj < 0.05 & pw_test$direction == "C_higher",  na.rm = TRUE)))
cat("\nTop 5 AB-enriched pathways:\n")
pw_test %>% filter(direction == "AB_higher", p_adj < 0.05) %>%
  slice_head(n = 5) %>% select(feature, log2FC, p_adj, sig) %>%
  print()
cat("\nTop 5 C-enriched pathways:\n")
pw_test %>% filter(direction == "C_higher", p_adj < 0.05) %>%
  slice_head(n = 5) %>% select(feature, log2FC, p_adj, sig) %>%
  print()
cat("\nOutput files:\n")
cat("  picrust2_pathway_ABvsC.csv\n")
cat("  picrust2_top_pathways_per_family.csv\n")
cat("  picrust2_family_pathway_test.csv\n")
cat("  picrust2_KO_functional_summary.csv\n")
cat("  picrust2_animal_vs_envshared_pathways.csv\n")
cat("  picrust2_pathway_heatmap.pdf/.png\n")
cat("  picrust2_gut_families_pathways.pdf/.png\n")
cat("  picrust2_C_specific_pathways.pdf/.png\n")
cat("  picrust2_RF_top_pathways.pdf/.png\n")
