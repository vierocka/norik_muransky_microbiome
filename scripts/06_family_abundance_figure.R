library(tidyverse)
library(patchwork)

setwd("/home/veve/Dropbox/kone/qiime2")

# -------------------------
# Load data
# -------------------------
contrib  <- read.csv("pca_family_contributions.csv", check.names = FALSE)
meta     <- read.csv("/home/veve/Dropbox/kone/kraken2_reports/metadata.csv",
                     check.names = FALSE) %>%
            mutate(across(everything(), trimws), age = as.numeric(age))
kw_topo4 <- read.csv("famabund_kw_topo4.csv",      check.names = FALSE)
kw_age2  <- read.csv("famabund_wilcoxon_age2.csv",  check.names = FALSE)
sp_age   <- read.csv("famabund_spearman_age.csv",   check.names = FALSE)

fam_cols <- setdiff(colnames(contrib), c("sample_id", "topology"))

# Shorten long unresolved taxonomy label
shorten_fam <- function(x) {
  x <- sub("d__Bacteria;p__Patescibacteria;c__Saccharimonadia;o__Saccharimonadales;__",
           "Saccharimonadales", x)
  x
}

# Long format with topo4 grouping (environment included)
dat_long <- contrib %>%
  pivot_longer(all_of(fam_cols), names_to = "family", values_to = "abundance") %>%
  mutate(abundance = as.numeric(abundance),
         family    = shorten_fam(family)) %>%
  left_join(meta %>% select(sample_id, age), by = "sample_id") %>%
  mutate(
    age2  = ifelse(age <= 7, "young (<=7)", "old (>=8)") %>%
            factor(levels = c("young (<=7)", "old (>=8)")),
    topo4 = case_when(
      topology %in% c("Left front pastern", "Muzzle")     ~ "A: ground contact",
      topology %in% c("Ventral abdomen",    "Udder")       ~ "B: near ground",
      topology %in% c("Dorsum", "Forehead", "Neck",
                      "Pectoral area")                      ~ "C: elevated",
      topology == "Environment"                             ~ "D: environment"
    ) %>% factor(levels = c("A: ground contact", "B: near ground",
                             "C: elevated",        "D: environment"))
  )

# Family classification: env-shared vs animal-only
env_row   <- dat_long %>% filter(topology == "Environment") %>%
             group_by(family) %>% summarise(env_ab = mean(abundance), .groups = "drop")
env_shared  <- env_row$family[env_row$env_ab >  0]
anim_only_f <- env_row$family[env_row$env_ab == 0]

# Mean abundance per family × topo4
topo4_means <- dat_long %>%
  group_by(family, topo4) %>%
  summarise(mean_ab = mean(abundance) * 100, .groups = "drop")

topo4_cols <- c(
  "A: ground contact" = "#2ca02c",
  "B: near ground"    = "#17becf",
  "C: elevated"       = "#d62728",
  "D: environment"    = "#7f7f7f"
)

# -------------------------
# Panel 1: env-shared families
# -------------------------
ord1 <- topo4_means %>%
  filter(family %in% env_shared) %>%
  group_by(family) %>% summarise(tot = sum(mean_ab), .groups = "drop") %>%
  arrange(tot) %>% pull(family)

p1 <- topo4_means %>%
  filter(family %in% env_shared) %>%
  mutate(family = factor(family, levels = ord1)) %>%
  ggplot(aes(x = mean_ab, y = family, colour = topo4)) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_colour_manual(values = topo4_cols, name = NULL) +
  scale_x_continuous(labels = scales::label_number(suffix = "%")) +
  theme_classic(base_size = 10) +
  theme(legend.position  = "none",
        axis.text.y      = element_text(size = 8),
        plot.title       = element_text(size = 9, face = "bold"),
        plot.subtitle    = element_text(size = 7)) +
  labs(title    = "Families shared with environment",
       subtitle = sprintf("n = %d families", length(env_shared)),
       x = "Mean relative abundance", y = NULL)

# -------------------------
# Panel 2: animal-only families
# -------------------------
ord2 <- topo4_means %>%
  filter(family %in% anim_only_f, topo4 != "D: environment") %>%
  group_by(family) %>% summarise(tot = sum(mean_ab), .groups = "drop") %>%
  arrange(tot) %>% pull(family)

p2 <- topo4_means %>%
  filter(family %in% anim_only_f, topo4 != "D: environment") %>%
  mutate(family = factor(family, levels = ord2)) %>%
  ggplot(aes(x = mean_ab, y = family, colour = topo4)) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_colour_manual(values = topo4_cols, name = NULL,
                      guide  = guide_legend(override.aes = list(size = 3))) +
  scale_x_continuous(labels = scales::label_number(suffix = "%")) +
  theme_classic(base_size = 10) +
  theme(legend.position  = c(0.72, 0.15),
        legend.text      = element_text(size = 7),
        legend.key.size  = unit(0.4, "cm"),
        axis.text.y      = element_text(size = 8),
        plot.title       = element_text(size = 9, face = "bold"),
        plot.subtitle    = element_text(size = 7)) +
  labs(title    = "Families absent from environment",
       subtitle = sprintf("n = %d families | animal body parts only", length(anim_only_f)),
       x = "Mean relative abundance", y = NULL)

# -------------------------
# Panel 3: topo4 significance (A/B/C + env)
# -------------------------
kw4_plot <- kw_topo4 %>%
  mutate(family      = shorten_fam(family),
         log10p      = -log10(p_adj),
         sig         = ifelse(p_adj < 0.05, "p_adj < 0.05", "ns"),
         family      = fct_reorder(family, log10p))

p3 <- ggplot(kw4_plot, aes(x = log10p, y = family, colour = sig)) +
  geom_segment(aes(xend = 0, yend = family), colour = "grey80", linewidth = 0.4) +
  geom_point(size = 2.5) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed",
             colour = "firebrick", linewidth = 0.6) +
  scale_colour_manual(values = c("p_adj < 0.05" = "#1f77b4", "ns" = "grey60"),
                      name = NULL) +
  theme_classic(base_size = 10) +
  theme(legend.position  = c(0.7, 0.1),
        legend.text      = element_text(size = 7),
        axis.text.y      = element_text(size = 7),
        plot.title       = element_text(size = 9, face = "bold"),
        plot.subtitle    = element_text(size = 7)) +
  labs(title    = "Significance: abundance ~ contact zone (A/B/C + environment)",
       subtitle = "Kruskal-Wallis, BH-adjusted | dashed line = p_adj = 0.05",
       x = expression(-log[10](p[adj])), y = NULL)

# -------------------------
# Panel 4: Age non-dependence — 2 facets
# -------------------------
age_plot <- bind_rows(
  sp_age %>%
    mutate(family    = shorten_fam(family),
           statistic = rho,
           test      = "Age as integer\n(Spearman rho)") %>%
    select(family, statistic, p_adj, test),
  kw_age2 %>%
    mutate(family    = shorten_fam(family),
           statistic = -log10(p_adj),
           test      = "Age as category\n(-log10 p_adj, Wilcoxon)") %>%
    select(family, statistic, p_adj, test)
) %>%
  mutate(sig = ifelse(p_adj < 0.05, "p_adj < 0.05", "ns"),
         test = factor(test, levels = c("Age as integer\n(Spearman rho)",
                                        "Age as category\n(-log10 p_adj, Wilcoxon)")))

# Order families by Spearman rho for consistency across both facets
fam_order_age <- age_plot %>%
  filter(test == "Age as integer\n(Spearman ρ)") %>%
  arrange(statistic) %>% pull(family)

age_plot <- age_plot %>% mutate(family = factor(family, levels = fam_order_age))

p4 <- ggplot(age_plot, aes(x = statistic, y = family, colour = sig)) +
  geom_segment(aes(xend = 0, yend = family), colour = "grey85", linewidth = 0.3) +
  geom_point(size = 1.8, alpha = 0.85) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey40") +
  facet_wrap(~test, scales = "free_x", nrow = 1) +
  scale_colour_manual(values = c("p_adj < 0.05" = "firebrick", "ns" = "grey50"),
                      name = NULL) +
  theme_classic(base_size = 10) +
  theme(legend.position  = "none",
        strip.text       = element_text(size = 8, face = "bold"),
        strip.background = element_rect(fill = "grey95", colour = NA),
        axis.text.y      = element_text(size = 6.5),
        axis.text.x      = element_text(size = 8),
        plot.title       = element_text(size = 9, face = "bold"),
        plot.subtitle    = element_text(size = 7)) +
  labs(title    = "No age dependence in family abundance",
       subtitle = "All families non-significant after BH correction",
       x = NULL, y = NULL)

# -------------------------
# Compose 4-panel figure
# -------------------------
fig <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title   = "Family-level relative abundance: topology and age effects",
    caption = "BH: Benjamini-Hochberg FDR | A: Muzzle+Pastern | B: Udder+Ventral abd. | C: Dorsum+Forehead+Neck+Pectoral",
    theme   = theme(plot.title   = element_text(size = 12, face = "bold"),
                    plot.caption = element_text(size = 7, colour = "grey40"))
  ) +
  plot_layout(heights = c(1.1, 1.3))

ggsave("family_abundance_4panel.pdf", fig, width = 16, height = 18)
ggsave("family_abundance_4panel.png", fig, width = 16, height = 18, dpi = 300)
cat("Saved family_abundance_4panel.pdf/.png\n")
