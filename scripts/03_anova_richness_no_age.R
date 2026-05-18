library(multcomp, exclude = "select")   # load before tidyverse so dplyr::select wins
library(tidyverse)

setwd("/home/veve/Dropbox/kone/qiime2")

# Age dropped: non-significant in both full models (see anova_richness_topology.R)
# ANOVA 1: p_age=0.459 | ANOVA 2: p_age=0.561

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
  mutate(topology = factor(topology),
         horse    = sub("-.*", "", sample_id),
         topo3    = case_when(
           topology %in% c("Left front pastern", "Muzzle")     ~ "A_ground_contact",
           topology %in% c("Ventral abdomen",    "Udder")       ~ "B_near_ground",
           topology %in% c("Dorsum", "Forehead", "Neck",
                           "Pectoral area")                      ~ "C_elevated"
         ) %>% factor(levels = c("A_ground_contact", "B_near_ground", "C_elevated")))

topology_cols <- c(
  "Dorsum"             = "#1f77b4", "Forehead"           = "#ff7f0e",
  "Left front pastern" = "#2ca02c", "Muzzle"             = "#9467bd",
  "Neck"               = "#8c564b", "Pectoral area"      = "#7f7f7f",
  "Udder"              = "#bcbd22", "Ventral abdomen"    = "#17becf"
)
topo3_cols <- c(
  "A_ground_contact" = "#2ca02c",
  "B_near_ground"    = "#17becf",
  "C_elevated"       = "#d62728"
)
topo3_labels <- c(
  "A_ground_contact" = "A: ground contact\n(Muzzle, Pastern)",
  "B_near_ground"    = "B: near ground\n(Udder, Ventral abd.)",
  "C_elevated"       = "C: elevated\n(Dorsum, Forehead,\nNeck, Pectoral)"
)

# -------------------------
# Model 1: richness ~ topology (8 sites)
# -------------------------
m1 <- aov(richness ~ topology, data = dat)
cat("=== ANOVA 1: richness ~ topology (8 sites) ===\n")
print(summary(m1))

a1_tbl <- broom::tidy(m1) %>%
  mutate(across(where(is.numeric), \(x) round(x, 4)),
         sig = case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**",
                         p.value < 0.05  ~ "*",   TRUE ~ "ns"))
write.csv(a1_tbl, "anova1_topology_noage.csv", row.names = FALSE)
cat("Saved anova1_topology_noage.csv\n")

# Compact letter display
cld1     <- cld(glht(m1, linfct = mcp(topology = "Tukey")),
                level = 0.05, decreasing = TRUE)
letters1 <- data.frame(topology = names(cld1$mcletters$Letters),
                       letter   = trimws(cld1$mcletters$Letters))

topo_sum <- dat %>%
  group_by(topology) %>%
  summarise(med  = median(richness),
            ymax = max(richness),
            mean_rich = round(mean(richness), 1),
            sd_rich   = round(sd(richness), 1),
            n         = n(), .groups = "drop") %>%
  left_join(letters1, by = "topology") %>%
  arrange(desc(med)) %>%
  mutate(topology = factor(topology, levels = topology))

dat <- dat %>%
  mutate(topology = factor(topology, levels = levels(topo_sum$topology)))

cat("\n=== Tukey compact letters (topology) ===\n")
print(topo_sum %>% select(topology, n, mean_rich, sd_rich, letter))
write.csv(topo_sum %>% select(topology, n, mean_rich, sd_rich, letter),
          "anova1_richness_by_topology_summary.csv", row.names = FALSE)
cat("Saved anova1_richness_by_topology_summary.csv\n")

p1 <- dat %>%
  ggplot(aes(x = topology, y = richness, fill = topology)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, width = 0.6) +
  geom_jitter(aes(colour = topology), width = 0.15, size = 1.8, alpha = 0.6) +
  geom_text(data = topo_sum,
            aes(x = topology, y = ymax + 1.5, label = letter),
            inherit.aes = FALSE, size = 4.5, fontface = "bold") +
  scale_fill_manual(values   = topology_cols) +
  scale_colour_manual(values = topology_cols) +
  theme_classic(base_size = 12) +
  theme(axis.text.x   = element_text(angle = 40, hjust = 1),
        legend.position = "none") +
  labs(title    = "Family richness by body site",
       subtitle = sprintf("One-way ANOVA F(%d,%d)=%.2f, p<0.001 | Tukey letters above boxes",
                          summary(m1)[[1]]$Df[1],
                          summary(m1)[[1]]$Df[2],
                          summary(m1)[[1]]$`F value`[1]),
       x = NULL, y = "Families present")

ggsave("anova1_richness_by_topology.pdf", p1, width = 8, height = 5)
ggsave("anova1_richness_by_topology.png", p1, width = 8, height = 5, dpi = 300)
cat("Saved anova1_richness_by_topology.pdf/.png\n")

# -------------------------
# Model 2: richness ~ topo3 (3 groups)
# -------------------------
m2 <- aov(richness ~ topo3, data = dat)
cat("\n=== ANOVA 2: richness ~ reduced topology (3 groups) ===\n")
print(summary(m2))

a2_tbl <- broom::tidy(m2) %>%
  mutate(across(where(is.numeric), \(x) round(x, 4)),
         sig = case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**",
                         p.value < 0.05  ~ "*",   TRUE ~ "ns"))
write.csv(a2_tbl, "anova2_topo3_noage.csv", row.names = FALSE)
cat("Saved anova2_topo3_noage.csv\n")

cat("\n=== Tukey post-hoc: topo3 ===\n")
print(TukeyHSD(m2, which = "topo3"))

cld2     <- cld(glht(m2, linfct = mcp(topo3 = "Tukey")),
                level = 0.05, decreasing = TRUE)
letters2 <- data.frame(topo3  = names(cld2$mcletters$Letters),
                       letter = trimws(cld2$mcletters$Letters))

topo3_sum <- dat %>%
  group_by(topo3) %>%
  summarise(med       = median(richness),
            ymax      = max(richness),
            mean_rich = round(mean(richness), 1),
            sd_rich   = round(sd(richness), 1),
            n         = n(), .groups = "drop") %>%
  left_join(letters2, by = "topo3")

cat("\n=== Mean richness per reduced topology group ===\n")
print(topo3_sum %>% select(topo3, n, mean_rich, sd_rich, letter))

p2 <- dat %>%
  ggplot(aes(x = topo3, y = richness, fill = topo3)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, width = 0.5) +
  geom_jitter(aes(colour = topo3), width = 0.15, size = 2, alpha = 0.65) +
  geom_text(data = topo3_sum,
            aes(x = topo3, y = ymax + 1.5, label = letter),
            inherit.aes = FALSE, size = 5, fontface = "bold") +
  scale_fill_manual(values = topo3_cols,   labels = topo3_labels) +
  scale_colour_manual(values = topo3_cols, labels = topo3_labels) +
  scale_x_discrete(labels = topo3_labels) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none") +
  labs(title    = "Family richness by environmental contact zone",
       subtitle = sprintf("One-way ANOVA F(%d,%d)=%.2f, p=%.3f | Tukey letters above boxes",
                          summary(m2)[[1]]$Df[1],
                          summary(m2)[[1]]$Df[2],
                          summary(m2)[[1]]$`F value`[1],
                          summary(m2)[[1]]$`Pr(>F)`[1]),
       x = NULL, y = "Families present")

ggsave("anova2_richness_by_topo3.pdf", p2, width = 7, height = 5)
ggsave("anova2_richness_by_topo3.png", p2, width = 7, height = 5, dpi = 300)
cat("Saved anova2_richness_by_topo3.pdf/.png\n")
