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
  mutate(topology = factor(topology),
         horse    = sub("-.*", "", sample_id))

# -------------------------
# Reduced topology (3 levels based on environmental contact)
# -------------------------
dat <- dat %>%
  mutate(topo3 = case_when(
    topology %in% c("Left front pastern", "Muzzle")        ~ "GCtS_A",
    topology %in% c("Ventral abdomen",    "Udder")          ~ "GCtS_B",
    topology %in% c("Dorsum", "Forehead", "Neck",
                    "Pectoral area")                         ~ "EtS"
  ) %>% factor(levels = c("GCtS_A", "GCtS_B", "EtS")))

cat("=== Topology grouping (reduced) ===\n")
print(dat %>% distinct(topology, topo3) %>% arrange(topo3))

# -------------------------
# Helper: type II marginal F-tests via drop1
# -------------------------
drop1_table <- function(model) {
  d <- drop1(model, test = "F")
  as.data.frame(d) %>%
    rownames_to_column("term") %>%
    filter(term != "<none>") %>%
    rename(df = Df, SS = `Sum of Sq`, F = `F value`, p = `Pr(>F)`) %>%
    mutate(across(c(SS, F), \(x) round(x, 3)),
           p   = round(p, 4),
           sig = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**",
                           p < 0.05  ~ "*",   p < 0.1  ~ ".",
                           TRUE ~ "ns")) %>%
    select(term, df, SS, F, p, sig)
}

# -------------------------
# ANOVA 1: richness ~ topology (8 levels) + age
# -------------------------
m1     <- lm(richness ~ topology + age, data = dat)
a1_tbl <- drop1_table(m1)

cat("\n=== ANOVA 1: richness ~ topology (8 sites) + age | Type II marginal F ===\n")
print(a1_tbl)
write.csv(a1_tbl, "anova1_topology_age.csv", row.names = FALSE)
cat("Saved anova1_topology_age.csv\n")

# -------------------------
# ANOVA 2: richness ~ topo3 (3 levels) + age
# -------------------------
m2     <- lm(richness ~ topo3 + age, data = dat)
a2_tbl <- drop1_table(m2)

cat("\n=== ANOVA 2: richness ~ reduced topology (3 groups) + age | Type II marginal F ===\n")
print(a2_tbl)
write.csv(a2_tbl, "anova2_topo3_age.csv", row.names = FALSE)
cat("Saved anova2_topo3_age.csv\n")

# Post-hoc Tukey for topo3 if significant
if (a2_tbl$p[a2_tbl$term == "topo3"] < 0.05) {
  cat("\n=== Post-hoc: pairwise Tukey for topo3 ===\n")
  ph <- TukeyHSD(aov(richness ~ topo3 + age, data = dat), which = "topo3")
  print(ph)
}

# -------------------------
# Summary table: mean richness per topo3 group
# -------------------------
cat("\n=== Mean richness per reduced topology group ===\n")
sum_tbl <- dat %>%
  group_by(topo3) %>%
  summarise(n         = n(),
            mean_rich = round(mean(richness), 1),
            sd_rich   = round(sd(richness),   1),
            .groups   = "drop")
print(sum_tbl)
