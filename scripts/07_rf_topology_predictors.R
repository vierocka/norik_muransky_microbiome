library(tidyverse)
library(randomForest)

set.seed(42)
setwd("/home/veve/Dropbox/kone/qiime2")

contrib <- read.csv("pca_family_contributions.csv", check.names = FALSE)
meta    <- read.csv("/home/veve/Dropbox/kone/kraken2_reports/metadata.csv",
                    check.names = FALSE) %>%
           mutate(across(everything(), trimws))

fam_cols <- setdiff(colnames(contrib), c("sample_id", "topology"))

# Shorten long taxonomy label
shorten_fam <- function(x) sub(
  "d__Bacteria;p__Patescibacteria;c__Saccharimonadia;o__Saccharimonadales;__",
  "Saccharimonadales", x)

# -------------------------
# Build matrix: animal samples only, A/B/C groups
# -------------------------
dat <- contrib %>%
  filter(topology != "Environment") %>%
  mutate(topo3 = case_when(
    topology %in% c("Left front pastern", "Muzzle")     ~ "GCtS_A",
    topology %in% c("Ventral abdomen",    "Udder")       ~ "GCtS_B",
    topology %in% c("Dorsum", "Forehead", "Neck",
                    "Pectoral area")                      ~ "EtS"
  ) %>% factor(levels = c("GCtS_A", "GCtS_B", "EtS")))

X <- dat %>% select(all_of(fam_cols)) %>%
     mutate(across(everything(), as.numeric)) %>%
     rename_with(shorten_fam)
y <- dat$topo3

# -------------------------
# Identify env-shared vs animal-only families
# -------------------------
env_row <- contrib %>% filter(topology == "Environment") %>%
           select(all_of(fam_cols)) %>%
           mutate(across(everything(), as.numeric))
env_shared  <- shorten_fam(fam_cols[as.numeric(env_row[1,]) >  0])
anim_only_f <- shorten_fam(fam_cols[as.numeric(env_row[1,]) == 0])

# -------------------------
# Random Forest: 1000 trees, importance = TRUE
# -------------------------
rf <- randomForest(x = X, y = y, ntree = 1000, importance = TRUE)

cat("=== Random Forest: OOB confusion matrix ===\n")
print(rf$confusion)
cat(sprintf("\nOOB error rate: %.1f%%\n", rf$err.rate[1000, "OOB"] * 100))

# Variable importance (Mean Decrease Accuracy — most reliable)
imp <- importance(rf, type = 1) %>%
  as.data.frame() %>%
  rownames_to_column("family") %>%
  rename(MDA = MeanDecreaseAccuracy) %>%
  mutate(source = ifelse(family %in% env_shared, "env-shared", "animal-only")) %>%
  arrange(desc(MDA))

cat("\n=== Variable importance (Mean Decrease Accuracy, ranked) ===\n")
print(as.data.frame(imp), row.names = FALSE)
write.csv(imp, "rf_family_importance.csv", row.names = FALSE)
cat("Saved rf_family_importance.csv\n")

cat(sprintf("\nTop 10 predictors — animal-only: %d / 10\n",
            sum(imp$source[1:10] == "animal-only")))
cat(sprintf("Top 20 predictors — animal-only: %d / 20\n",
            sum(imp$source[1:20] == "animal-only")))
cat(sprintf("All 46 families  — animal-only: %d / 46\n", length(anim_only_f)))

# -------------------------
# Plot: variable importance coloured by env origin
# -------------------------
top_n <- 20
imp_plot <- imp %>%
  slice_head(n = top_n) %>%
  mutate(family = factor(family, levels = rev(family)))

src_cols <- c("animal-only" = "#d62728", "env-shared" = "#7f7f7f")

p_imp <- ggplot(imp_plot, aes(x = MDA, y = family, fill = source)) +
  geom_col(width = 0.7, alpha = 0.85) +
  scale_fill_manual(values = src_cols,
                    labels = c("animal-only"  = "Absent from environment",
                               "env-shared"   = "Present in environment"),
                    name   = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position  = c(0.7, 0.15),
        legend.text      = element_text(size = 9),
        plot.title       = element_text(size = 11, face = "bold"),
        plot.subtitle    = element_text(size = 8)) +
  labs(title    = sprintf("Top %d family predictors of body-site contact zone (A/B/C)", top_n),
       subtitle = sprintf("Random Forest (1000 trees) | OOB accuracy = %.1f%% | Mean Decrease Accuracy",
                          (1 - rf$err.rate[1000, "OOB"]) * 100),
       x = "Mean Decrease Accuracy", y = NULL)

ggsave("rf_family_importance.pdf", p_imp, width = 8, height = 6)
ggsave("rf_family_importance.png", p_imp, width = 8, height = 6, dpi = 300)
cat("Saved rf_family_importance.pdf/.png\n")

# -------------------------
# Per-class MDA: which families drive each group
# -------------------------
imp_class <- importance(rf, scale = FALSE) %>%
  as.data.frame() %>%
  rownames_to_column("family") %>%
  mutate(source = ifelse(family %in% env_shared, "env-shared", "animal-only")) %>%
  arrange(desc(GCtS_A))

cat("\n=== Per-class MDA (top 10 per group) ===\n")
for (grp in c("GCtS_A", "GCtS_B", "EtS")) {
  cat(sprintf("\n%s:\n", grp))
  print(imp_class %>% arrange(desc(.data[[grp]])) %>%
        select(family, all_of(grp), source) %>% slice_head(n = 10))
}
write.csv(imp_class, "rf_family_importance_perclass.csv", row.names = FALSE)
cat("\nSaved rf_family_importance_perclass.csv\n")
