###############################################################################
# 01_Alpha_Diversity.R
# Alpha diversity analysis for gut metagenome data
# - Shannon and Pielou indices
# - MaAsLin2 for covariate-adjusted comparisons (Gender, BMI; Family as random effect)
# - Visualization: half-violin + half-boxplot + half-point (Figure 1d)
#
# Required input files:
#   data/clean/metaphlan/MicrobiomeData.rds
#
# Required R packages:
#   here, dplyr, MicrobiomeStat, Maaslin2, ggplot2, gghalves, patchwork
###############################################################################

# --- 0. Load libraries -------------------------------------------------------
library(here)
library(dplyr)
library(MicrobiomeStat)
library(Maaslin2)
library(ggplot2)
library(gghalves)
library(patchwork)

# --- 1. Data loading ----------------------------------------------------------
data.obj <- readRDS(here("data", "metaphlan", "MicrobiomeData.rds"))

# --- 2. Calculate alpha diversity indices -------------------------------------
alpha.obj <- mStat_calculate_alpha_diversity(
    data.obj$feature.tab,
    c("shannon", "simpson", "observed_species", "chao1", "ace", "pielou")
)

alpha_tab <- cbind(
    alpha.obj$shannon,
    alpha.obj$simpson,
    alpha.obj$observed_species,
    alpha.obj$pielou
)

# --- 3. Prepare metadata with group labels ------------------------------------
alpha_meta <- data.obj$meta.dat %>%
    mutate(Group = case_when(
        Disease == "Psoriasis" ~ "Psoriasis",
        Disease == "Healthy" & Dataset == "Train" ~ "Relatives",
        Disease == "Healthy" & Dataset == "Validate" ~ "Control"
    )) %>%
    dplyr::select(Gender, Age, BMI, Group, Family)

alpha_meta$Group <- factor(alpha_meta$Group, levels = c("Psoriasis", "Relatives", "Control"))
alpha_meta$Gender <- as.factor(alpha_meta$Gender)

# --- 4. MaAsLin2: Three-group overall comparison -----------------------------
#   Adjusts for Gender and BMI as fixed effects, Family as random effect
Maaslin2(
    input_data = alpha_tab,
    input_metadata = alpha_meta,
    output = here("outputs", "alpha_diversity", "Maaslin", "All_samples"),
    min_abundance = 0.0,
    min_prevalence = 0.1,
    min_variance = 0.0,
    normalization = "none",
    transform = "none",
    analysis_method = "LM",
    max_significance = 0.05,
    random_effects = "Family",
    fixed_effects = c("Gender", "BMI", "Group"),
    correction = "BH",
    standardize = FALSE,
    cores = 8,
    plot_heatmap = FALSE,
    plot_scatter = FALSE,
    heatmap_first_n = 50,
    reference = NULL
)

# --- 5. MaAsLin2: Pairwise comparisons ---------------------------------------

## 5a. Control vs Psoriasis
input_meta_cp <- alpha_meta %>% dplyr::filter(Group != "Relatives")
Maaslin2(
    input_data = alpha_tab,
    input_metadata = input_meta_cp,
    output = here("outputs", "alpha_diversity", "Maaslin", "98Pso_vs_17Con"),
    min_abundance = 0.0,
    min_prevalence = 0.1,
    min_variance = 0.0,
    normalization = "none",
    transform = "none",
    analysis_method = "LM",
    max_significance = 0.05,
    random_effects = NULL,
    fixed_effects = c("Gender", "BMI", "Group"),
    correction = "BH",
    standardize = FALSE,
    cores = 8,
    plot_heatmap = FALSE,
    plot_scatter = FALSE,
    heatmap_first_n = 50,
    reference = NULL
)

## 5b. Relatives vs Psoriasis
input_meta_rp <- alpha_meta %>% dplyr::filter(Group != "Control")
Maaslin2(
    input_data = alpha_tab,
    input_metadata = input_meta_rp,
    output = here("outputs", "alpha_diversity", "Maaslin", "98Pso_vs_28Rel"),
    min_abundance = 0.0,
    min_prevalence = 0.1,
    min_variance = 0.0,
    normalization = "none",
    transform = "none",
    analysis_method = "LM",
    max_significance = 0.05,
    random_effects = "Family",
    fixed_effects = c("Gender", "BMI", "Group"),
    correction = "BH",
    standardize = FALSE,
    cores = 8,
    plot_heatmap = FALSE,
    plot_scatter = FALSE,
    heatmap_first_n = 50,
    reference = NULL
)

## 5c. Control vs Relatives
input_meta_cr <- alpha_meta %>% dplyr::filter(Group != "Psoriasis")
Maaslin2(
    input_data = alpha_tab,
    input_metadata = input_meta_cr,
    output = here("outputs", "alpha_diversity", "Maaslin", "28Rel_vs_17Con"),
    min_abundance = 0.0,
    min_prevalence = 0.1,
    min_variance = 0.0,
    normalization = "none",
    transform = "none",
    analysis_method = "LM",
    max_significance = 0.05,
    random_effects = NULL,
    fixed_effects = c("Gender", "BMI", "Group"),
    correction = "BH",
    standardize = FALSE,
    cores = 8,
    plot_heatmap = FALSE,
    plot_scatter = FALSE,
    heatmap_first_n = 50,
    reference = NULL
)

# --- 6. Visualization: Half-violin + boxplot (Figure 1d) ---------------------
ordercolors <- c("indianred", "darkslateblue", "steelblue")
pld <- cbind(alpha_meta, alpha_tab)

# Shannon index
psh <- ggplot(data = pld, aes(x = Group, y = shannon, fill = Group)) +
    geom_half_violin(side = "r", color = NA, alpha = 0.35) +
    geom_half_boxplot(side = "r", errorbar.draw = FALSE, width = 0.2, linewidth = 0.5) +
    geom_half_point_panel(side = "l", shape = 21, size = 3, color = "white") +
    scale_fill_manual(values = ordercolors) +
    labs(y = "Shannon Index", x = NULL) +
    theme_classic() +
    theme(
        legend.position = "bottom",
        axis.title = element_text(size = 16, color = "black"),
        axis.text = element_text(size = 13, color = "black")
    )

# Pielou evenness index
ppi <- ggplot(data = pld, aes(x = Group, y = pielou, fill = Group)) +
    geom_half_violin(side = "r", color = NA, alpha = 0.35) +
    geom_half_boxplot(side = "r", errorbar.draw = FALSE, width = 0.2, linewidth = 0.5) +
    geom_half_point_panel(side = "l", shape = 21, size = 3, color = "white") +
    scale_fill_manual(values = ordercolors) +
    labs(y = "Pielou Index", x = NULL) +
    theme_classic() +
    theme(
        legend.position = "bottom",
        axis.title = element_text(size = 16, color = "black"),
        axis.text = element_text(size = 13, color = "black")
    )

# Combine and save
ggsave(psh + ppi,
    filename = here("outputs", "alpha_diversity", "Figure1d.pdf"),
    width = 10, height = 8
)
