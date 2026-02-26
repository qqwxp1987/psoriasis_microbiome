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
dir.create(here("outputs", "alpha_diversity", "Maaslin", "All_samples"), recursive = TRUE, showWarnings = FALSE)
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
dir.create(here("outputs", "alpha_diversity", "Maaslin", "98Pso_vs_17Con"), recursive = TRUE, showWarnings = FALSE)
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
dir.create(here("outputs", "alpha_diversity", "Maaslin", "98Pso_vs_28Rel"), recursive = TRUE, showWarnings = FALSE)
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
dir.create(here("outputs", "alpha_diversity", "Maaslin", "28Rel_vs_17Con"), recursive = TRUE, showWarnings = FALSE)
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

# --- 6. Visualization: violin + boxplot + jitter (Figure 1d) ------------------
ordercolors <- c("indianred", "steelblue", "darkslateblue")
pld <- cbind(alpha_meta, alpha_tab)

# Shannon index
psh <- ggplot(data = pld, aes(x = Group, y = shannon, fill = Group)) +
    geom_violin(alpha = 0.35, color = NA) +
    geom_boxplot(width = 0.2, linewidth = 0.5, outlier.shape = NA) +
    geom_jitter(shape = 21, size = 2, color = "white", width = 0.1) +
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
    geom_violin(alpha = 0.35, color = NA) +
    geom_boxplot(width = 0.2, linewidth = 0.5, outlier.shape = NA) +
    geom_jitter(shape = 21, size = 2, color = "white", width = 0.1) +
    scale_fill_manual(values = ordercolors) +
    labs(y = "Pielou Index", x = NULL) +
    theme_classic() +
    theme(
        legend.position = "bottom",
        axis.title = element_text(size = 16, color = "black"),
        axis.text = element_text(size = 13, color = "black")
    )

# Combine and save
dir.create(here("outputs", "alpha_diversity"), recursive = TRUE, showWarnings = FALSE)
ggsave(psh + ppi,
    filename = here("outputs", "alpha_diversity", "alpha_diversity_shannon_pielou.pdf"),
    width = 10, height = 8
)
