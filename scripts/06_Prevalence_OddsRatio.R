###############################################################################
# 06_Prevalence_OddsRatio.R
# Prevalence analysis for Fimenecus sp000432435 (s__GGB51647_SGB4348)
# - Detection rate comparison across Pso, RC, HC groups
# - Fisher's exact test (overall + pairwise)
# - Odds Ratio, Sensitivity, Specificity, PPV, NPV
# - Prevalence bar chart + waffle plot visualisation
# - Non-zero abundance ANCOVA comparison
#
# Required input files:
#   data/clean/metaphlan/phyloseq_taxa_20240401.rds
#
# Required R packages:
#   here, phyloseq, dplyr, ggplot2, ggpubr, rstatix, patchwork,
#   waffle, microbiomeMarker, emmeans
###############################################################################

# --- 0. Load libraries -------------------------------------------------------
library(here)
library(phyloseq)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(patchwork)
library(waffle)
library(emmeans)

# --- 1. Data loading and preparation -----------------------------------------
phylo.obj <- readRDS(here("data", "clean", "metaphlan", "phyloseq_taxa_20240401.rds"))
data.obj <- phylo.obj

Species <- "k__Bacteria;p__Firmicutes;c__CFGB4806;o__OFGB4806;f__FGB4806;g__GGB51647;s__GGB51647_SGB4348"

# Normalise to relative abundance (TSS)
data.obj <- microbiomeMarker::normalize(data.obj, method = "TSS")
otu <- data.obj@otu_table@.Data %>% as.data.frame()
Abundance <- otu[Species, ] %>% unlist()

meta <- data.obj@sam_data
data <- cbind(Abundance, meta)
data <- data %>%
    mutate(Group = case_when(
        Disease == "Healthy" & Dataset == "Train" ~ "RC",
        Disease == "Psoriasis" ~ "Pso",
        Disease == "Healthy" & Dataset == "Validate" ~ "HC"
    )) %>%
    mutate(Presence = ifelse(Abundance > 0, "Present", "Absent")) %>%
    select(Abundance, Presence, Group, Gender, BMI)

data$Group <- factor(data$Group, levels = c("Pso", "RC", "HC"))
data$Gender <- as.factor(data$Gender)

# Subset: non-zero abundance samples for downstream ANCOVA
data_nonzero <- data %>% filter(Abundance > 0)

# --- 2. Prevalence (detection rate) analysis ----------------------------------

## 2a. Overall Fisher's exact test
contingency_table <- table(data$Group, data$Presence)
cat("Contingency table:\n")
print(contingency_table)

fisher_overall <- fisher.test(contingency_table)
cat("\nOverall Fisher's exact test:\n")
print(fisher_overall)

## 2b. Pairwise Fisher's tests
pairwise_fisher <- pairwise_fisher_test(contingency_table, p.adjust.method = "fdr")
cat("\nPairwise Fisher's tests (FDR-corrected):\n")
print(pairwise_fisher)

# --- 3. Odds Ratio, Sensitivity, Specificity ---------------------------------
# Combine RC and HC as "Control" for binary classification metrics
# Matrix: rows = Present/Absent, cols = Psoriasis/Control
TP <- sum(data$Group == "Pso" & data$Presence == "Present")
FP <- sum(data$Group != "Pso" & data$Presence == "Present")
FN <- sum(data$Group == "Pso" & data$Presence == "Absent")
TN <- sum(data$Group != "Pso" & data$Presence == "Absent")

matrix_data <- matrix(c(TP, FP, FN, TN), nrow = 2, byrow = TRUE)
rownames(matrix_data) <- c("Present", "Absent")
colnames(matrix_data) <- c("Psoriasis", "Control")

result <- fisher.test(matrix_data)
OR_value <- result$estimate
CI_lower <- result$conf.int[1]
CI_upper <- result$conf.int[2]
P_value <- result$p.value

cat(sprintf(
    "\nOR: %.2f (95%% CI: %.2f - %.2f), p = %.4e\n",
    OR_value, CI_lower, CI_upper, P_value
))

Sensitivity <- TP / (TP + FN)
Specificity <- TN / (TN + FP)
PPV <- TP / (TP + FP)
NPV <- TN / (TN + FN)

cat(sprintf(
    "Sensitivity: %.4f\nSpecificity: %.4f\nPPV: %.4f\nNPV: %.4f\n",
    Sensitivity, Specificity, PPV, NPV
))

# --- 4. Prevalence visualisation: bar chart -----------------------------------
prevalence_data <- data %>%
    group_by(Group) %>%
    summarise(
        Total = n(),
        Present_Count = sum(Presence == "Present"),
        Prevalence_Percent = (Present_Count / Total) * 100,
        .groups = "drop"
    )

group_colors <- c("Pso" = "indianred", "RC" = "steelblue", "HC" = "darkslateblue")

barplot_prev <- ggplot(prevalence_data, aes(x = Group, y = Prevalence_Percent, fill = Group)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = paste0(round(Prevalence_Percent, 1), "%")), vjust = -0.5) +
    scale_fill_manual(values = group_colors) +
    scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.1))) +
    labs(title = "(A) Prevalence", x = "", y = "Prevalence (%)") +
    theme_classic() +
    theme(legend.position = "none")

ggsave(barplot_prev,
    filename = here(
        "outputs", "Differential",
        "Bartplot_Prevalence_s__GGB51647_SGB4348.pdf"
    )
)

# --- 5. Waffle chart ----------------------------------------------------------
pso_counts <- contingency_table["Pso", ]
rc_counts <- contingency_table["RC", ]
hc_counts <- contingency_table["HC", ]

pso_waffle <- waffle(rev(pso_counts),
    rows = 10,
    colors = c("Present" = "indianred", "Absent" = "grey90"),
    title = "Pso"
) + theme(legend.position = "none")
rc_waffle <- waffle(rev(rc_counts),
    rows = 10,
    colors = c("Present" = "steelblue", "Absent" = "grey90"),
    title = "RC"
) + theme(legend.position = "none")
hc_waffle <- waffle(rev(hc_counts),
    rows = 10,
    colors = c("Present" = "darkslateblue", "Absent" = "grey90"),
    title = "HC"
) + theme(legend.position = "none")

plot_a_waffle <- (pso_waffle | rc_waffle | hc_waffle) +
    plot_annotation(title = "Prevalence of s__GGB51647_SGB4348")

ggsave(plot_a_waffle,
    filename = here(
        "outputs", "Differential",
        "Waffle_Prevalence_s__GGB51647_SGB4348.pdf"
    )
)

# --- 6. Non-zero abundance analysis (ANCOVA) ----------------------------------
data_nonzero_log <- data_nonzero %>%
    mutate(log10_Abundance = log10(Abundance))

# ANCOVA adjusting for Gender + BMI
ancova_model <- lm(log10_Abundance ~ Group + Gender + BMI, data = data_nonzero_log)
cat("\nOverall ANCOVA result on log-transformed abundance:\n")
print(anova(ancova_model))

# Pairwise comparisons using emmeans
pairwise_ancova <- emmeans(ancova_model, pairwise ~ Group)
cat("\nPairwise comparisons (adjusted for Gender and BMI):\n")
print(pairwise_ancova$contrasts)

