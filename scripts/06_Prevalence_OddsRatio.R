# --- 0. Load libraries -------------------------------------------------------
library(here)
library(phyloseq)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(patchwork)
library(emmeans)

# --- 1. Data loading and preparation -----------------------------------------
phylo.obj <- readRDS(here("data", "metaphlan", "phyloseq_taxa.rds"))
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

dir.create(here("outputs", "differential_abundance"), recursive = TRUE, showWarnings = FALSE)
ggsave(barplot_prev,
    filename = here(
        "outputs", "differential_abundance",
        "Bartplot_Prevalence_s__GGB51647_SGB4348.pdf"
    )
)

# --- 5. Custom ggplot2 Waffle chart (no waffle package dependency) --------------
create_waffle <- function(counts, colors, title) {
    # counts is a named vector, e.g. c("Present" = X, "Absent" = Y)
    # Total should be roughly 100 (since % prevalence)
    total <- sum(counts)
    if (total == 0) {
        return(ggplot() +
            theme_void())
    }

    # Create the sequence of categories
    # The original waffle plot rev() the counts.
    categories <- rep(names(counts), counts)

    # Pad or truncate to 100 to fit a 10x10 matrix (assuming round(counts) aligns)
    # (Here we expect total == nrow of group)
    # We will compute percentage and round to make 100 squares
    percents <- round((counts / total) * 100)

    # Adjust for rounding errors to ensure exactly 100 tiles
    diff <- 100 - sum(percents)
    if (diff != 0) percents[1] <- percents[1] + diff

    cat_seq <- factor(rep(names(percents), percents), levels = names(colors))

    # 10x10 grid (x and y coords)
    df <- data.frame(
        x = rep(1:10, each = 10),
        y = rep(10:1, times = 10),
        category = cat_seq
    )

    ggplot(df, aes(x = x, y = y, fill = category)) +
        geom_tile(color = "white", linewidth = 1) +
        coord_equal() +
        scale_fill_manual(values = colors) +
        labs(title = title) +
        theme_void() +
        theme(
            legend.position = "none",
            plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
        )
}

pso_counts <- contingency_table["Pso", ]
rc_counts <- contingency_table["RC", ]
hc_counts <- contingency_table["HC", ]

pso_waffle <- create_waffle(rev(pso_counts),
    colors = c("Present" = "indianred", "Absent" = "grey90"),
    title = "Pso"
)
rc_waffle <- create_waffle(rev(rc_counts),
    colors = c("Present" = "steelblue", "Absent" = "grey90"),
    title = "RC"
)
hc_waffle <- create_waffle(rev(hc_counts),
    colors = c("Present" = "darkslateblue", "Absent" = "grey90"),
    title = "HC"
)

plot_a_waffle <- (pso_waffle | rc_waffle | hc_waffle) +
    plot_annotation(title = "Prevalence of s__GGB51647_SGB4348")

ggsave(plot_a_waffle,
    filename = here(
        "outputs", "differential_abundance",
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
