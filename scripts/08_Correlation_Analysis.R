###############################################################################
# 08_Correlation_Analysis.R
# Partial Spearman correlation analysis (adjusting for BMI)
# - SGB4348 vs clinical/immune markers (Figure 6b)
# - Species vs PASI (Supplementary)
#
# Required input files:
#   data/metaphlan/phyloseq_taxa.rds
#   data/metadata.txt  (with immune/clinical data)
#
# Required R packages:
#   here, phyloseq, ppcor, dplyr, microbiomeMarker
###############################################################################

# --- 0. Load libraries -------------------------------------------------------
library(here)
library(phyloseq)
library(ppcor)
library(dplyr)

# --- 1. Data loading and normalisation ----------------------------------------
phylo.obj <- readRDS(here("data", "metaphlan", "phyloseq_taxa.rds"))

# Subset to psoriasis patients only
pso.obj <- subset_samples(phylo.obj, Disease == "Psoriasis")

# Filter: remove taxa absent in all samples
pso_filt.obj <- prune_taxa(taxa_sums(pso.obj) > 0, pso.obj)
# Filter: remove taxa present in fewer than half of samples
pso_filt.obj <- filter_taxa(pso_filt.obj, function(x) sum(x > 0) > (0.5 * length(x)), TRUE)

# CLR normalisation (add pseudocount 1 first to handle zeros)
pso_filt.obj <- transform_sample_counts(pso_filt.obj, function(x) x + 1)
data.obj <- microbiomeMarker::normalize(pso_filt.obj, method = "CLR")


# =============================================================================
# PART 1: SGB4348 vs clinical/immune markers (Figure 6b)
# Partial Spearman correlation adjusting for BMI
# =============================================================================

method <- "spearman"

# Extract species abundance
otu_mat <- otu_table(data.obj)
otu_df <- as.data.frame(otu_mat) %>%
    t() %>%
    as.data.frame()

# Read extended metadata with immune/clinical data
sample_data_ext <- read.table(here("data", "metadata.txt"),
    sep = "\t", header = TRUE, row.names = 1
)
sample_data_ext <- sample_data_ext %>% filter(Disease == "Psoriasis")

# Ensure numeric types
sample_data_ext$BMI <- as.numeric(sample_data_ext$BMI)
sample_data_ext$PASI <- as.numeric(sample_data_ext$PASI)
sample_data_ext$IL17Foxp3ratio <- sample_data_ext$IL.17 / sample_data_ext$FOXP3

# --- 1.1 Compute partial correlations ----------------------------------------
# Target species: Fimenecus sp000432435 (SGB4348)
target_species <- "k__Bacteria;p__Firmicutes;c__CFGB4806;o__OFGB4806;f__FGB4806;g__GGB51647;s__GGB51647_SGB4348"
phenotypes <- c(
    "Lymphocytes", "CD4", "CD8", "FOXP3", "IFNg",
    "TNFalpha", "IL.17", "BSA", "PASI", "VAS", "IL17Foxp3ratio"
)

results <- data.frame(
    Species = character(), Index = character(),
    Correlation = numeric(), P_value = numeric(),
    stringsAsFactors = FALSE
)

for (index in phenotypes) {
    species_data <- otu_df[, target_species, drop = FALSE]
    combined_data <- cbind(species_data, sample_data_ext[, c("BMI", index)])
    combined_data <- na.omit(combined_data)

    corr_test <- pcor.test(combined_data[, 1], combined_data[, 3],
        combined_data[, 2],
        method = method
    )

    results <- rbind(results, data.frame(
        Species     = "SGB4348",
        Index       = index,
        Correlation = corr_test$estimate,
        P_value     = corr_test$p.value
    ))
}

# --- 1.2 Save results --------------------------------------------------------
out_dir <- here("outputs", "partial_correlation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(results,
    file.path(out_dir, "SGB4348_Phenotype_spearman_adjustBMI_all.csv"),
    row.names = FALSE
)

cat("\n=== SGB4348 vs Clinical/Immune Markers (Partial Spearman, adjusting BMI) ===\n")
print(results)


# =============================================================================
# PART 2: All species vs PASI (adjusting for BMI)
# =============================================================================

# Use the same CLR-normalised data
sample_meta <- sample_data(data.obj)
sample_meta$BMI <- as.numeric(sample_meta$BMI)
sample_meta$PASI <- as.numeric(sample_meta$PASI)

results_pasi <- data.frame(
    Species = character(), Correlation = numeric(),
    P_value = numeric(), stringsAsFactors = FALSE
)

for (species in colnames(otu_df)) {
    species_data <- otu_df[, species, drop = FALSE]
    combined_data <- cbind(species_data, sample_meta[, c("BMI", "PASI")])
    combined_data <- na.omit(combined_data)

    corr_test <- pcor.test(combined_data[, 1], combined_data[, 3],
        combined_data[, 2],
        method = method
    )

    if (corr_test$p.value < 0.05) {
        results_pasi <- rbind(results_pasi, data.frame(
            Species     = species,
            Correlation = corr_test$estimate,
            P_value     = corr_test$p.value
        ))
    }
}

write.csv(results_pasi,
    file.path(out_dir, "Species_PASI_pSpearman_all.csv"),
    row.names = FALSE
)

cat("\n=== Species significantly correlated with PASI (p < 0.05) ===\n")
print(results_pasi)
