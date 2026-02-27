###############################################################################
# 00_Data_Preparation.R
# Preprocessing script to format MetaPhlAn outputs and metadata into R objects.
#
# Inputs:
#   - data/metaphlan/metaphlan_taxonomic_profiles.tsv
#   - data/metadata.txt
# Outputs (to data/processed/):
#   - phyloseq_taxa.rds
#   - MicrobiomeData.rds
###############################################################################

# --- 0. Load necessary libraries ----------------------------------------------
library(here)
library(dplyr)
library(tidyr)
library(stringr)
library(phyloseq)

# --- 1. Load raw data ---------------------------------------------------------
metaphlan_file <- here("data", "metaphlan", "metaphlan_taxonomic_profiles.tsv")
metadata_file <- here("data", "metadata.txt")

# Read MetaPhlAn profile (disable check.names to keep full taxa strings)
metaphlan <- read.table(metaphlan_file,
  header = TRUE, check.names = FALSE,
  sep = "\t", row.names = 1
)

# Read metadata
metadata <- read.table(metadata_file,
  header = TRUE, check.names = FALSE,
  sep = "\t", row.names = 1
)

# --- 2. Function: Convert MetaPhlAn to phyloseq object ------------------------
metaphlan_to_phyloseq <- function(all_level_profile, phenotype) {
  taxa <- all_level_profile

  # Replace "|" with ";" to avoid regex issues downstream
  rownames(taxa) <- gsub("\\|", ";", rownames(taxa))

  # Filter for Species-level resolution and drop Strain-level ("t__")
  otu <- taxa %>%
    filter(grepl("s__", rownames(taxa)) & !grepl("t__", rownames(taxa)))

  # Generate taxonomic annotation table
  otu_anno <- data.frame(Taxonomy = rownames(otu), row.names = rownames(otu)) %>%
    separate(
      col = Taxonomy,
      into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
      sep = ";", fill = "right"
    )

  # Convert character columns to factors
  otu_anno <- otu_anno %>% mutate_if(is.character, as.factor)
  rownames(otu_anno) <- rownames(otu)

  # Assemble phyloseq object
  phylo_otu <- otu_table(as.matrix(otu), taxa_are_rows = TRUE)
  phylo_anno <- tax_table(as.matrix(otu_anno))
  phylo_phen <- sample_data(phenotype)

  physeq <- phyloseq(phylo_otu, phylo_anno, phylo_phen)
  return(physeq)
}

# --- 3. Process and save phyloseq_taxa.rds ------------------------------------
phylo_obj <- metaphlan_to_phyloseq(metaphlan, metadata)

# Create output directory mapping if not exists
out_dir <- here("data", "processed")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(phylo_obj, file.path(out_dir, "phyloseq_taxa.rds"))
cat("Saved phyloseq object: data/processed/phyloseq_taxa.rds\n")


# --- 4. Process and save MicrobiomeData.rds (for MicrobiomeStat analyses) -----
# Extract components from phyloseq
feature_tab <- as.data.frame(otu_table(phylo_obj))
meta_dat <- as.data.frame(as.matrix(sample_data(phylo_obj)))
meta_dat$Subject <- rownames(meta_dat)

# Convert appropriate metadata columns to numeric or factors
numeric_vars <- c(
  "Age", "BMI", "BSA", "PASI", "duration",
  "CD4IFN-g", "CD4IL-17"
) # Matches variables broadly requested
# Select strictly the existing numeric columns from metadata
cols_to_numeric <- intersect(numeric_vars, names(meta_dat))
for (col in cols_to_numeric) {
  meta_dat[[col]] <- as.numeric(meta_dat[[col]])
}

factor_vars <- c("Disease", "Gender", "Smoke", "Drink", "Dataset")
cols_to_factor <- intersect(factor_vars, names(meta_dat))
for (col in cols_to_factor) {
  meta_dat[[col]] <- as.factor(meta_dat[[col]])
}

feature_ann <- as.data.frame(tax_table(phylo_obj))

# Replace generic row names with final Species name string to be clean
# Remove "s__" prefix for aesthetics or keep it; we'll keep it exactly as-is in $Species
rownames(feature_ann) <- feature_ann$Species
rownames(feature_tab) <- feature_ann$Species

# Assemble into list structure expected by MicrobiomeStat
MicrobiomeData_metaphlan <- list(
  feature.tab = as.matrix(feature_tab),
  meta.dat = meta_dat,
  feature.ann = as.matrix(feature_ann)
)

saveRDS(MicrobiomeData_metaphlan, file.path(out_dir, "MicrobiomeData.rds"))
cat("Saved sequence list object: data/processed/MicrobiomeData.rds\n")

cat("00_Data_Preparation.R execution completed successfully.\n")
