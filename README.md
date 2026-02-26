# Psoriasis Gut Metagenome Analysis

Analysis scripts and data for the manuscript: *"Disentangling Environmental and Disease-Specific Signatures in the Gut Microbiome of Psoriasis: Discovery of Fimenecus sp. as a Novel Biomarker and Characterization of the Gut Virome"*

## Directory Structure

```
psoriasis_microbiome/
├── README.md
├── scripts/
│   ├── 01_Alpha_Diversity.R          # Shannon & Pielou indices, MaAsLin2
│   ├── 02_Beta_Diversity.R           # PERMANOVA, JSD, Bray-Curtis distances
│   ├── 03_Network_Analysis.R         # Co-occurrence networks (ggClusterNet)
│   ├── 04_Differential_Abundance.R   # LEfSe, LinDA, MaAsLin2 (species/KO/Pathway)
│   ├── 05_RandomForest_Classifier.R  # Boruta + RF, ROC/PRC, confusion matrix
│   ├── 06_Prevalence_OddsRatio.R     # Prevalence, Fisher test, OR, waffle chart
│   ├── 07_Genomic_Annotation.R       # KEGG enrichment, EggNOG-mapper barplot
│   └── 08_Correlation_Analysis.R     # Partial Spearman correlation (ppcor)
└── data/
    ├── metadata.txt
    ├── humann/
    │   ├── filter_ko_relab.tsv
    │   └── filter_path_relab.tsv
    └── metaphlan/
        ├── MicrobiomeData.rds
        └── phyloseq_taxa.rds
```

## Prerequisites

### R (≥ 4.2)

```r
# CRAN packages
install.packages(c(
  "here", "dplyr", "ggplot2", "patchwork", "gghalves", "vegan",
  "mice", "reshape2", "ppcor", "caret", "Boruta", "pROC", "ROCR",
  "ggpubr", "rstatix", "waffle", "emmeans", "rio", "tidyr",
  "scales", "ggsci", "igraph", "sna", "tidyfst", "ggnewscale",
  "ggstance", "stringr"
))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c(
  "phyloseq", "microbiomeMarker", "MicrobiomeStat", "Maaslin2",
  "microbiome", "MicrobiomeProfiler"
))

# GitHub packages
devtools::install_github("taowenmicro/ggClusterNet")
```

## Script–Figure Mapping

| Script | Analysis | Manuscript Figure |
|--------|----------|-------------------|
| `01_Alpha_Diversity.R` | Alpha diversity indices & statistical testing | Figure 1d |
| `02_Beta_Diversity.R` | Beta diversity (PERMANOVA, PCoA) | Figure 1a, 1b |
| `03_Network_Analysis.R` | Co-occurrence network construction & stability | Figure 3 |
| `04_Differential_Abundance.R` | Differential abundance (species, KO, pathway) | Figure 2, 4 |
| `05_RandomForest_Classifier.R` | Random Forest classification & validation | Figure 6 |
| `06_Prevalence_OddsRatio.R` | Prevalence analysis & diagnostic metrics | Figure 5 |
| `07_Genomic_Annotation.R` | Functional annotation enrichment | Figure 5a |
| `08_Correlation_Analysis.R` | Partial correlation with clinical phenotypes | Figure 7 |

## Notes

- BMI missing values are imputed using the `mice` package (predictive mean matching, seed = 123).
- Family structure is accounted for as a random effect where applicable.
- LEfSe analysis uses the `microbiomeMarker` package implementation.
- Gender and BMI are adjusted as fixed-effect confounders in all differential analyses.
