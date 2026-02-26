###############################################################################
# 07_Genomic_Annotation.R
# Functional annotation and enrichment analysis
# - MicrobiomeProfiler: KEGG module/pathway enrichment from differential KOs
# - EggNOG-mapper annotation barplot for t__SGB4348 genome
#
# Required input files:
#   Differential KO list from LEfSe results
#
# Required R packages:
#   here, ggplot2, MicrobiomeProfiler
###############################################################################

# --- 0. Load libraries -------------------------------------------------------
library(here)
library(ggplot2)

# --- 1. KEGG enrichment analysis with MicrobiomeProfiler ---------------------
# The differential KOs from LEfSe (script 06) or MaAsLin2 are used as input.
# MicrobiomeProfiler provides a Shiny-based interactive interface.
#
# Differential KOs enriched in Psoriasis:
#   K02919, K01784, K03569, K00773, K02108
#
# Differential KOs enriched in Healthy:
#   K01190, K05349, K02500, K07114, K01187, K03497, K01812, K21571,
#   K01785, K06142, K06871, K01950, K04762, K20265, K00991, K01689,
#   K09807, K03564, K03781, K01681, K03281, K01192, K00764, K04564, K01923, K08961

library(MicrobiomeProfiler)
# Launch the interactive interface to perform enrichment
# run_MicrobiomeProfiler()

# Results summary (from interactive analysis):
# Enriched KEGG pathways (p.adjust < 0.05):
#   map00052  Galactose metabolism          4/19  p=0.006
#   map04212  Longevity regulating pathway  2/19  p=0.020
#   map00511  Other glycan degradation      2/19  p=0.020
#   map04146  Peroxisome                    2/19  p=0.020

cat("
=================================================================
KEGG Enrichment Results (MicrobiomeProfiler)
=================================================================
Pathway                              GeneRatio  BgRatio    p.adjust
map00052 Galactose metabolism        4/19       48/3034    0.006
map04212 Longevity regulating (worm) 2/19       10/3034    0.020
map00511 Other glycan degradation    2/19       12/3034    0.020
map04146 Peroxisome                  2/19       12/3034    0.020
=================================================================
")


# --- 2. EggNOG-mapper annotation summary for t__SGB4348 ----------------------
# Annotation pipeline: Prodigal -> eggNOG-mapper -> KEGG modules
# Data represents gene annotation steps for t__GGB51647_SGB4348

steps <- c(
    "Predicted Genes", "Genes with COGs",
    "Genes with KO Annotations", "KEGG Modules"
)
gene_counts <- c(187, 179, 114, 33)

df <- data.frame(Step = factor(steps, levels = steps), GeneCount = gene_counts)

bar_plot <- ggplot(df) +
    aes(x = Step, y = GeneCount, fill = Step) +
    geom_col() +
    scale_fill_brewer(palette = "YlGnBu", direction = -1) +
    theme_minimal() +
    geom_text(aes(label = GeneCount), vjust = -0.5, size = 5) +
    labs(
        title = "Gene Annotation Process for t__GGB51647_SGB4348",
        y = "Number"
    )

out_dir <- here("outputs", "t__SGB4348_annotation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(bar_plot,
    filename = file.path(out_dir, "GenesNo_stat.pdf"),
    width = 10
)

