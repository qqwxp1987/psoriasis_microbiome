# --- 0. Load Libraries -------------------------------------------------------
library(here)
library(dplyr)
library(ggplot2)

# --- 1. Load Data ------------------------------------------------------------
# The KO annotations for the Fimenecus sp000432435 (SGB4348) genome
# were generated via the eggNOG-mapper pipeline.
# were generated via the eggNOG-mapper pipeline.
# To support reproducible non-interactive execution, we directly load the
# KEGG enrichment results previously generated via MicrobiomeProfiler.
res_dir <- here("CodeAvailiable", "outputs", "genomic_annotation")
res_file <- file.path(res_dir, "Enrichment_Res.csv")
enrich_res <- read.csv(res_file, row.names = 1)

# Filter for significant pathways (pvalue < 0.05)
# This results in the 15 KEGG pathways shown in Figure 4
sig_res <- enrich_res %>%
    filter(pvalue < 0.05) %>%
    arrange(pvalue) %>%
    head(15)

# Format the GeneRatio column from "x/y" to a numeric proportion
sig_res <- sig_res %>%
    mutate(GeneRatio_Num = sapply(strsplit(as.character(GeneRatio), "/"), function(x) {
        as.numeric(x[1]) / as.numeric(x[2])
    }))

# Make Description a factor for ordered plotting
sig_res$Description <- factor(sig_res$Description, levels = rev(sig_res$Description))

# --- 2. Generate Figure 4 (Dot Plot) -----------------------------------------
p_dot <- ggplot(sig_res, aes(x = GeneRatio_Num, y = Description)) +
    geom_point(aes(size = Count, color = pvalue)) +
    scale_color_gradient(low = "red", high = "blue") +
    theme_bw() +
    labs(
        title = "KEGG Enrichment of Fimenecus sp000432435 Genome",
        x = "Gene Ratio",
        y = "",
        color = "pvalue",
        size = "Count"
    ) +
    theme(
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5)
    )

# --- 3. Save Outputs ---------------------------------------------------------
out_dir <- here("CodeAvailiable", "outputs", "genomic_annotation")
if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
}

ggsave(
    filename = file.path(out_dir, "Figure4_SGB4348_KEGG_Enrichment_DotPlot.pdf"),
    plot = p_dot, width = 8, height = 6
)

# Print a summary to the console
cat("\n=== Top 15 Significant KEGG Pathways (pvalue < 0.05) ===\n")
print(sig_res[, c("Description", "Count", "pvalue")])
cat("\nPlot saved to:", file.path(out_dir, "Figure4_SGB4348_KEGG_Enrichment_DotPlot.pdf"), "\n")
