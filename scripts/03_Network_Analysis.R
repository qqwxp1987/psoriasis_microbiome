###############################################################################
# 03_Network_Analysis.R
# Co-occurrence network analysis using ggClusterNet
# - Three-group network construction (Psoriasis, Family_control, Healthy)
# - Network property comparison (ZiPi, degree, etc.)
# - Network stability: module comparison, random removal, targeted removal
#
# Required input files:
#   data/clean/metaphlan/phyloseq_taxa_20240401.rds
#   scripts/ggClusterNetFig.R  (custom modified network.i function)
#
# Required R packages:
#   here, dplyr, phyloseq, ggClusterNet, igraph, sna, ggplot2,
#   microbiome, ggnewscale, tibble, tidyverse, tidyfst
###############################################################################

# --- 0. Load libraries -------------------------------------------------------
library(here)
library(dplyr)
library(igraph)
library(sna)
library(phyloseq)
library(ggClusterNet)
library(ggplot2)
library(microbiome)
library(ggnewscale)
library(tibble)

# --- 1. Data loading ----------------------------------------------------------
ps0 <- readRDS(here("data", "metaphlan", "phyloseq_taxa_20240401.rds"))

# Assign group labels: Family_control for matched healthy relatives
sample_data(ps0)$Group <- case_when(
    sample_data(ps0)$Disease == "Healthy" & sample_data(ps0)$Dataset == "Train" ~ "Family_control",
    .default = as.character(sample_data(ps0)$Disease)
)
ps <- ps0

# Collapse rare phyla into "Others"
tax_table_df <- as.data.frame(tax_table(ps))
common_phyla <- c("p__Actinobacteria", "p__Bacteroidetes", "p__Firmicutes", "p__Proteobacteria")
tax_table_df <- tax_table_df %>%
    mutate(Phylum = ifelse(Phylum %in% common_phyla, Phylum, "Others"))
tax_table(ps) <- as.matrix(tax_table_df)

# --- 2. Custom theme ----------------------------------------------------------
mytheme1 <- theme_classic() +
    theme(
        panel.background   = element_blank(),
        panel.grid         = element_blank(),
        legend.position    = "right",
        legend.title       = element_blank(),
        legend.background  = element_blank(),
        legend.key         = element_blank(),
        plot.title         = element_text(vjust = -8.5, hjust = 0.1),
        axis.title.y       = element_text(size = 15, face = "bold", colour = "black"),
        axis.title.x       = element_text(size = 15, face = "bold", colour = "black"),
        axis.text          = element_text(size = 10, face = "bold"),
        axis.text.x        = element_text(colour = "black", size = 10),
        axis.text.y        = element_text(colour = "black", size = 10),
        legend.text        = element_text(size = 10, face = "bold")
    )

# --- 3. Load custom network function -----------------------------------------
# This file contains a modified `network.i.mod()` function
# with customised color scheme
source(here("scripts", "ggClusterNetFig.R"))

# --- 4. Data transformation --------------------------------------------------
ps.log10p <- microbiome::transform(ps, "log10p")
netpath <- here("outputs", "network_igraph", "20240623")
dir.create(netpath, recursive = TRUE, showWarnings = FALSE)
gnum <- 3

# --- 5. Network construction: 143 samples, 3 groups --------------------------
result <- network.i.mod(
    ps              = ps.log10p,
    N               = 0,
    r.threshold     = 0.8,
    big             = TRUE,
    select_layout   = FALSE,
    method          = "pearson",
    scale           = FALSE,
    layout_net      = "model_igraph2",
    p.threshold     = 0.05,
    label           = FALSE,
    group           = "Group",
    fill            = "Phylum",
    size            = "igraph.degree",
    path            = netpath,
    maxnode         = 5,
    ncol            = 3,
    nrow            = 1,
    zipi            = TRUE,
    order           = NULL,
    ncpus           = 8
)

# Save plots and tables
p1 <- result[[1]]
dat <- result[[2]]
p <- result[[5]]

write.csv(dat, paste0(netpath, "/co-occurrence_Globel_net.csv"))
ggsave(paste0(netpath, "/network_all.pdf"), p1, width = 3 * gnum, height = 2.5, limitsize = FALSE)
ggsave(paste0(netpath, "/network_all2.pdf"), p, width = 6 * gnum, height = 6, limitsize = FALSE)

# --- 6. Network construction: 28 pairs only ----------------------------------
ps28pairs.log10p <- subset_samples(ps.log10p, Dataset == "Train")

result_28 <- network.i.mod(
    ps              = ps28pairs.log10p,
    N               = 0,
    r.threshold     = 0.8,
    big             = TRUE,
    select_layout   = TRUE,
    method          = "pearson",
    scale           = FALSE,
    layout_net      = "model_igraph2",
    p.threshold     = 0.05,
    label           = FALSE,
    path            = netpath,
    maxnode         = 5,
    ncol            = 2,
    nrow            = 1,
    zipi            = TRUE,
    order           = NULL
)

p1_28 <- result_28[[1]]
dat_28 <- result_28[[2]]
p_28 <- result_28[[5]]

write.csv(dat_28, paste0(netpath, "/co-occurrence_28pairs_net.csv"))
ggsave(paste0(netpath, "/network_28pairs.pdf"), p1_28, width = 6, height = 2.5, limitsize = FALSE)
ggsave(paste0(netpath, "/network_28pairs2.pdf"), p_28, width = 12, height = 6, limitsize = FALSE)

# --- 7. Network stability analysis -------------------------------------------
library(tidyverse)
library(tidyfst)

stabpath <- here("outputs", "network_stab", "20240605")
dir.create(stabpath, recursive = TRUE, showWarnings = FALSE)

## 7a. Module comparison
res_mod <- module.compare.m(
    ps          = ps.log10p,
    Top         = 200,
    degree      = TRUE,
    zipi        = TRUE,
    r.threshold = 0.8,
    p.threshold = 0.05,
    method      = "pearson",
    padj        = FALSE,
    n           = 3
)

p_mod <- res_mod[[1]]
dat_mod <- res_mod[[2]]
dat2 <- res_mod[[3]]

dat2$m1 <- dat2$module1 %>%
    strsplit("model") %>%
    sapply(`[`, 1)
dat2$m2 <- dat2$module2 %>%
    strsplit("model") %>%
    sapply(`[`, 1)
dat2$cross <- paste(dat2$m1, dat2$m2, sep = "_Vs_")

p_sim <- ggplot(dat2) +
    geom_bar(aes(x = cross, fill = cross)) +
    labs(x = "", y = "numbers.of.similar.modules") +
    theme_classic()

ggsave(paste0(stabpath, "/module.compare.groups.pdf"), p_mod, width = 10, height = 10)
ggsave(paste0(stabpath, "/numbers.of.similar.modules.pdf"), p_sim, width = 8, height = 8)
write.csv(dat_mod, paste0(stabpath, "/module.otu.csv"), quote = FALSE)
write.csv(dat2, paste0(stabpath, "/module.compare.groups.csv"), quote = FALSE)

## 7b. Random node removal (Robustness)
res_rr <- Robustness.Random.removal(
    ps          = ps,
    Top         = 500,
    r.threshold = 0.8,
    p.threshold = 0.05,
    method      = "pearson"
)
p_rr <- res_rr[[1]]
dat_rr <- res_rr[[2]]

rr_path <- paste0(stabpath, "/Robustness_Random_removal/")
dir.create(rr_path, recursive = TRUE, showWarnings = FALSE)
write.csv(dat_rr, paste0(rr_path, "random_removal_network.csv"))
ggsave(paste0(rr_path, "random_removal_network.pdf"), p_rr, width = 8, height = 4)

## 7c. Targeted node removal (Robustness)
res_tr <- Robustness.Targeted.removal(
    ps          = ps,
    Top         = 500,
    degree      = TRUE,
    zipi        = FALSE,
    r.threshold = 0.8,
    p.threshold = 0.05,
    method      = "pearson"
)
p_tr <- res_tr[[1]]
dat_tr <- res_tr[[2]]

tr_path <- paste0(stabpath, "/Robustness_Targeted_removal/")
dir.create(tr_path, recursive = TRUE, showWarnings = FALSE)
write.csv(dat_tr, paste0(tr_path, "Robustness_Targeted_removal_network.csv"))
ggsave(paste0(tr_path, "Robustness_Targeted_removal_network.pdf"), p_tr, width = 8, height = 4)

## 7d. Network vulnerability
res_vul <- Vulnerability.micro(
    ps          = ps,
    Top         = 500,
    degree      = TRUE,
    zipi        = FALSE,
    r.threshold = 0.8,
    p.threshold = 0.05,
    method      = "spearman"
)
p_vul <- res_vul[[1]] + theme_bw()

message("05_Network_Analysis.R completed successfully.")
