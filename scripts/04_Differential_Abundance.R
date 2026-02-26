# --- 0. Load libraries -------------------------------------------------------
library(here)
library(dplyr)
library(phyloseq)
library(microbiomeMarker)
library(MicrobiomeStat)
library(Maaslin2)
library(mice)
library(ggplot2)

# =============================================================================
# PART 1: Species-level LEfSe on 28 matched pairs (Discovery set)
# =============================================================================

# --- 1.1 Data loading --------------------------------------------------------
phylo.obj <- readRDS(here("data", "metaphlan", "phyloseq_taxa.rds"))
phylo28pairs.obj <- subset_samples(phylo.obj, Dataset == "Train")

# --- 1.2 BMI imputation ------------------------------------------------------
meta.dat <- phylo28pairs.obj@sam_data %>%
    as.matrix() %>%
    as.data.frame() %>%
    select(Disease, Family, Gender, Age, BMI)

meta.dat$Disease <- as.factor(meta.dat$Disease)
meta.dat$Gender <- as.factor(meta.dat$Gender)
meta.dat$Age <- as.numeric(meta.dat$Age)
meta.dat$BMI <- as.numeric(meta.dat$BMI)

imputed_data <- mice(meta.dat %>% select(-Family), m = 5, method = "pmm", maxit = 50, seed = 123)
data_imputed <- complete(imputed_data)
data_imputed <- data_imputed %>% mutate(family = meta.dat$Family)
phylo28pairs.obj@sam_data$BMI <- data_imputed$BMI

ordercolors <- c("steelblue", "indianred")

# --- 1.3 Run LEfSe -----------------------------------------------------------
dir.create(here("outputs", "differential_abundance"), recursive = TRUE, showWarnings = FALSE)
res <- run_lefse(
    ps = phylo28pairs.obj,
    group = "Disease",
    subgroup = NULL,
    taxa_rank = "all",
    transform = "identity",
    norm = "CPM",
    kw_cutoff = 0.05,
    lda_cutoff = 2,
    bootstrap_n = 30,
    bootstrap_fraction = 2 / 3,
    wilcoxon_cutoff = 0.05,
    multigrp_strat = FALSE,
    strict = "0",
    sample_min = 10,
    only_same_subgrp = FALSE,
    curv = FALSE
)

# --- 1.4 LEfSe plots ---------------------------------------------------------
p1 <- plot_ef_dot(res)
ggsave(p1,
    filename = here("outputs", "differential_abundance", "LEfsE_28pairs_ef_dot.pdf"),
    width = 20, height = 10
)

p2 <- plot_abundance(res, group = "Disease")
ggsave(p2,
    filename = here("outputs", "differential_abundance", "LEfsE_28pairs_abundance.pdf"),
    width = 20, height = 10
)

p3 <- plot_cladogram(res, color = ordercolors)
ggsave(p3,
    filename = here("outputs", "differential_abundance", "LEfsE_28pairs_cladogram.pdf"),
    width = 25, height = 20
)

# --- 1.5 Extract markers -----------------------------------------------------
res.marker <- subset_marker(res)
marker.tab <- res.marker@marker_table %>% as.matrix()
write.csv(marker.tab, file = here("outputs", "differential_abundance", "LefsE_28pairs_markers_res.csv"))

marker.profile <- res.marker@otu_table@.Data %>%
    as.data.frame() %>%
    dplyr::filter(rownames(.) %in% marker.tab[, 1])
write.csv(marker.profile, file = here("outputs", "differential_abundance", "LefsE_28pairs_markers_profile.csv"))


# =============================================================================
# PART 2: Species-level LinDA (MicrobiomeStat) on 28 pairs
# =============================================================================

# --- 2.1 Data loading --------------------------------------------------------
phylo28pairs.obj2 <- subset_samples(phylo.obj, Dataset == "Train")

# Convert to MicrobiomeStat object
data.obj2 <- mStat_convert_phyloseq_to_data_obj(phylo28pairs.obj2)

# Use Species-level names as row names for clarity
rownames(data.obj2$feature.tab) <- data.obj2$feature.ann[, 7]
rownames(data.obj2$feature.ann) <- data.obj2$feature.ann[, 7]

# --- 2.2 BMI imputation -------------------------------------------------------
meta.dat2 <- data.obj2$meta.dat %>% select(Disease, Family, Gender, Age, BMI)
meta.dat2$Disease <- as.factor(meta.dat2$Disease)
meta.dat2$Gender <- as.factor(meta.dat2$Gender)
meta.dat2$Age <- as.numeric(meta.dat2$Age)
meta.dat2$BMI <- as.numeric(meta.dat2$BMI)

imputed2 <- mice(meta.dat2 %>% select(-Family), m = 5, method = "pmm", maxit = 50, seed = 123)
data_imputed2 <- complete(imputed2)
data_imputed2 <- data_imputed2 %>% mutate(family = meta.dat2$Family)
data.obj2$meta.dat <- data_imputed2

# --- 2.3 LinDA with Gender + BMI adjustment -----------------------------------
test.list <- generate_taxa_test_single(
    data.obj         = data.obj2,
    group.var        = "Disease",
    adj.vars         = c("Gender", "BMI"),
    prev.filter      = 0.1,
    abund.filter     = 0.0001,
    feature.level    = "Species",
    feature.dat.type = "count"
)

# --- 2.4 Volcano plot ---------------------------------------------------------
dir.create(here("outputs", "differential_abundance", "LinDA"), recursive = TRUE, showWarnings = FALSE)
volcano_plots <- generate_taxa_volcano_single(
    data.obj          = data.obj2,
    group.var         = "Disease",
    test.list         = test.list,
    feature.sig.level = 0.1,
    palette           = c("#F9F871", "#F4A261", "#FF6347"),
    feature.mt.method = "fdr"
)

ggsave(volcano_plots$Species$`Psoriasis vs Healthy (Reference)`,
    filename = here(
        "outputs", "differential_abundance", "LinDA",
        "Adjust_Gender_BMI_28pairs.pdf"
    ),
    width = 10, height = 6
)


# =============================================================================
# PART 3: Species MaAsLin2 validation (28 pairs, adjusting Gender + BMI)
# =============================================================================

# --- 3.1 Prepare input -------------------------------------------------------
input_data3 <- phylo28pairs.obj2@otu_table@.Data %>% as.data.frame()
taxtab3 <- phylo28pairs.obj2@tax_table@.Data %>% as.data.frame()
rownames(input_data3) <- taxtab3$Species
input_data3 <- input_data3 %>% dplyr::filter(rowSums(.) > 0)
input_data3 <- prop.table(as.matrix(input_data3), margin = 2) %>% as.data.frame()

input_meta3 <- phylo28pairs.obj2@sam_data %>%
    as.matrix() %>%
    as.data.frame()
input_meta3$Age <- as.numeric(input_meta3$Age)
input_meta3$BMI <- as.numeric(input_meta3$BMI)
imputed3 <- mice(input_meta3 %>% select(Disease, Gender, Age, BMI),
    m = 5, method = "pmm", maxit = 50, seed = 123
)
data_imputed3 <- complete(imputed3)
data_imputed3 <- data_imputed3 %>%
    mutate(Family = input_meta3$Family) %>%
    arrange(Family)
input_meta3 <- data_imputed3

# --- 3.2 Run MaAsLin2 --------------------------------------------------------
dir.create(here("outputs", "differential_abundance", "MaAsLin2"), recursive = TRUE, showWarnings = FALSE)
fit_data3 <- Maaslin2(
    input_data = input_data3,
    input_metadata = input_meta3,
    output = here("outputs", "differential_abundance", "MaAsLin2"),
    min_abundance = 0.0,
    min_prevalence = 0.1,
    min_variance = 0.0,
    normalization = "NONE",
    transform = "NONE",
    analysis_method = "LM",
    max_significance = 0.1,
    random_effects = NULL,
    fixed_effects = c("Gender", "BMI", "Disease"),
    correction = "BH",
    standardize = TRUE,
    cores = 2,
    plot_heatmap = FALSE,
    plot_scatter = TRUE,
    heatmap_first_n = 50,
    reference = "Disease,Healthy"
)


# =============================================================================
# PART 4: Functional features - KO & Pathway LEfSe
# =============================================================================

# --- 4.1 Load metadata for 28 pairs ------------------------------------------
meta_all <- read.delim(here("data", "metadata.txt"), row.names = 1)
meta_28 <- meta_all %>%
    tibble::rownames_to_column("SampleID") %>%
    dplyr::filter(Dataset == "Train")

# --- 4.2 Pathway LEfSe -------------------------------------------------------
path <- read.delim(here("data", "humann", "filter_path_relab.tsv"), row.names = 1)
# Remove UNMAPPED and UNINTEGRATED rows if present
path <- path[!grepl("^UNMAPPED$|^UNINTEGRATED", rownames(path)), ]

path28pairs <- path[, meta_28$SampleID]

phylo.otu_p <- as.matrix(path28pairs) %>% otu_table(taxa_are_rows = TRUE)
tax.tab_p <- data.frame(
    row.names = rownames(path28pairs),
    Genus = rep("Pathway", nrow(path28pairs)),
    Species = rownames(path28pairs)
)
phylo.tax_p <- as.matrix(tax.tab_p) %>% tax_table()
phylo.phen_p <- meta_28 %>%
    tibble::column_to_rownames("SampleID") %>%
    sample_data()
ps_path <- phyloseq(phylo.otu_p, phylo.tax_p, phylo.phen_p)

res_path <- run_lefse(
    ps = ps_path, group = "Disease", subgroup = NULL, taxa_rank = "Species",
    transform = "identity", norm = "CPM", kw_cutoff = 0.05, lda_cutoff = 2,
    bootstrap_n = 30, bootstrap_fraction = 2 / 3, wilcoxon_cutoff = 0.05,
    multigrp_strat = FALSE, strict = "0", sample_min = 10,
    only_same_subgrp = FALSE, curv = FALSE
)

write.csv(subset_marker(res_path)@marker_table %>% as.matrix(),
    file = here("outputs", "differential_abundance", "LefsE_28pairs_pathway_markers_res.csv")
)

# --- 4.3 KO LEfSe ------------------------------------------------------------
KO <- read.delim(here("data", "humann", "filter_ko_relab.tsv"), row.names = 1)
# Remove UNGROUPED row if present
KO <- KO[!grepl("^UNGROUPED$", rownames(KO)), ]

ko28pairs <- KO[, meta_28$SampleID]

phylo.otu_k <- as.matrix(ko28pairs) %>% otu_table(taxa_are_rows = TRUE)
tax.tab_k <- data.frame(
    row.names = rownames(ko28pairs),
    Genus = rep("KO", nrow(ko28pairs)),
    Species = rownames(ko28pairs)
)
phylo.tax_k <- as.matrix(tax.tab_k) %>% tax_table()
ps_ko <- phyloseq(phylo.otu_k, phylo.tax_k, phylo.phen_p)

res_ko <- run_lefse(
    ps = ps_ko, group = "Disease", subgroup = NULL, taxa_rank = "Species",
    transform = "identity", norm = "CPM", kw_cutoff = 0.05, lda_cutoff = 2,
    bootstrap_n = 30, bootstrap_fraction = 2 / 3, wilcoxon_cutoff = 0.05,
    multigrp_strat = FALSE, strict = "0", sample_min = 10,
    only_same_subgrp = FALSE, curv = FALSE
)

write.csv(subset_marker(res_ko)@marker_table %>% as.matrix(),
    file = here("outputs", "differential_abundance", "LefsE_28pairs_KO_markers_res.csv")
)


# =============================================================================
# PART 5: Figure 4A - LEfSe combined KO + Pathway bar plot
# =============================================================================
# Note: This part requires a manually curated input file that combines
# LEfSe results from KO and Pathway analyses. The file is generated from
# the marker tables produced in PART 4 above.
# If the file does not exist, this section will be skipped.

lefse_fig_path <- here("outputs", "differential_abundance", "Function_lefse_res_input4figure.xlsx")
if (file.exists(lefse_fig_path)) {
    dir.create(here("outputs", "differential_abundance"), recursive = TRUE, showWarnings = FALSE)
    data_fig4a <- rio::import(lefse_fig_path)
    data_fig4a$LDA <- ifelse(data_fig4a$enrich_group == "Healthy",
        data_fig4a$ef_lda * (-1), data_fig4a$ef_lda
    )
    data_fig4a <- data_fig4a[order(data_fig4a$Type, -data_fig4a$LDA), ]
    data_fig4a$Name <- factor(data_fig4a$Name, levels = data_fig4a$Name)

    p_fig4a <- ggplot(data_fig4a) +
        aes(x = Name, y = LDA, fill = enrich_group) +
        geom_col() +
        scale_fill_manual(values = c(Healthy = "steelblue", Psoriasis = "indianred")) +
        coord_flip() +
        theme_minimal()

    ggsave(p_fig4a,
        filename = here("outputs", "differential_abundance", "Lefse_path_KO.pdf"),
        width = 20, height = 10
    )
} else {
    message("Skipping PART 5: Function_lefse_res_input4figure.xlsx not found.")
}


# =============================================================================
# PART 6: Differential KO boxplots (free y-axis)
# =============================================================================

library(patchwork)
library(reshape2)

ko_rel <- read.delim(here("data", "humann", "filter_ko_relab.tsv"), row.names = 1)
ko_rel <- ko_rel[!grepl("^UNGROUPED$", rownames(ko_rel)), ]

samp_incl <- meta_28$SampleID
ko_rel_incl <- ko_rel %>% dplyr::select(all_of(samp_incl))
ko_rel_incl <- proportions(as.matrix(ko_rel_incl), margin = 2) %>% as.data.frame()

diff_ko <- c(
    "K02919", "K02500", "K20265", "K01950", "K03564",
    "K00764", "K03281", "K03781", "K21571", "K07114"
)
pld_ko <- ko_rel_incl[diff_ko, ] %>%
    t() %>%
    as.data.frame()
pld_ko <- pld_ko * 10e4
pld_ko$Group <- rep(c("Control", "Psoriasis"), each = 28)

pld_ko_long <- melt(pld_ko)

KO_plot <- ggplot(data = pld_ko_long, aes(x = Group, y = value, fill = Group)) +
    geom_violin(alpha = 0.35, color = NA) +
    geom_boxplot(width = 0.2, linewidth = 0.5, outlier.shape = NA) +
    geom_jitter(shape = 21, size = 2, color = "white", width = 0.1) +
    scale_fill_manual(values = c("steelblue", "indianred")) +
    scale_x_discrete(labels = c("Healthy", "Psoriasis")) +
    labs(y = "Relative Abundance (10e-4)", x = NULL) +
    theme_classic() +
    theme(
        legend.position = "bottom",
        axis.title = element_text(size = 16, color = "black"),
        axis.text = element_text(size = 13, color = "black")
    )

KO_plot_col2 <- KO_plot + facet_wrap(vars(variable), scales = "free_y", ncol = 2L)
dir.create(here("outputs", "differential_abundance"), recursive = TRUE, showWarnings = FALSE)
ggsave(KO_plot_col2,
    filename = here("outputs", "differential_abundance", "KOgenes_boxplot_col2_freey.pdf"),
    height = 15
)


# =============================================================================
# PART 7: 8 differential species scatter plot
# =============================================================================

phylo28p <- subset_samples(phylo.obj, Dataset == "Train")

input_data7 <- phylo28p@otu_table@.Data %>% as.data.frame()
taxtab7 <- phylo28p@tax_table@.Data %>% as.data.frame()
rownames(input_data7) <- taxtab7$Species
input_data7 <- input_data7 %>% dplyr::filter(rowSums(.) > 0)
input_data7 <- prop.table(as.matrix(input_data7), margin = 2) %>% as.data.frame()

diff_species <- c(
    "s__Prevotella_copri_clade_D", "s__GGB51647_SGB4348",
    "s__Lawsonibacter_asaccharolyticus", "s__Clostridiaceae_bacterium",
    "s__Clostridium_fessum", "s__Enterocloster_bolteae",
    "s__Megamonas_funiformis", "s__Parasutterella_SGB9260"
)
df7 <- input_data7[diff_species, ] %>%
    t() %>%
    as.data.frame()
df7$Group <- phylo28p@sam_data$Disease %>% as.character()

library(tidyr)
df7$ID <- rownames(df7)
long_data <- df7 %>%
    gather(key = "Species", value = "Abundance", -ID, -Group)

summary_data <- long_data %>%
    group_by(Group, Species) %>%
    summarise(
        NonZeroCount = sum(Abundance > 0),
        MeanAbundance = ifelse(NonZeroCount > 0, mean(Abundance[Abundance > 0]), 0),
        .groups = "drop"
    )

summary_data <- summary_data %>%
    mutate(FacetGroup = case_when(
        Species %in% c("s__Prevotella_copri_clade_D", "s__Clostridium_fessum") ~ "Prevalence-driven",
        Species %in% c(
            "s__Lawsonibacter_asaccharolyticus", "s__Clostridiaceae_bacterium",
            "s__Enterocloster_bolteae"
        ) ~ "Abundance-driven",
        Species %in% c(
            "s__GGB51647_SGB4348", "s__Megamonas_funiformis",
            "s__Parasutterella_SGB9260"
        ) ~ "Both"
    ))

dir.create(here("outputs", "differential_abundance"), recursive = TRUE, showWarnings = FALSE)
p_8sp <- ggplot(summary_data, aes(x = NonZeroCount, y = MeanAbundance, color = Species)) +
    geom_point(aes(shape = Group), size = 4) +
    geom_line(aes(group = Species), linetype = "dashed", color = "gray") +
    scale_x_continuous(breaks = seq(0, 28, by = 4), limits = c(0, 28)) +
    labs(x = "Number of Non-Zero Values", y = "Mean Abundance (Non-Zero)") +
    facet_wrap(~FacetGroup, scales = "free_y") +
    theme_minimal() +
    theme(legend.position = "right")

ggsave(p_8sp,
    filename = here("outputs", "differential_abundance", "eight_diff_species.pdf"),
    width = 15
)

message("04_Differential_Abundance.R completed successfully.")
