# --- 0. Load libraries -------------------------------------------------------
library(here)
library(dplyr)
library(vegan)
library(mice)
library(ape)
library(reshape2)
library(phyloseq)
library(ggplot2)
# library(gghalves)  # removed due to ggplot2 v4.5+ incompatibility
library(patchwork)

# --- 1. Data loading ----------------------------------------------------------
data.obj <- readRDS(here("data", "metaphlan", "MicrobiomeData.rds"))

# --- 2. Metadata preparation with BMI imputation -----------------------------
meta <- data.obj$meta.dat %>%
    mutate(Group = case_when(
        Disease == "Psoriasis" ~ "Psoriasis",
        Disease == "Healthy" & Dataset == "Train" ~ "Relatives",
        Disease == "Healthy" & Dataset == "Validate" ~ "Control"
    )) %>%
    dplyr::select(Gender, Age, BMI, Group, Family)

meta$Group <- factor(meta$Group, levels = c("Psoriasis", "Relatives", "Control"))
meta$Gender <- as.factor(meta$Gender)

# Impute missing BMI values using mice
imputed_data <- mice(meta %>% select(-Family), m = 5, method = "pmm", maxit = 50, seed = 123)
data_imputed <- complete(imputed_data)
data_imputed <- data_imputed %>% mutate(Family = meta$Family)

# --- 3. PERMANOVA analysis (adonis2) ------------------------------------------
feature_tab <- data.obj$feature.tab
profile <- prop.table(feature_tab, 2) * 100
profile_t <- t(profile) %>% as.data.frame()

adonis_result <- adonis2(
    profile_t ~ Gender + BMI + Group,
    data_imputed,
    permutations = 999,
    by = "margin",
    distance = "bray"
)
dir.create(here("outputs", "beta_diversity"), recursive = TRUE, showWarnings = FALSE)
write.csv(adonis_result,
    file = here("outputs", "beta_diversity", "PERMANOVA_res.csv"),
    quote = FALSE, na = ""
)
message("PERMANOVA results saved.")

# --- 4. Pairwise inter-group distance calculation -----------------------------
# Jensen-Shannon Divergence
JSD <- function(InMat, pseudocount = 0.000001) {
    kld <- function(x, y) sum(x * log(x / y))
    jsd <- function(x, y) sqrt(0.5 * kld(x, (x + y) / 2) + 0.5 * kld(y, (x + y) / 2))
    ncol1 <- ncol(InMat)
    colname <- colnames(InMat)
    resultMatrix <- matrix(0, ncol1, ncol1)
    InMat <- apply(InMat, 1:2, function(x) ifelse(x == 0, pseudocount, x))
    for (i in 1:ncol1) {
        for (j in 1:ncol1) {
            resultMatrix[i, j] <- jsd(as.vector(InMat[, i]), as.vector(InMat[, j]))
        }
    }
    colnames(resultMatrix) <- colname
    rownames(resultMatrix) <- colname
    resultMatrix <- as.dist(resultMatrix)
    attr(resultMatrix, "method") <- "dist"
    return(resultMatrix)
}

profile_rel <- prop.table(feature_tab, 2) * 100

Jsd_mat <- JSD(profile_rel) %>%
    as.matrix() %>%
    melt()
profile_df <- as.data.frame(t(profile_rel))
Bc_mat <- vegdist(profile_df, method = "bray") %>%
    as.matrix() %>%
    melt()
Jaccard_mat <- vegdist(profile_df, method = "jaccard") %>%
    as.matrix() %>%
    melt()

dist <- merge(Bc_mat, Jaccard_mat, by = c("Var1", "Var2")) %>%
    merge(Jsd_mat, by = c("Var1", "Var2"))
colnames(dist) <- c("Sample1", "Sample2", "BC", "jaccard", "JSD")
dist$Sample1 <- as.character(dist$Sample1)
dist$Sample2 <- as.character(dist$Sample2)
dist <- dist[dist$Sample1 < dist$Sample2, ]

# --- 5. Assign group labels to pairwise distances ----------------------------
dat2 <- meta %>%
    tibble::rownames_to_column(var = "Sample") %>%
    select(Sample, Group) %>%
    group_split(Group)

dist <- dist %>%
    mutate(Group = case_when(
        Sample1 %in% dat2[[1]]$Sample & Sample2 %in% dat2[[1]]$Sample ~
            paste(unique(dat2[[1]]$Group), unique(dat2[[1]]$Group), sep = "-"),
        Sample1 %in% dat2[[2]]$Sample & Sample2 %in% dat2[[2]]$Sample ~
            paste(unique(dat2[[2]]$Group), unique(dat2[[2]]$Group), sep = "-"),
        Sample1 %in% dat2[[3]]$Sample & Sample2 %in% dat2[[3]]$Sample ~
            paste(unique(dat2[[3]]$Group), unique(dat2[[3]]$Group), sep = "-"),
        Sample1 %in% dat2[[1]]$Sample & Sample2 %in% dat2[[2]]$Sample ~
            paste(unique(dat2[[1]]$Group), unique(dat2[[2]]$Group), sep = "-"),
        Sample1 %in% dat2[[1]]$Sample & Sample2 %in% dat2[[3]]$Sample ~
            paste(unique(dat2[[1]]$Group), unique(dat2[[3]]$Group), sep = "-"),
        Sample1 %in% dat2[[2]]$Sample & Sample2 %in% dat2[[1]]$Sample ~
            paste(unique(dat2[[2]]$Group), unique(dat2[[1]]$Group), sep = "-"),
        Sample1 %in% dat2[[2]]$Sample & Sample2 %in% dat2[[3]]$Sample ~
            paste(unique(dat2[[2]]$Group), unique(dat2[[3]]$Group), sep = "-"),
        Sample1 %in% dat2[[3]]$Sample & Sample2 %in% dat2[[1]]$Sample ~
            paste(unique(dat2[[3]]$Group), unique(dat2[[1]]$Group), sep = "-"),
        Sample1 %in% dat2[[3]]$Sample & Sample2 %in% dat2[[2]]$Sample ~
            paste(unique(dat2[[3]]$Group), unique(dat2[[2]]$Group), sep = "-")
    ))

write.csv(dist, file = here("outputs", "beta_diversity", "dist_tab.csv"), row.names = FALSE)

# --- 6. Figure 1B: Between-group Bray-Curtis distance boxplot ----------------
data_subset <- dist %>%
    filter(Group %in% c("Relatives-Control", "Relatives-Psoriasis", "Control-Psoriasis"))

pal <- c("#8b4d73", "#4760a3", "#8a6f8c")

p_beta <- data_subset %>%
    ggplot(aes(x = Group, y = BC)) +
    geom_violin(aes(fill = Group), alpha = 0.35, color = NA) +
    geom_boxplot(aes(fill = Group), width = 0.15, outlier.shape = NA) +
    scale_fill_manual(values = pal, guide = "none") +
    labs(y = "Bray-Curtis Distance") +
    theme_classic()

ggsave(here("outputs", "beta_diversity", "beta_diversity_BC_boxplot.pdf"), plot = p_beta)

# --- 7. Hierarchical clustering dendrogram + stacked barplot (Figure 1a) ------
library(ggtree)
library(ggstance)
library(ggsci)
library(scales)
library(ggClusterNet)

ps0 <- readRDS(here("data", "metaphlan", "phyloseq_taxa.rds"))

# Assign Group labels
ps <- ps0
sample_data(ps)$Group <- case_when(
    sample_data(ps)$Disease == "Psoriasis" ~ "Psoriasis",
    sample_data(ps)$Disease == "Healthy" & sample_data(ps)$Dataset == "Train" ~ "Relatives",
    sample_data(ps)$Disease == "Healthy" & sample_data(ps)$Dataset == "Validate" ~ "Control"
)

ps1_rela <- transform_sample_counts(ps, function(x) x / sum(x))

# Aggregate OTU by group
otu <- as.data.frame(ggClusterNet::vegan_otu(ps1_rela))
group.split <- split(otu, as.factor(sample_data(ps1_rela)$Group))
group.apply <- lapply(group.split, colMeans)
group.combine <- do.call(rbind, group.apply)
otuG <- t(group.combine)

ps_grp <- phyloseq(otu_table(otuG, taxa_are_rows = TRUE), tax_table(ps1_rela))

# Hierarchical clustering
hc <- ps_grp %>%
    distance(method = "bray") %>%
    hclust(method = "complete")

clus <- cutree(hc, k = 3)
d <- data.frame(label = names(clus), member = factor(clus))
map <- data.frame(
    ID = unique(sample_data(ps1_rela)$Group),
    row.names = unique(sample_data(ps1_rela)$Group),
    Group = unique(sample_data(ps1_rela)$Group)
)
dd <- merge(d, map, by = "row.names", all = FALSE)
rownames(dd) <- dd$Row.names
dd$Row.names <- NULL

p1 <- ggtree(hc) %<+% dd +
    geom_tippoint(size = 5, shape = 21, aes(fill = member, x = x)) +
    geom_tiplab(aes(color = member, x = x * 1.2), hjust = 1) +
    scale_color_manual(values = c("steelblue", "indianred", "darkslateblue"))

# Stacked barplot at Phylum level
j <- "Phylum"
Top <- 6
psdata <- ggClusterNet::tax_glom_wt(ps = ps1_rela, ranks = j)
psdata <- psdata %>% transform_sample_counts(function(x) x / sum(x))

otu_bar <- otu_table(psdata)
tax_bar <- tax_table(psdata)

for (i in 1:nrow(tax_bar)) {
    if (rownames(tax_bar)[i] %in% names(sort(rowSums(otu_bar), decreasing = TRUE)[1:Top])) {
        tax_bar[i, j] <- tax_bar[i, j]
    } else {
        tax_bar[i, j] <- "Other"
    }
}
tax_table(psdata) <- tax_bar

Taxonomies <- psdata %>% psmelt()
grotax <- Taxonomies %>%
    group_by(Group, !!sym(j)) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop")

table_data <- data.frame()
for (i in seq_along(unique(sample_data(psdata)$Group))) {
    a <- as.data.frame(table(sample_data(psdata)$Group))[i, 1]
    b <- as.data.frame(table(sample_data(psdata)$Group))[i, 2]
    c_data <- grotax %>% filter(Group == a)
    c_data$Abundance <- c_data$Abundance / b
    table_data <- rbind(table_data, c_data)
}

pal_d3 <- pal_d3("category20")(14)
color <- c("steelblue", "indianred", "darkslateblue", rev(pal_d3))

p2 <- facet_plot(p1,
    panel = "Stacked Barplot",
    data = table_data,
    geom = ggstance::geom_barh,
    mapping = aes(x = Abundance, fill = !!sym(j)),
    stat = "identity"
) +
    scale_fill_manual(values = color)

ggsave(p2,
    filename = here("outputs", "beta_diversity", "beta_diversity_cluster_barplot.pdf"),
    width = 25, height = 6
)
