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
ps0 <- readRDS(here("data", "metaphlan", "phyloseq_taxa.rds"))

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

# --- 3. Custom network function (modified from ggClusterNet::network.bindif) --
#   Customised color scheme for phylum-level network plots
network.i.mod <- function(otu = NULL, tax = NULL, map = NULL, ps = NULL, N = 0,
                          big = FALSE, select_layout = FALSE, layout_net = "model_igraph2",
                          r.threshold = 0.6, p.threshold = 0.05, maxnode = 2, method = "spearman",
                          label = FALSE, lab = "elements", group = "Group", path = "./",
                          fill = "Phylum", size = "igraph.degree", scale = TRUE, zipi = FALSE,
                          clu_method = "cluster_fast_greedy", step = 100, yourmem = theme_void(),
                          ncol = 3, nrow = 1, R = 10, ncpus = 1, order = NULL) {
    ps <- inputMicro(otu, tax, map, tree, ps, group = group)
    if (scale) {
        ps_rela <- scale_micro(ps = ps, method = "rela")
    } else {
        ps_rela <- ps
    }
    mapping <- as.data.frame(sample_data(ps_rela))
    y <- matrix(1, nrow = 16, ncol = length(unique(mapping$Group)))
    if (is.null(order)) {
        layouts <- as.character(unique(mapping$Group))
    } else {
        layouts <- order
    }
    mapping$ID <- row.names(mapping)
    plots <- list()
    plots1 <- list()
    plots2 <- list()
    aa <- 1
    for (layout in layouts) {
        mapi <- mapping[mapping$Group == layout, ]
        psi <- phyloseq(
            otu_table(ps_rela), tax_table(ps_rela),
            sample_data(mapi)
        ) %>%
            filter_OTU_ps(Top = N) %>%
            filter_taxa(function(x) sum(x) > 0, TRUE)
        print(layout)
        if (big == TRUE) {
            result <- cor_Big_micro2(
                ps = psi, N = 0, r.threshold = r.threshold,
                p.threshold = p.threshold, method = method, scale = FALSE
            )
            a <- 2
        } else if (big == FALSE) {
            result <- corMicro(
                ps = psi, N = 0, r.threshold = r.threshold,
                p.threshold = p.threshold, method = method, R = R,
                ncpus = ncpus
            )
            a <- 1
        }
        print("cor matrix culculating over")
        cor <- result[[1]]
        if (cor %>% as.vector() %>% max() == 0) {
            stop("The connect value in cor matrix all was zone")
        }
        if (FALSE) {
            node <- NULL
            node <- culculate_node_axis(
                cor.matrix = cor, layout = layout_net,
                seed = 1, group = NULL, model = FALSE, method = clu_method
            )
        } else if (layout_net == "model_Gephi.2") {
            result2 <- model_Gephi.2(
                cor = cor, method = clu_method,
                seed = 12
            )
            node <- result2[[1]]
        } else if (layout_net == "model_igraph2") {
            result2 <- model_igraph2(
                cor = cor, method = clu_method,
                seed = 12
            )
            node <- result2[[1]]
            head(node)
            dat <- result2[[2]]
            head(dat)
            tem <- data.frame(mod = dat$model, col = dat$color) %>%
                dplyr::distinct(mod, .keep_all = TRUE)
            col <- tem$col
            names(col) <- tem$mod
        }
        otu_table <- as.data.frame(t(vegan_otu(psi)))
        tax_table <- as.data.frame(vegan_tax(psi))
        nodes <- nodeadd(
            plotcord = node, otu_table = otu_table,
            tax_table = tax_table
        )
        edge <- edgeBuild(cor = cor, node = node)
        if (layout_net == "model_igraph2") {
            tem2 <- dat %>%
                dplyr::select(OTU, model, color) %>%
                dplyr::right_join(edge, by = c(OTU = "OTU_1")) %>%
                dplyr::rename(OTU_1 = OTU, model1 = model, color1 = color)
            head(tem2)
            tem3 <- dat %>%
                dplyr::select(OTU, model, color) %>%
                dplyr::right_join(edge, by = c(OTU = "OTU_2")) %>%
                dplyr::rename(OTU_2 = OTU, model2 = model, color2 = color)
            head(tem3)
            tem4 <- tem2 %>% inner_join(tem3)
            head(tem4)
            edge2 <- tem4 %>% mutate(color = ifelse(model1 ==
                model2, as.character(model1), "across"), manual = ifelse(model1 ==
                model2, as.character(color1), "#C1C1C1"))
            head(edge2)
            col_edge <- edge2 %>%
                as.tibble() %>%
                dplyr::distinct(color,
                    .keep_all = TRUE
                ) %>%
                dplyr::select(color, manual)
            col0 <- col_edge$manual
            names(col0) <- col_edge$color
            p1 <- ggplot() +
                geom_segment(
                    aes(
                        x = X1, y = Y1,
                        xend = X2, yend = Y2, color = color
                    ),
                    data = edge2,
                    size = 1
                ) +
                scale_colour_manual(values = col0)
            p2 <- p1 + new_scale_color() + geom_point(aes(X1,
                X2,
                color = model
            ), data = dat, size = 4) + scale_colour_manual(values = col0) +
                scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
                theme(panel.background = element_blank()) + theme(
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank()
                ) + theme(legend.background = element_rect(colour = NA)) +
                theme(panel.background = element_rect(
                    fill = "white",
                    colour = NA
                )) + theme(
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank()
                ) + labs(title = paste(layout,
                    "links:", dim(edge2)[1], "nodes:", dim(dat)[1],
                    sep = ""
                ))
            p2
            plotname <- paste(path, "/network_igraph2", layout,
                ".pdf",
                sep = ""
            )
            ggsave(plotname, p2, width = 16, height = 14)
        }
        edge_Gephi <- data.frame(
            source = edge$OTU_1, target = edge$OTU_2,
            correlation = edge$weight, direct = "undirected",
            cor = edge$cor
        )
        node_Gephi <- data.frame(
            ID = nodes$elements, nodes[4:dim(nodes)[2]],
            Label = nodes$elements
        )
        idedge <- c(as.character(edge_Gephi$source), as.character(edge_Gephi$target))
        idedge <- unique(idedge)
        row.names(node_Gephi) <- as.character(node_Gephi$ID)
        node_Gephi1 <- node_Gephi[idedge, ]
        write.csv(edge_Gephi, paste(path, "/", layout, "_Gephi_edge.csv",
            sep = ""
        ), row.names = FALSE, quote = FALSE)
        write.csv(node_Gephi, paste(path, "/", layout, "_Gephi_allnode.csv",
            sep = ""
        ), row.names = FALSE, quote = FALSE)
        write.csv(node_Gephi1, paste(path, "/", layout, "_Gephi_edgenode.csv",
            sep = ""
        ), row.names = FALSE, quote = FALSE)
        igraph <- igraph::graph_from_data_frame(nodeEdge(cor = cor)[[1]],
            directed = FALSE, vertices = nodeEdge(cor = cor)[[2]]
        )
        nodepro <- node_properties(igraph)
        write.csv(nodepro, paste(path, "/", layout, "_node_properties.csv",
            sep = ""
        ), row.names = TRUE)
        nodeG <- merge(nodes, nodepro, by = "row.names", all.x = TRUE)
        row.names(nodeG) <- nodeG$Row.names
        nodeG$Row.names <- NULL
        numna <- (dim(nodeG)[2] - 3):dim(nodeG)[2]
        nodeG[, numna][is.na(nodeG[, numna])] <- 0
        head(nodeG)
        pnet <- ggplot() +
            geom_segment(
                data = edge, aes(
                    x = X1, y = Y1, xend = X2,
                    yend = Y2, color = cor
                ), size = 0.5,
                alpha = 0.8
            ) +
            geom_point(data = nodeG, aes(X1, X2,
                fill = !!sym(fill),
                size = !!sym(size)
            ), pch = 21, color = "gray40") +
            labs(title = paste(layout, "network", sep = "_")) +
            scale_colour_manual(values = c("#4b0082", "#ffd700")) +
            scale_fill_manual(values = c(
                "p__Actinobacteria" = "#f8766d",
                "p__Bacteroidetes" = "#e68613",
                "p__Firmicutes" = "#00c19a",
                "p__Proteobacteria" = "#00a9ff",
                "Others" = "#7e7e7e"
            )) +
            scale_size(range = c(0.8, maxnode)) +
            scale_x_continuous(breaks = NULL) +
            scale_y_continuous(breaks = NULL) +
            theme(
                panel.background = element_blank(),
                plot.title = element_text(hjust = 0.5)
            ) +
            theme(
                axis.title.x = element_blank(),
                axis.title.y = element_blank()
            ) +
            theme(legend.background = element_rect(colour = NA)) +
            theme(panel.background = element_rect(
                fill = "white",
                colour = NA
            )) +
            theme(
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank()
            )
        pnet1 <- ggplot() +
            geom_curve(
                data = edge, aes(
                    x = X1, y = Y1,
                    xend = X2, yend = Y2, color = cor
                ),
                size = 0.5, alpha = 0.5, curvature = -0.2
            ) +
            geom_point(
                data = nodeG, aes(X1, X2, fill = !!sym(fill), size = !!sym(size)),
                pch = 21, color = "gray40"
            ) +
            labs(title = paste(layout, "network", sep = "_")) +
            scale_colour_manual(values = c("#4b0082", "#ffd700")) +
            scale_fill_manual(values = c(
                "p__Actinobacteria" = "#f8766d",
                "p__Bacteroidetes" = "#e68613",
                "p__Firmicutes" = "#00c19a",
                "p__Proteobacteria" = "#00a9ff",
                "Others" = "#7e7e7e"
            )) +
            scale_size(range = c(0.8, maxnode)) +
            scale_x_continuous(breaks = NULL) +
            scale_y_continuous(breaks = NULL) +
            theme(panel.background = element_blank(), plot.title = element_text(hjust = 0.5)) +
            theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
            theme(legend.background = element_rect(colour = NA)) +
            theme(panel.background = element_rect(
                fill = "white",
                colour = NA
            )) +
            theme(
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank()
            )
        pnet1
        if (label == TRUE) {
            pnet <- pnet + ggrepel::geom_text_repel(aes(X1, X2,
                label = !!sym(lab)
            ), size = 4, data = nodeG)
            pnet1 <- pnet1 + ggrepel::geom_text_repel(aes(X1,
                X2,
                label = !!sym(lab)
            ), size = 4, data = nodeG)
        }
        plotname <- paste(path, "/network", layout, ".pdf", sep = "")
        ggsave(plotname, pnet, width = 11, height = 9)
        plotname <- paste(path, "/network", layout, "_cover.pdf",
            sep = ""
        )
        ggsave(plotname, pnet1, width = 11, height = 9)
        plots[[aa]] <- pnet
        plots1[[aa]] <- pnet1
        if (layout_net == "model_igraph2") {
            plots2[[aa]] <- p2
        }
        if (zipi) {
            res <- ZiPiPlot(igraph = igraph, method = clu_method)
            p <- res[[1]]
            ggsave(paste(path, "/", layout, "_ZiPi.pdf", sep = ""),
                p,
                width = 12, height = 10
            )
            ZiPi <- res[[2]]
            write.csv(ZiPi, paste(path, "/", layout, "ZiPi.csv",
                sep = ""
            ), row.names = FALSE)
        }
        netpro_result <- net_properties.2(igraph)
        colnames(netpro_result) <- layout
        y <- as.data.frame(y)
        colnames(y) <- layouts
        y[layout] <- netpro_result[, 1]
        row.names(y) <- row.names(netpro_result)
        aa <- aa + 1
    }
    plotname <- paste(path, "/network_all.pdf", sep = "")
    p <- ggpubr::ggarrange(
        plotlist = plots, common.legend = TRUE,
        legend = "right", ncol = ncol, nrow = nrow
    )
    p1 <- ggpubr::ggarrange(
        plotlist = plots1, common.legend = TRUE,
        legend = "right", ncol = ncol, nrow = nrow
    )
    if (layout_net == "model_igraph2") {
        p2 <- ggpubr::ggarrange(
            plotlist = plots2, common.legend = TRUE,
            legend = "right", ncol = ncol, nrow = nrow
        )
    }
    if (length(layouts) == 1) {
        p <- pnet
        p1 <- pnet1
    }
    return(list(p, y, p1, cor, p2))
}

# --- 4. Data transformation --------------------------------------------------
ps.log10p <- microbiome::transform(ps, "log10p")
netpath <- here("outputs", "network_igraph")
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

stabpath <- here("outputs", "network_stab")
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
ggsave(paste0(stabpath, "/vulnerability.pdf"), p_vul, width = 8, height = 4)

message("03_Network_Analysis.R completed successfully.")
