#!/usr/bin/env Rscript
# scRNAseq_pipeline_seurat_GPTCelltype.R
# Purpose: Clean, analyze and annotate 10X single-cell RNA-seq data with Seurat + GPTCelltype.
# Author: Bruno Marques — adapted for publication on GitHub
# Date: 2025-11-03
# NOTE: This script is formatted for clarity and reproducibility. 

#########################
# 0. Setup / arguments  #
#########################

suppressPackageStartupMessages({
  library(optparse)
  library(glue)
})

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Comma-separated list of 10X H5 files or a directory containing them.", metavar = "character"),
  make_option(c("-n", "--names"), type = "character", default = NULL,
              help = "Comma-separated sample names in the same order as files.", metavar = "character"),
  make_option(c("-o", "--outdir"), type = "character", default = "outputs",
              help = "Output directory [default: %default]"),
  make_option(c("--openai_key_env"), type = "character", default = "OPENAI_API_KEY",
              help = "Environment variable name that stores your OpenAI API key [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input) || is.null(opt$names)) {
  message("Usage: Rscript scRNAseq_pipeline_seurat_GPTCelltype.R -i <file1.h5,file2.h5,...> -n <name1,name2,...> [-o outputs]")
  quit(status = 1)
}

# Prepare paths and variables
outdir <- opt$outdir
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

files <- strsplit(opt$input, ",")[[1]]
sample_names <- strsplit(opt$names, ",")[[1]]
if (length(files) != length(sample_names)) stop("Number of files must match number of sample names")

#########################
# 1. Load libraries     #
#########################
suppressPackageStartupMessages({
  # Core single-cell analysis
  library(Seurat)
  library(hdf5r)   # Read10X_h5 dependency
  library(harmony)
  
  # Data handling / plotting
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  
  # Additional analysis
  library(writexl)
  library(readxl)
  library(pheatmap)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  
  # Cell-cell communication
  library(CellChat)
  
  # GPTCelltype (installed from GitHub). The README shows dependency on openai package.
  # installation: install.packages('openai'); remotes::install_github('Winnie09/GPTCelltype')
  library(GPTCelltype)
})

#########################
# 2. Utility functions  #
#########################

log_msg <- function(...) {
  message(sprintf("[%s] %s", Sys.time(), glue::glue(...)))
}

safe_save_plot <- function(plot_obj, filename, width = 8, height = 4) {
  ggsave(filename = file.path(outdir, filename), plot = plot_obj, width = width, height = height)
}

#########################
# 3. Read 10X H5 files  #
#########################
log_msg("Reading %d files...", length(files))
seurat.list <- lapply(seq_along(files), function(i) {
  counts <- Read10X_h5(files[i])
  so <- CreateSeuratObject(counts = counts, project = sample_names[i])
  so <- RenameCells(so, add.cell.id = sample_names[i])
  so$sample <- sample_names[i]
  return(so)
})
names(seurat.list) <- sample_names

# Merge into a single object
combined <- merge(x = seurat.list[[1]], y = seurat.list[-1], add.cell.ids = sample_names, project = "Project")
combined$group <- ifelse(grepl("Adherent", combined$sample), "Adherent", "Floating")

#########################
# 4. QC and filtering   #
#########################
log_msg("Calculating QC metrics and filtering...")
combined["percent.mt"] <- PercentageFeatureSet(combined, pattern = "^MT-")
combined["percent.ribo"] <- PercentageFeatureSet(combined, pattern = "^RPS|^RPL")

# Basic filtering thresholds (adjustable)
min_features <- 500
max_features <- 7000
max_mito <- 15
combined <- subset(combined, subset = nFeature_RNA > min_features & nFeature_RNA < max_features & percent.mt < max_mito)

# Quick QC violin plots
meta <- combined@meta.data

theme_base <- theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7))

p1 <- ggplot(meta, aes(x = sample, y = nCount_RNA, fill = sample)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.8, color = "black") +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", color = "gray30") +
  geom_hline(yintercept = 1000, color = "red", size = 0.7) +
  scale_y_log10(breaks = 10^seq(1, 5), labels = trans_format("log10", math_format(10^.x))) +
  labs(title = "nCount_RNA", y = "# UMI / célula", x = NULL) + theme_base

p2 <- ggplot(meta, aes(x = sample, y = nFeature_RNA, fill = sample)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.8, color = "black") +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", color = "gray30") +
  geom_hline(yintercept = 200, color = "red", size = 0.7) +
  scale_y_log10(breaks = 10^seq(1, 5), labels = trans_format("log10", math_format(10^.x))) +
  labs(title = "nFeature_RNA", y = "# genes / célula", x = NULL) + theme_base

p3 <- ggplot(meta, aes(x = sample, y = percent.mt, fill = sample)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.8, color = "black") +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", color = "gray30") +
  geom_hline(yintercept = 10, color = "red", size = 0.7) + scale_y_continuous(limits = c(0, NA)) +
  labs(title = "percent.mt", y = "% mitocondrial", x = NULL) + theme_base

p4 <- ggplot(meta, aes(x = sample, y = percent.ribo, fill = sample)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.8, color = "black") +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", color = "gray30") +
  geom_hline(yintercept = 20, color = "red", size = 0.7) + scale_y_continuous(limits = c(0, NA)) +
  labs(title = "percent.ribo", y = "% ribossomal", x = NULL) + theme_base

qc_plot <- (p1 | p2) / (p3 | p4)

safe_save_plot(qc_plot, "QC_violin_minimal.png", width = 12, height = 8)

#########################
# 5. Normalization / PCA / UMAP
#########################
set.seed(123)
combined <- combined %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(verbose = FALSE)

# Run UMAP on PCA
combined <- RunUMAP(combined, dims = 1:20)

# Harmony integration by sample
combined <- RunHarmony(combined, group.by.vars = "sample")
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:30)
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:30)
combined <- FindClusters(combined, resolution = 1)

# UMAP plots: group, sample, clusters
p_umap_group <- DimPlot(combined, group.by = "group", reduction = "umap") + ggtitle("UMAP per group")
p_umap_sample <- DimPlot(combined, group.by = "sample", reduction = "umap") + ggtitle("UMAP per sample")
p_umap_clusters <- DimPlot(combined, group.by = "seurat_clusters", label = TRUE, reduction = "umap") + ggtitle("UMAP per clusters")

safe_save_plot(p_umap_group | p_umap_sample | p_umap_clusters, "UMAP_triplo.png", width = 15, height = 5)

#########################
# 6. Markers: conserved across groups
#########################
log_msg("Finding conserved markers per cluster...")
cluster_ids <- levels(as.factor(Idents(combined)))
conserved_list <- lapply(cluster_ids, function(cluster) {
  FindConservedMarkers(combined, ident.1 = cluster, grouping.var = "group", only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25)
})
names(conserved_list) <- cluster_ids

# Save conserved markers to xlsx files
for (cluster in names(conserved_list)) {
  df <- conserved_list[[cluster]]
  if (nrow(df) > 0) {
    df$gene <- rownames(df)
    write_xlsx(df, file.path(outdir, paste0("Conserved_Markers_cluster_", cluster, ".xlsx")))
  }
}

#########################
# 7. Cell type annotation using GPTCelltype
#########################
# GPTCelltype expects a DataFrame of marker genes (output from FindAllMarkers or FindConservedMarkers).
# IMPORTANT: GPTCelltype uses OpenAI's GPT-4 (or other models) and requires an API key set as an environment variable.
# Set your key before running, e.g.: Sys.setenv(OPENAI_API_KEY = 'sk-...') OR export OPENAI_API_KEY in your shell.

api_env_name <- opt$openai_key_env
if (is.na(Sys.getenv(api_env_name, unset = NA))) {
  log_msg("WARNING: %s is not set. GPTCelltype will fail without an OpenAI API key. Set it and re-run.", api_env_name)
} else {
  log_msg("Running GPTCelltype to annotate clusters (this will call OpenAI API: costs may apply)...")
  
  # Prepare a markers table to feed GPTCelltype. The package vignette suggests passing FindAllMarkers output.
  # We'll build a 'markers' table that contains top markers per cluster (conserved across groups if present,
  # otherwise use FindAllMarkers).
  
  # Combine conserved markers into a single data.frame with cluster column
  combined_markers <- do.call(rbind, lapply(names(conserved_list), function(cl) {
    df <- conserved_list[[cl]]
    if (nrow(df) == 0) return(NULL)
    df$cluster <- cl
    df$gene <- rownames(df)
    return(df)
  }))
  
  # If no conserved markers (unlikely), fallback to FindAllMarkers
  if (is.null(combined_markers) || nrow(combined_markers) == 0) {
    log_msg("No conserved markers detected; falling back to FindAllMarkers() to generate markers for GPTCelltype.")
    all_markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    all_markers$cluster <- all_markers$cluster
    all_markers$gene <- rownames(all_markers)
    markers_for_gpt <- all_markers
  } else {
    markers_for_gpt <- combined_markers
  }
  
  # Optionally keep top N markers per cluster to reduce cost and improve clarity
  topN <- 20
  markers_for_gpt <- markers_for_gpt %>%
    group_by(cluster) %>%
    arrange(cluster, desc(avg_log2FC)) %>%
    slice_head(n = topN) %>%
    ungroup()
  
  # Run GPTCelltype: this returns a named vector/list where names correspond to cluster ids
  # The 'tissuename' parameter can improve accuracy; change if appropriate
  Sys.setenv(OPENAI_API_KEY = Sys.getenv(api_env_name))
  gpt_res <- tryCatch({
    gptcelltype(markers_for_gpt, model = 'gpt-4', tissuename = 'olfactory_epithelium')
  }, error = function(e) {
    log_msg("GPTCelltype failed: %s", e$message)
    return(NULL)
  })
  
  if (!is.null(gpt_res)) {
    # Map GPTCelltype annotation back to cells in the Seurat object
    # The package vignette suggests: obj@meta.data$celltype <- as.factor(res[as.character(Idents(obj))])
    ann_vector <- gpt_res[as.character(Idents(combined))]
    combined@meta.data$celltype <- as.factor(ann_vector)
    
    # Save annotation table
    ann_df <- data.frame(cell = colnames(combined), celltype = combined@meta.data$celltype, sample = combined@meta.data$sample)
    write_xlsx(ann_df, file.path(outdir, "GPTCelltype_annotation_per_cell.xlsx"))
    
    # UMAP plot colored by celltype
    p_celltype_umap <- DimPlot(combined, reduction = "umap", group.by = "celltype", label = TRUE, label.size = 3)
    safe_save_plot(p_celltype_umap, "Umap_celltype.png", width = 8, height = 6)
    
    log_msg("GPTCelltype annotation finished and saved.")
  }
}

#########################
# 8. Additional analyses
#########################
# The rest of the analyses (dotplots, differential testing between groups within a cell type,
# enrichment analysis, CellChat workflows, etc.) can be modularized into functions and
# called here. For brevity we provide a minimal example of differential testing for one cell type.

# Example: differential testing in 'Immature Neurons' between Adherent and Floating (if present)
if ("celltype" %in% colnames(combined@meta.data)) {
  if ("Immature Neurons" %in% levels(as.factor(combined@meta.data$celltype))) {
    combined_subset <- subset(combined, subset = celltype == "Immature Neurons")
    Idents(combined_subset) <- "group"
    markers_in <- FindMarkers(combined_subset, ident.1 = "Adherent", ident.2 = "Floating")
    markers_in$gene <- rownames(markers_in)
    write_xlsx(markers_in, file.path(outdir, "DE_genes_ImmatureNeurons.xlsx"))
  }
}


#########################
# 9. Visualization: UMAP faceted by sample and celltype
#########################
log_msg("Creating UMAP faceted plots by sample and celltype...")

# Ensure celltype exists
if ("celltype" %in% colnames(combined@meta.data)) {
  umap_df <- as.data.frame(Embeddings(combined, "umap"))
  umap_df$celltype <- combined@meta.data$celltype
  umap_df$sample <- combined@meta.data$sample
  
  # If celltype is a factor, preserve levels
  if (!is.factor(umap_df$celltype)) umap_df$celltype <- factor(umap_df$celltype)
  
  # Faceted UMAP by sample colored by celltype
  p_umap_facets <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = celltype)) +
    geom_point(size = 0.4, alpha = 0.8) +
    facet_wrap(~ sample, ncol = 2) +
    theme_classic() +
    labs(color = "Cell Type") +
    theme(strip.text = element_text(size = 12, face = "bold"))
  
  safe_save_plot(p_umap_facets, "Umap_celltype_patients.png", width = 10, height = 6)
}

#########################
# 10. Feature plots for specific markers (subset approach)
#########################
log_msg("Generating feature plots for selected markers...")

DefaultAssay(combined) <- "RNA"
markers_feature <- c("TUBB3", "SOX2", "CALD1", "CTSV")
for (m in markers_feature) {
  p <- tryCatch({
    FeaturePlot(combined, features = m, reduction = "umap", split.by = "group")
  }, error = function(e) {
    log_msg("FeaturePlot failed for %s: %s", m, e$message)
    return(NULL)
  })
  if (!is.null(p)) safe_save_plot(p, paste0("featureplot_", m, ".png"), width = 10, height = 4)
}

#########################
# 11. FindAllMarkers / DotPlots
#########################
log_msg("Finding all markers and creating dotplots...")
all_markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write_xlsx(all_markers, file.path(outdir, "FindAllMarkers_all_clusters.xlsx"))

# Top 5 per cluster
library(dplyr)
top5_per_cluster <- all_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

# Remove unwanted genes
genes_excluir <- c("ENSG00000289901", "MT-ND2", "MT-ND1", "MT-ND5", "MT-ND4L", "MT-ND6", "LINC00472")
genes_filtrados <- setdiff(unique(top5_per_cluster$gene), genes_excluir)

dotplottotal <- DotPlot(combined, features = genes_filtrados, group.by = "celltype") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
safe_save_plot(dotplottotal, "dotplot_total.png", width = 14, height = 5)

dotplotgrupo <- DotPlot(combined, features = genes_filtrados, split.by = "group", group.by = "celltype") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
safe_save_plot(dotplotgrupo, "dotplot_group.png", width = 16, height = 6)

#########################
# 12. Proportional testing with speckle (propeller)
#########################
log_msg("Running propeller (speckle) to test cell type proportions between groups...")

if (!requireNamespace("speckle", quietly = TRUE)) {
  log_msg("Package 'speckle' is not installed. Skipping propeller analysis.")
} else {
  library(speckle)
  # Build groups mapping from 'sample' to 'group' if necessary
  sample_levels <- unique(combined@meta.data$sample)
  # If the user provided a different mapping, adapt here. We'll attempt to use current combined metadata.
  combined@meta.data$groups <- combined@meta.data$group
  
  celltypes <- factor(combined@meta.data$celltype)
  samples <- factor(combined@meta.data$sample)
  groups <- factor(combined@meta.data$groups)
  
  propeller_results <- tryCatch({
    speckle::propeller(clusters = celltypes, sample = samples, group = groups, transform = "asin")
  }, error = function(e) {
    log_msg("propeller failed: %s", e$message)
    return(NULL)
  })
  
  if (!is.null(propeller_results)) {
    write_xlsx(as.data.frame(propeller_results), file.path(outdir, "propeller_results.xlsx"))
    
    # Create heatmap of mean proportions
    heatmap_df <- propeller_results %>% dplyr::select(BaselineProp.clusters, starts_with("PropMean."))
    if (ncol(heatmap_df) >= 2) {
      heatmap_matrix <- heatmap_df %>% column_to_rownames(var = "BaselineProp.clusters") %>% as.matrix()
      pheatmap::pheatmap(heatmap_matrix, cluster_rows = FALSE, cluster_cols = FALSE,
                         filename = file.path(outdir, "heatmap_propeller.png"), width = 8, height = 6)
    }
    
    # Plot proportions (speckle helper)
    prop_plot <- tryCatch({
      speckle::plotCellTypeProps(clusters = celltypes, sample = groups) + theme_classic()
    }, error = function(e) NULL)
    if (!is.null(prop_plot)) safe_save_plot(prop_plot, "proportion_celltype.png", width = 8, height = 4)
  }
}

#########################
# 13. DE within a cell type (example: Immature Neurons) + Volcano
#########################
log_msg("Differential expression within Immature Neurons and volcano plotting...")
if ("celltype" %in% colnames(combined@meta.data) && "Immature Neurons" %in% combined@meta.data$celltype) {
  combined_subset <- subset(combined, subset = celltype == "Immature Neurons")
  # Ensure the grouping variable is set as identity
  Idents(combined_subset) <- combined_subset$group
  
  markers3 <- FindMarkers(combined_subset, ident.1 = "Adherent", ident.2 = "Floating")
  if (nrow(markers3) > 0) {
    markers3$gene <- rownames(markers3)
    markers3 <- markers3[!grepl("^ENSG", markers3$gene), ]
    write_xlsx(markers3, file.path(outdir, "DE_genes_ImmatureNeurons.xlsx"))
    
    # Volcano plot using EnhancedVolcano
    if (requireNamespace("EnhancedVolcano", quietly = TRUE)) {
      library(EnhancedVolcano)
      volcano <- EnhancedVolcano(markers3, lab = markers3$gene, x = 'avg_log2FC', y = 'p_val_adj',
                                 title = 'Adherent vs Floating (Immature Neurons)', pCutoff = 0.01, FCcutoff = 1,
                                 pointSize = 2.0, labSize = 3.0)
      safe_save_plot(volcano, "Volcano_Immature_Neurons.png", width = 12, height = 8)
    } else {
      log_msg("EnhancedVolcano not available; skipping volcano plot.")
    }
    
    # Top genes for dotplot
    top20_adherent <- markers3 %>% filter(avg_log2FC > 0) %>% arrange(desc(avg_log2FC)) %>% slice_head(n = 20)
    top20_floating <- markers3 %>% filter(avg_log2FC < 0) %>% arrange(avg_log2FC) %>% slice_head(n = 20)
    top_genes <- unique(c(rownames(top20_adherent), rownames(top20_floating)))
    
    # DotPlot for Immature Neurons
    dp_in <- tryCatch({
      DotPlot(combined_subset, features = top_genes, group.by = "group") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }, error = function(e) {
      log_msg("DotPlot failed: %s", e$message); NULL
    })
    if (!is.null(dp_in)) safe_save_plot(dp_in, "dotplot_ImmatureNeurons.png", width = 12, height = 6)
  }
}

#########################
# 14. Functional enrichment (GO + Reactome) for significant genes
#########################
log_msg("Running functional enrichment analyses (GO and Reactome)...")
if (exists("markers3") && nrow(markers3) > 0) {
  gene_list <- rownames(markers3[markers3$p_val_adj < 0.05 & markers3$avg_log2FC > 0.1, ])
  if (length(gene_list) > 0) {
    ego_BP <- tryCatch({
      enrichGO(gene = gene_list, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)
    }, error = function(e) { log_msg("enrichGO BP failed: %s", e$message); NULL })
    
    if (!is.null(ego_BP) && nrow(ego_BP) > 0) {
      p1 <- dotplot(ego_BP)
      safe_save_plot(p1, "dotplot_enrichment_ImmatureNeurons_BP.png", width = 12, height = 8)
    }
    
    ego_CC <- tryCatch({
      enrichGO(gene = gene_list, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05)
    }, error = function(e) { log_msg("enrichGO CC failed: %s", e$message); NULL })
    
    if (!is.null(ego_CC) && nrow(ego_CC) > 0) {
      p2 <- dotplot(ego_CC, showCategory = 20)
      safe_save_plot(p2, "dotplot_enrichment_ImmatureNeurons_CC.png", width = 12, height = 8)
    }
    
    # Reactome
    gene_entrez <- tryCatch({ bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID }, error = function(e) NULL)
    if (!is.null(gene_entrez) && length(gene_entrez) > 0) {
      ereact <- tryCatch({ enrichPathway(gene = gene_entrez, organism = "human", pvalueCutoff = 0.05, readable = TRUE) }, error = function(e) NULL)
      if (!is.null(ereact) && nrow(ereact@result) > 0) {
        p_react <- dotplot(ereact, showCategory = 10) + ggtitle("Reactome Pathways Enrichment")
        safe_save_plot(p_react, "reactome_immature_neurons.png", width = 12, height = 9)
        
        # Termsim and network plots
        ereact <- tryCatch({ pairwise_termsim(ereact) }, error = function(e) { log_msg("pairwise_termsim failed: %s", e$message); NULL })
        if (!is.null(ereact)) {
          tryCatch({ emapplot(ereact); ggsave(file.path(outdir, "emapplot_reactome.png"), width = 10, height = 8) }, error = function(e) NULL)
          tryCatch({ cnetplot(ereact, showCategory = 10); ggsave(file.path(outdir, "cnet_reactome.png"), width = 12, height = 8) }, error = function(e) NULL)
        }
        
        # Export term descriptions
        desc_df <- data.frame(Description = ereact@result$Description)
        write_xlsx(desc_df, file.path(outdir, "ereact_Description.xlsx"))
      }
    }
  } else {
    log_msg("No significant genes found for enrichment (gene_list length = 0).")
  }
}

#########################
# 15. CellChat pipeline (cell-cell communication)
#########################
log_msg("Running CellChat analysis (this can be memory/time intensive)...")

if (requireNamespace("CellChat", quietly = TRUE)) {
  library(CellChat)
  tryCatch({
    # Subset by group
    adherent <- subset(combined, subset = group == "Adherent")
    floating <- subset(combined, subset = group == "Floating")
    
    # Prepare data and meta
    data.adherent <- GetAssayData(adherent, assay = "RNA", slot = "data")
    meta.adherent <- data.frame(labels = adherent$celltype, row.names = colnames(adherent))
    
    data.floating <- GetAssayData(floating, assay = "RNA", slot = "data")
    meta.floating <- data.frame(labels = floating$celltype, row.names = colnames(floating))
    
    cellchat.adherent <- createCellChat(object = data.adherent, meta = meta.adherent, group.by = "labels")
    cellchat.floating <- createCellChat(object = data.floating, meta = meta.floating, group.by = "labels")
    
    process_cellchat <- function(cellchat) {
      cellchat@DB <- CellChatDB.human
      cellchat <- subsetData(cellchat)
      cellchat <- identifyOverExpressedGenes(cellchat)
      cellchat <- identifyOverExpressedInteractions(cellchat)
      cellchat <- computeCommunProb(cellchat)
      cellchat <- filterCommunication(cellchat, min.cells = 10)
      cellchat <- computeCommunProbPathway(cellchat)
      cellchat <- aggregateNet(cellchat)
      return(cellchat)
    }
    
    cellchat.adherent <- process_cellchat(cellchat.adherent)
    cellchat.floating <- process_cellchat(cellchat.floating)
    
    # Merge and compare
    object.list <- list(Adherent = cellchat.adherent, Floating = cellchat.floating)
    cellchat.merged <- mergeCellChat(object.list, add.names = names(object.list))
    
    # Number of interactions comparison
    png(file.path(outdir, "number_of_interactions.png"), width = 800, height = 600)
    compareInteractions(cellchat.merged, show.legend = FALSE, group = c(1,2))
    dev.off()
    
    # Circle plots
    png(file.path(outdir, "visual_circle_adherent.png"), width = 800, height = 800)
    netVisual_circle(cellchat.adherent@net$count, vertex.weight = TRUE, title.name = "Adherent", vertex.size = 25, weight.scale = TRUE)
    dev.off()
    
    png(file.path(outdir, "visual_circle_floating.png"), width = 800, height = 800)
    netVisual_circle(cellchat.floating@net$count, vertex.weight = TRUE, title.name = "Floating", vertex.size = 25, weight.scale = TRUE)
    dev.off()
    
    # Assign UMAP coordinates to CellChat objects
    umap_all <- Embeddings(combined, "umap")
    cellchat.adherent@dr$UMAP <- umap_all[colnames(adherent), , drop = FALSE]
    cellchat.floating@dr$UMAP <- umap_all[colnames(floating), , drop = FALSE]
    
    # Centrality heatmaps
    cellchat.adherent <- netAnalysis_computeCentrality(cellchat.adherent)
    ht <- netAnalysis_signalingRole_heatmap(cellchat.adherent)
    png(file.path(outdir, "signalingrole_heatmap_adherent.png"), width = 2500, height = 2000, res = 300)
    ComplexHeatmap::draw(ht)
    dev.off()
    
    cellchat.floating <- netAnalysis_computeCentrality(cellchat.floating)
    ht2 <- netAnalysis_signalingRole_heatmap(cellchat.floating)
    png(file.path(outdir, "signalingrole_heatmap_floating.png"), width = 2500, height = 2000, res = 300)
    ComplexHeatmap::draw(ht2)
    dev.off()
    
    # Chord plots for selected pathways (if present)
    pathways_to_plot <- c("LAMININ", "CEACAM", "CLDN", "NOTCH")
    for (pw in pathways_to_plot) {
      if (pw %in% names(cellchat.adherent@netP$pathways)) {
        png(file.path(outdir, paste0("chord_gene_adherent_", pw, ".png")), width = 2000, height = 1600)
        netVisual_chord_gene(cellchat.adherent, signaling = pw, slot.name = "net")
        dev.off()
      }
      if (pw %in% names(cellchat.floating@netP$pathways)) {
        png(file.path(outdir, paste0("chord_gene_floating_", pw, ".png")), width = 2000, height = 1600)
        netVisual_chord_gene(cellchat.floating, signaling = pw, slot.name = "net")
        dev.off()
      }
    }
    
    # Communication patterns (river/river plots)
    tryCatch({
      cellchat.adherent <- identifyCommunicationPatterns(cellchat.adherent, pattern = "outgoing", k = 3)
      cellchat.adherent <- identifyCommunicationPatterns(cellchat.adherent, pattern = "incoming", k = 3)
      cellchat.floating <- identifyCommunicationPatterns(cellchat.floating, pattern = "outgoing", k = 3)
      cellchat.floating <- identifyCommunicationPatterns(cellchat.floating, pattern = "incoming", k = 3)
      
      png(file.path(outdir, "river_incoming_adherent.png"), width = 1200, height = 900)
      netAnalysis_river(cellchat.adherent, pattern = "incoming")
      dev.off()
      
      png(file.path(outdir, "river_incoming_floating.png"), width = 1200, height = 900)
      netAnalysis_river(cellchat.floating, pattern = "incoming")
      dev.off()
    }, error = function(e) { log_msg("Communication pattern analysis failed: %s", e$message) })
    
  }, error = function(e) {
    log_msg("CellChat analysis failed: %s", e$message)
  })
} else {
  log_msg("CellChat not installed; skip CellChat analysis.")
}

#########################
# 16. Wrap-up: save Seurat object and session info
#########################
log_msg("Saving processed Seurat object and session info...")
saveRDS(combined, file = file.path(outdir, "combined_seurat_processed.rds"))
writeLines(capture.output(sessionInfo()), con = file.path(outdir, "sessionInfo.txt"))

log_msg("All done. Outputs are in %s", normalizePath(outdir))

# End of extended analyses

log_msg("Pipeline finished. All outputs are in: %s", normalizePath(outdir))
