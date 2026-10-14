#!/usr/bin/env Rscript
# Exploratory ATAC differential accessibility and peak-to-gene associations.
# There is one biological sample per condition; cell-level p values are not
# evidence of a reproducible LDHA effect.
source("config/paths_config.R")
set.seed(SEED)
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(ComplexHeatmap)
})
final_cache <- file.path(OUT_DIR, "rds", "05_geneact.rds")
if (file.exists(final_cache)) {
  combined <- readRDS(final_cache)
} else {
combined <- readRDS(file.path(OUT_DIR, "rds", "04_chromvar_prefootprint.rds"))

DefaultAssay(combined) <- "ATAC"
da_path <- file.path(OUT_DIR, "tables", "05_DE_peaks.csv")
if (file.exists(da_path)) {
  da_peaks <- read.csv(da_path)
  stopifnot(nrow(da_peaks) > 0, all(c("peak", "p_val_adj", "avg_log2FC") %in% names(da_peaks)))
} else {
  da_peaks <- FindMarkers(
    combined, ident.1 = COND_OE, ident.2 = COND_CTRL,
    group.by = "Condition", only.pos = FALSE, test.use = DA_TEST,
    min.pct = 0.05, latent.vars = DA_LATENT
  )
  da_peaks$peak <- rownames(da_peaks)
  write.csv(da_peaks, da_path, row.names = FALSE)
}

# Restrict links to RNA variable genes to keep computation bounded and
# explicitly record the tested gene set.
genes_to_link <- intersect(VariableFeatures(combined[["RNA"]]), rownames(combined[["RNA"]]))
writeLines(genes_to_link, file.path(OUT_DIR, "tables", "05_linked_genes.txt"))
links_cache <- file.path(OUT_DIR, "rds", "05_linked.rds")
if (file.exists(links_cache)) {
  combined <- readRDS(links_cache)
} else {
  combined <- LinkPeaks(combined, peak.assay = "ATAC", expression.assay = "RNA",
                        genes.use = genes_to_link, n_sample = 200)
  saveRDS(combined, links_cache)
}
peak_links <- as.data.frame(Links(combined[["ATAC"]]))
stopifnot(all(c("peak", "gene") %in% colnames(peak_links)))
da_peak_gene <- dplyr::inner_join(da_peaks, peak_links, by = "peak")
write.csv(da_peak_gene, file.path(OUT_DIR, "tables", "05_DE_peak_gene.csv"), row.names = FALSE)

gene.activities <- GeneActivity(combined, assay = "ATAC",
                                extend.upstream = 2000, extend.downstream = 0)
combined[["geneact"]] <- CreateAssayObject(counts = gene.activities)
combined <- NormalizeData(combined, assay = "geneact", verbose = FALSE)
saveRDS(combined, final_cache)
}

make_heatmap <- function(assay, label, top_n = 15L) {
  DefaultAssay(combined) <- assay
  markers <- presto::wilcoxauc(combined, group_by = "celltype",
                               assay = "data", seurat_assay = assay)
  top <- markers |>
    dplyr::filter(padj < 0.05, logFC > 0, pct_in > 10) |>
    dplyr::group_by(group) |>
    dplyr::slice_max(order_by = auc, n = top_n, with_ties = FALSE) |>
    dplyr::pull(feature) |>
    unique()
  if (length(top) < 2) {
    top <- markers |>
      dplyr::group_by(group) |>
      dplyr::slice_max(order_by = auc, n = top_n, with_ties = FALSE) |>
      dplyr::pull(feature) |>
      unique()
  }
  top <- intersect(top, rownames(combined[[assay]]))
  if (assay == "chromvar") {
    # chromVAR data are deviation Z scores, not log expression. Average
    # them directly; AverageExpression() would incorrectly exponentiate them.
    values <- GetAssayData(combined, assay = assay, layer = "data")
    groups <- factor(combined$celltype)
    mat <- do.call(cbind, lapply(levels(groups), function(g) {
      Matrix::rowMeans(values[top, which(groups == g), drop = FALSE])
    }))
    rownames(mat) <- top
    colnames(mat) <- levels(groups)
  } else {
    avg <- AverageExpression(combined, assays = assay, features = top,
                             group.by = "celltype", layer = "data",
                             verbose = FALSE)[[assay]]
    mat <- as.matrix(avg)
  }
  mat <- t(scale(t(mat)))
  mat[!is.finite(mat)] <- 0
  labels <- if (assay == "chromvar") {
    make.unique(ConvertMotifID(combined, id = rownames(mat), assay = "ATAC"))
  } else rownames(mat)
  feature_table <- data.frame(feature = rownames(mat), label = labels)
  rownames(mat) <- labels
  pdf(file.path(FIG_DIR, paste0("05_heatmap_", label, ".pdf")),
      width = 12, height = max(7, min(20, nrow(mat) * 0.28)))
  draw(Heatmap(mat, name = "row Z", show_row_names = nrow(mat) <= 80,
               row_names_gp = grid::gpar(fontsize = 8),
               column_title = paste(assay, "average by CD8 state")))
  dev.off()
  write.csv(feature_table,
            file.path(OUT_DIR, "tables", paste0("05_heatmap_", label, "_features.csv")),
            row.names = FALSE)
}
make_heatmap("geneact", "geneact", top_n = 8L)
make_heatmap("RNA", "gex", top_n = 8L)
make_heatmap("chromvar", "motif", top_n = 5L)

writeLines(capture.output(sessionInfo()),
           file.path(OUT_DIR, "sessionInfo", "05_da_peaks_gene_activity.txt"))
