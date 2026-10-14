#!/home/chengennong/tools/mambaforge/envs/atac/bin/Rscript
# Pilot Pando GRNs with a bounded number of cells and target genes.
# With one sample per condition, network contrasts remain exploratory.
source("config/paths_config.R")
set.seed(SEED)
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(Pando)
  library(JASPAR2020)
  library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg38)
})
BiocParallel::register(BiocParallel::SerialParam())
combined <- readRDS(file.path(OUT_DIR, "rds", "06_scored.rds"))
core_genes <- c("LDHA", "HIF1A", "NFKB1", "REL", "JUN", "FOS", "TBX21",
                "KLF2", "LEF1", "TCF7", "TOX2", "BATF", "GZMB", "PRF1")
target_genes <- unique(c(head(SeuratObject::VariableFeatures(combined, assay = "RNA"), N_PANDO_GENES),
                         core_genes))
target_genes <- intersect(target_genes, rownames(combined[["RNA"]]))
writeLines(target_genes, file.path(OUT_DIR, "tables", "07_target_genes.txt"))
pfm <- getMatrixSet(JASPAR2020, opts = list(species = 9606, all_versions = FALSE))

run_pando <- function(obj, label) {
  cache <- file.path(OUT_DIR, "rds", paste0("07_pando_", label, ".rds"))
  if (file.exists(cache)) {
    muo <- readRDS(cache)
  } else {
    if (ncol(obj) > N_PANDO) {
      obj <- obj[, sample(colnames(obj), N_PANDO)]
    }
    message("[", label, "] cells=", ncol(obj), " target genes=", length(target_genes))
    muo <- initiate_grn(obj, peak_assay = "ATAC", rna_assay = "RNA")
    muo <- find_motifs(muo, pfm = pfm, genome = BSgenome.Hsapiens.UCSC.hg38)
    muo <- infer_grn(muo, genes = target_genes, method = "glm", parallel = FALSE)
    muo <- find_modules(muo)
    network_name <- DefaultNetwork(muo)
    coefs <- muo@grn@networks[[network_name]]@coefs
    saveRDS(coefs, file.path(OUT_DIR, "big", paste0("07_Pando_", label, "_all.rds")))
    if ("padj" %in% colnames(coefs)) {
      coefs <- coefs[order(coefs$padj), , drop = FALSE]
    }
    write.csv(head(coefs, 10000),
              file.path(OUT_DIR, "tables", paste0("07_Pando_", label, ".csv")),
              row.names = FALSE)
    modules <- muo@grn@networks[[network_name]]@modules@meta
    write.csv(modules,
              file.path(OUT_DIR, "tables", paste0("07_Pando_modules_", label, ".csv")),
              row.names = FALSE)
    saveRDS(muo, cache)
  }
  error_path <- file.path(OUT_DIR, "tables", paste0("07_plot_error_", label, ".txt"))
  tryCatch({
    muo <- get_network_graph(muo, umap_method = "none")
    p <- plot_network_graph(muo, layout = "fr")
    ggplot2::ggsave(file.path(FIG_DIR, paste0("07_pando_", label, ".pdf")),
                    p, width = 8, height = 7)
    saveRDS(muo, cache)
    if (file.exists(error_path)) file.remove(error_path)
  }, error = function(e) {
    writeLines(conditionMessage(e), error_path)
  })
  muo
}

for (label in c("combined", COND_CTRL, COND_OE)) {
  obj <- if (label == "combined") combined else
    subset(combined, cells = colnames(combined)[combined$Condition == label])
  run_pando(obj, label)
}
writeLines(capture.output(sessionInfo()),
           file.path(OUT_DIR, "sessionInfo", "07_grn_pando.txt"))
