#!/usr/bin/env Rscript
# Reads separately delivered RNA and ATAC matrices, checks paired barcodes,
# filters cells, then requantifies both fragment files on one common peak set.
source("config/paths_config.R")
set.seed(SEED)
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(EnsDb.Hsapiens.v86)
})

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
GenomeInfoDb::seqlevels(annotation) <- sub("^chrMT$", "chrM",
                                           paste0("chr", GenomeInfoDb::seqlevels(annotation)))

read_sample <- function(label, info) {
  message("Reading ", label)
  rna <- Read10X(info$rna, gene.column = 2)
  atac <- Read10X(info$atac, gene.column = 2)
  stopifnot(identical(colnames(rna), colnames(atac)))
  keep <- grepl("^chr([1-9]|1[0-9]|2[0-2]|X|Y):", rownames(atac))
  atac <- atac[keep, , drop = FALSE]
  obj <- CreateSeuratObject(rna, assay = "RNA", project = label)
  obj[["ATAC"]] <- CreateChromatinAssay(
    counts = atac, sep = c(":", "-"), fragments = info$fragments,
    annotation = annotation, validate.fragments = FALSE
  )
  obj$sample_id <- label
  obj$Condition <- info$condition
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, assay = "RNA", pattern = "^MT-")
  obj
}
obj_list <- Map(read_sample, names(samples), samples)
names(obj_list) <- names(samples)

qc_raw <- do.call(rbind, lapply(names(obj_list), function(s) {
  m <- obj_list[[s]][[]]
  data.frame(sample = s, nCount_RNA = m$nCount_RNA,
             nCount_ATAC = m$nCount_ATAC, percent.mt = m$percent.mt)
}))
write.csv(qc_raw, file.path(OUT_DIR, "tables", "01_qc_metrics_before.csv"), row.names = FALSE)
pdf(file.path(FIG_DIR, "01_qc_violin.pdf"), width = 12, height = 5)
print(VlnPlot(obj_list[[1]], features = c("nCount_RNA", "nCount_ATAC", "percent.mt"),
              pt.size = 0, ncol = 3) + ggplot2::ggtitle(names(obj_list)[1]))
print(VlnPlot(obj_list[[2]], features = c("nCount_RNA", "nCount_ATAC", "percent.mt"),
              pt.size = 0, ncol = 3) + ggplot2::ggtitle(names(obj_list)[2]))
dev.off()

qc_summary <- do.call(rbind, lapply(names(obj_list), function(s) {
  obj <- obj_list[[s]]
  m <- obj[[]]
  keep <- m$nCount_ATAC >= QC$nCount_ATAC_min &
    m$nCount_ATAC <= QC$nCount_ATAC_max &
    m$nCount_RNA >= QC$nCount_RNA_min &
    m$nCount_RNA <= QC$nCount_RNA_max &
    m$percent.mt <= QC$percent_mt_max
  obj_list[[s]] <<- subset(obj, cells = rownames(m)[keep])
  data.frame(sample = s, cells_before = ncol(obj),
             cells_after = sum(keep), retention = mean(keep))
}))
write.csv(qc_summary, file.path(OUT_DIR, "tables", "01_qc_summary.csv"), row.names = FALSE)
stopifnot(all(qc_summary$cells_after > 100))

# The two original peak lists have different coordinates. Using a union matrix
# would make condition-specific absent peaks look like zero accessibility.
peaks <- GenomicRanges::reduce(Reduce(c, lapply(obj_list, function(x) granges(x[["ATAC"]]))))
peaks <- peaks[width(peaks) > 20 & width(peaks) < 10000]
message("Common peaks: ", length(peaks))
for (s in names(obj_list)) {
  obj <- obj_list[[s]]
  cache <- file.path(OUT_DIR, "rds", paste0("01_common_peaks_", s, ".rds"))
  key_file <- paste0(cache, ".key")
  key <- digest::digest(list(peaks = as.character(peaks), cells = colnames(obj)))
  use_cache <- file.exists(cache) && file.exists(key_file) &&
    identical(readLines(key_file, warn = FALSE), key)
  if (use_cache) {
    counts <- readRDS(cache)
    use_cache <- nrow(counts) == length(peaks) &&
      identical(colnames(counts), colnames(obj))
  }
  if (!use_cache) {
    counts <- FeatureMatrix(fragments = Fragments(obj[["ATAC"]]),
                            features = peaks, cells = colnames(obj), verbose = TRUE)
    saveRDS(counts, cache)
    writeLines(key, key_file)
  }
  stopifnot(identical(colnames(counts), colnames(obj)))
  obj[["ATAC"]] <- CreateChromatinAssay(
    counts = counts, ranges = peaks, fragments = Fragments(obj[["ATAC"]]),
    annotation = annotation, min.cells = 0,
    min.features = 0
  )
  obj_list[[s]] <- obj
}

combined <- merge(obj_list[[1]], y = obj_list[-1],
                  add.cell.ids = names(obj_list), merge.data = FALSE)
combined <- JoinLayers(combined, assay = "RNA")
combined <- subset(combined, cells = colnames(combined)[
  combined$nCount_ATAC >= QC$nCount_ATAC_min &
    combined$nCount_ATAC <= QC$nCount_ATAC_max
])
qc_summary$cells_after_common_peaks <- as.integer(table(
  factor(combined$sample_id, levels = qc_summary$sample)))
write.csv(qc_summary, file.path(OUT_DIR, "tables", "01_qc_summary.csv"), row.names = FALSE)
stopifnot(length(Layers(combined[["RNA"]], search = "counts")) == 1)
stopifnot(ncol(combined) == sum(qc_summary$cells_after_common_peaks))
saveRDS(combined, file.path(OUT_DIR, "rds", "01_combined_qc.rds"))
writeLines(capture.output(sessionInfo()),
           file.path(OUT_DIR, "sessionInfo", "01_preprocess_qc.txt"))
message("QC complete: ", ncol(combined), " retained cells")
