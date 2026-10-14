#!/home/chengennong/tools/mambaforge/envs/atac/bin/Rscript
# ATAC fragment coverage by condition at selected CD8 CAR-T marker loci.
source("config/paths_config.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(patchwork)
  library(ggplot2)
})

combined <- readRDS(file.path(OUT_DIR, "rds", "05_geneact.rds"))
DefaultAssay(combined) <- "ATAC"
stopifnot(length(Fragments(combined[["ATAC"]])) == 2L,
          length(Annotation(combined[["ATAC"]])) > 0L,
          all(c(COND_CTRL, COND_OE) %in% unique(as.character(combined$Condition))))

markers <- data.frame(
  category = c(rep("Effector", 4L), rep("Proliferation", 3L),
               rep("Resting_memory", 3L)),
  gene = c("GZMB", "GZMA", "PRF1", "IFNG", "MKI67", "TOP2A", "PCNA",
           "TCF7", "IL7R", "KLF2"),
  source = "initial_coverage_panel", tf = "", motif_id = "",
  linked_peak = "", link_score = NA_real_, link_pvalue = NA_real_,
  stringsAsFactors = FALSE
)
# Include every gene shown in the phase-03 dot/feature plots, phase-05 RNA
# and gene-activity heatmaps, and labelled phase-09 volcano plots. Motif
# complexes (e.g. FOS::JUN) have no single RNA locus, so only simple TF
# names with an RNA feature are included from the motif heatmap/volcano.
table_dir <- file.path(OUT_DIR, "tables")
marker_sources <- list(
  state_dotplot = unique(read.csv(file.path(table_dir,
                                       "03_state_gene_markers.csv"))$feature),
  featureplot = c("KLF2", "RGS1", "MKI67", "JUN", "NFKB1", "LDHA"),
  RNA_heatmap = read.csv(file.path(table_dir,
                                   "05_heatmap_gex_features.csv"))$feature,
  geneact_heatmap = read.csv(file.path(table_dir,
                                       "05_heatmap_geneact_features.csv"))$feature,
  motif_heatmap_TF = read.csv(file.path(table_dir,
                                        "05_heatmap_motif_features.csv"))$label,
  RNA_volcano_label = read.csv(file.path(table_dir,
                                         "09_volcano_labels_RNA.csv"))$feature_label,
  geneact_volcano_label = read.csv(file.path(table_dir,
                                             "09_volcano_labels_geneact.csv"))$feature_label,
  motif_volcano_TF = read.csv(file.path(table_dir,
                                        "09_volcano_labels_motif.csv"))$feature_label
)
excluded_coverage <- c("AAVS1", "LDHA")
for (source_name in names(marker_sources)) {
  for (gene in unique(as.character(marker_sources[[source_name]]))) {
    if (is.na(gene) || gene %in% excluded_coverage ||
        !(gene %in% rownames(combined[["RNA"]]) ||
          (source_name == "geneact_heatmap" &&
           gene %in% rownames(combined[["geneact"]])))) next
    j <- match(gene, markers$gene)
    if (is.na(j)) {
      markers <- rbind(markers, data.frame(
        category = "Marker", gene = gene, source = source_name,
        tf = "", motif_id = "", linked_peak = "",
        link_score = NA_real_, link_pvalue = NA_real_))
    } else {
      markers$source[j] <- paste(markers$source[j], source_name, sep = ";")
    }
  }
}

# Candidate targets require a motif hit in a peak linked to the gene by
# Signac LinkPeaks. Keep only linked peaks at the gene locus (within 5 kb of
# the first annotated transcript), so the peak appears in the coverage plot.
# These are candidate regulatory links, not experimentally validated targets.
tf_motifs <- c(ELF4 = "MA0641.1", MEF2D = "MA0773.1", ETV6 = "MA0645.1",
               PPARD = "MA1550.1", KLF3 = "MA1516.1", PBX2 = "MA1113.2")
motif_matrix <- GetMotifData(combined, assay = "ATAC", slot = "data")
links <- as.data.frame(Links(combined[["ATAC"]]))
annotation <- Annotation(combined[["ATAC"]])
stopifnot(all(tf_motifs %in% colnames(motif_matrix)))
for (tf_name in names(tf_motifs)) {
  motif_id <- tf_motifs[[tf_name]]
  motif_peaks <- rownames(motif_matrix)[motif_matrix[, motif_id] != 0]
  candidates <- links[links$peak %in% motif_peaks &
                        links$gene %in% rownames(combined[["RNA"]]) &
                        is.finite(links$pvalue) & links$pvalue < 0.05 &
                        is.finite(links$score) & links$score > 0, ]
  gene_annotation <- annotation[match(candidates$gene, annotation$gene_name)]
  peak_ranges <- StringToGRanges(candidates$peak, sep = c("-", "-"))
  proximal <- !is.na(gene_annotation) &
    as.character(GenomicRanges::seqnames(peak_ranges)) ==
      as.character(GenomicRanges::seqnames(gene_annotation)) &
    GenomicRanges::start(peak_ranges) <= GenomicRanges::end(gene_annotation) + 5000L &
    GenomicRanges::end(peak_ranges) >= GenomicRanges::start(gene_annotation) - 5000L
  candidates <- candidates[proximal & !is.na(proximal), ]
  candidates <- candidates[order(-candidates$score, candidates$pvalue), ]
  candidates <- candidates[!duplicated(candidates$gene) &
                             !candidates$gene %in% markers$gene &
                             !candidates$gene %in% excluded_coverage, ]
  candidates <- head(candidates, 2L)
  if (nrow(candidates) < 2L) {
    stop("Fewer than two proximal motif-linked target candidates for ", tf_name)
  }
  markers <- rbind(markers, data.frame(
    category = "Target", gene = candidates$gene,
    source = "motif_hit_and_peak_gene_link", tf = tf_name,
    motif_id = motif_id, linked_peak = candidates$peak,
    link_score = candidates$score, link_pvalue = candidates$pvalue))
}
markers$rna_feature <- markers$gene
markers$rna_feature[!markers$rna_feature %in% rownames(combined[["RNA"]])] <- NA_character_
markers$region <- markers$gene
annotation_genes <- unique(as.character(Annotation(combined[["ATAC"]])$gene_name))
annotation_locus <- markers$gene
missing_annotation <- setdiff(annotation_locus, annotation_genes)
skipped <- markers[FALSE, ]
skipped$reason <- character()
if (length(missing_annotation)) {
  message("Genes absent from ATAC annotation: ",
          paste(missing_annotation, collapse = ", "))
  omitted <- markers[annotation_locus %in% missing_annotation, ]
  omitted$reason <- "Absent from ATAC gene annotation; no unambiguous coverage locus"
  skipped <- rbind(skipped, omitted)
  markers <- markers[!annotation_locus %in% missing_annotation, ]
}
mitochondrial <- grepl("^MT-", markers$gene)
if (any(mitochondrial)) {
  omitted <- markers[mitochondrial, ]
  omitted$reason <- "Mitochondrial locus has no usable ATAC fragment coverage"
  skipped <- rbind(skipped, omitted)
  markers <- markers[!mitochondrial, ]
}
write.csv(skipped, file.path(table_dir, "10_coverage_skipped.csv"),
          row.names = FALSE)
stopifnot(all(markers$gene %in% annotation_genes),
          all(is.na(markers$rna_feature) |
                markers$rna_feature %in% rownames(combined[["RNA"]])),
          !anyDuplicated(markers$gene))
markers$figure <- vapply(seq_len(nrow(markers)), function(i) {
  if (markers$category[i] == "Target") {
    sprintf("10_coverage_Target_%s_%s.pdf", markers$tf[i], markers$gene[i])
  } else sprintf("10_coverage_%s_%s.pdf", markers$category[i], markers$gene[i])
}, character(1))
message("Coverage plan: ", nrow(markers), " unique loci, including ",
        sum(markers$category == "Target"), " TF target candidates")
write.csv(markers, file.path(table_dir, "10_coverage_markers.csv"),
          row.names = FALSE)
if ("--plan-only" %in% commandArgs(trailingOnly = TRUE)) quit(status = 0)

for (i in seq_len(nrow(markers))) {
  gene <- markers$gene[i]
  category <- markers$category[i]
  filename <- markers$figure[i]
  if (file.exists(file.path(FIG_DIR, filename)) &&
      !(category == "Target" &&
        "--refresh-targets" %in% commandArgs(trailingOnly = TRUE))) next
  message("Coverage: ", category, " / ", gene)
  highlight <- if (category == "Target") {
    StringToGRanges(markers$linked_peak[i], sep = c("-", "-"))
  } else NULL
  caption <- if (category == "Target") {
    paste0("Candidate link: motif ", markers$motif_id[i], " at ",
           markers$linked_peak[i], " (shaded). n=1 per condition.")
  } else {
    if (is.na(markers$rna_feature[i])) {
      "No RNA feature for this gene; ATAC coverage only. One sample per condition."
    } else "CTRL and OE each represent one biological sample; tracks are descriptive."
  }
  p <- CoveragePlot(
    object = combined, region = markers$region[i], assay = "ATAC",
    features = if (is.na(markers$rna_feature[i])) NULL else markers$rna_feature[i],
    expression.assay = "RNA",
    group.by = "Condition", annotation = TRUE, peaks = TRUE,
    region.highlight = highlight,
    links = FALSE, window = 100, extend.upstream = 5000,
    extend.downstream = 5000
  ) +
    plot_annotation(
      title = paste(if (category == "Target") {
        paste(markers$tf[i], "→", gene, "candidate target:")
      } else paste(category, gene), if (is.na(markers$rna_feature[i])) {
        "ATAC coverage by condition"
      } else "ATAC coverage and RNA expression by condition"),
      caption = caption
    )
  ggsave(file.path(FIG_DIR, filename), p, width = 11, height = 6.8)
}

figures <- sort(basename(list.files(FIG_DIR, pattern = "\\.pdf$")))
write.csv(data.frame(figure = figures),
          file.path(OUT_DIR, "tables", "figure_index.csv"), row.names = FALSE)
writeLines(sub("[[:space:]]+$", "", capture.output(sessionInfo())),
           file.path(OUT_DIR, "sessionInfo", "10_condition_coverage.txt"))
