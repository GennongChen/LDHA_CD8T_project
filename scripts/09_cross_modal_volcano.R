#!/home/chengennong/tools/mambaforge/envs/atac/bin/Rscript
# Exploratory cell-level comparisons: one biological sample per condition.
source("config/paths_config.R")
set.seed(SEED)
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggrepel)
})

stopifnot(packageVersion("ggplot2") >= "3.5.0")
combined <- readRDS(file.path(OUT_DIR, "rds", "05_geneact.rds"))
stopifnot(all(c(COND_CTRL, COND_OE) %in% unique(as.character(combined$Condition))))

table_dir <- file.path(OUT_DIR, "tables")
top_n <- 5L
motif_cut <- 0.25 # chromVAR deviation-score units, not log2 fold change
palette <- c(Higher_OE = "#A23B72", Higher_CTRL = "#287D74", Other = "grey75")
subtitle <- "Exploratory cell-level test; one biological sample per condition"

get_de <- function(assay, filename) {
  path <- file.path(table_dir, filename)
  if (file.exists(path)) {
    de <- read.csv(path, check.names = FALSE)
  } else {
    message("Testing ", assay, " features")
    DefaultAssay(combined) <- assay
    de <- FindMarkers(combined, group.by = "Condition",
                      ident.1 = COND_OE, ident.2 = COND_CTRL,
                      assay = assay, test.use = "wilcox",
                      min.pct = 0.05, logfc.threshold = 0)
    de$gene <- rownames(de)
    write.csv(de, path, row.names = FALSE)
  }
  stopifnot(all(c("gene", "avg_log2FC", "p_val_adj") %in% names(de)),
            !anyDuplicated(de$gene))
  de
}

rna <- get_de("RNA", "09_RNA_DE_full.csv")
geneact <- get_de("geneact", "09_geneact_DE.csv")
motif <- read.csv(file.path(table_dir, "04_motifs_by_condition.csv"),
                  check.names = FALSE)
motif <- motif[motif$group == COND_OE, ]
stopifnot(nrow(motif) > 0L, !anyDuplicated(motif$feature),
          all(motif$feature %in% rownames(combined[["chromvar"]])))
# Use the actual signed difference in mean deviation Z score. ChromVAR scores
# can be negative, so a ratio or log2 fold change is undefined.
z <- GetAssayData(combined, assay = "chromvar", layer = "data")
oe <- which(combined$Condition == COND_OE)
ctrl <- which(combined$Condition == COND_CTRL)
delta <- Matrix::rowMeans(z[, oe, drop = FALSE]) -
         Matrix::rowMeans(z[, ctrl, drop = FALSE])
motif$delta_deviation <- unname(delta[motif$feature])
stopifnot(all(is.finite(motif$delta_deviation)))
write.csv(motif, file.path(table_dir, "09_motif_with_delta.csv"), row.names = FALSE)

volcano <- function(df, feature_col, effect_col, padj_col, label, x_title,
                    cut, filename, forced = character()) {
  df$feature_label <- as.character(df[[feature_col]])
  missing_label <- is.na(df$feature_label) | df$feature_label == ""
  if (any(missing_label)) df$feature_label[missing_label] <- df$feature[missing_label]
  df$effect <- df[[effect_col]]
  df$padj <- df[[padj_col]]
  df <- df[is.finite(df$effect) & is.finite(df$padj), ]
  df$neglog10_padj <- pmin(-log10(pmax(df$padj, .Machine$double.xmin)), 50)
  df$direction <- ifelse(df$padj < PADJ_CUT & df$effect >= cut, "Higher_OE",
                         ifelse(df$padj < PADJ_CUT & df$effect <= -cut,
                                "Higher_CTRL", "Other"))
  df$direction <- factor(df$direction, levels = names(palette))
  candidates <- df[df$direction != "Other", ]
  candidates <- candidates[order(candidates$padj, -abs(candidates$effect)), ]
  candidates <- candidates[!duplicated(candidates$feature_label), ]
  labels <- do.call(rbind, lapply(c("Higher_OE", "Higher_CTRL"), function(side) {
    head(candidates[as.character(candidates$direction) == side, ], top_n)
  }))
  labels <- rbind(labels, df[df$feature_label %in% forced, ])
  labels <- labels[!duplicated(labels$feature_label), ]
  write.csv(labels[, c("feature_label", "effect", "padj", "direction")],
            file.path(table_dir, paste0("09_volcano_labels_", label, ".csv")),
            row.names = FALSE)
  p <- ggplot(df, aes(effect, neglog10_padj)) +
    geom_hline(yintercept = -log10(PADJ_CUT), linetype = 2, color = "grey55") +
    geom_vline(xintercept = c(-cut, cut), linetype = 2, color = "grey55") +
    geom_point(aes(color = direction), size = 0.8, alpha = 0.7) +
    geom_text_repel(data = labels, aes(label = feature_label),
                    size = 2.8, max.overlaps = Inf, seed = SEED,
                    box.padding = 0.35, point.padding = 0.2,
                    segment.color = "grey55") +
    scale_y_continuous(limits = c(0, 62), breaks = seq(0, 50, 10),
                       expand = expansion(mult = c(0.01, 0.01))) +
    scale_color_manual(values = palette, drop = FALSE,
                       labels = c(Higher_OE = "Higher in OE",
                                  Higher_CTRL = "Higher in CTRL", Other = "Other")) +
    labs(title = paste(label, "volcano: OE versus CTRL"), subtitle = subtitle,
         x = x_title, y = "-log10 adjusted cell-level p value (capped at 50)",
         color = NULL, caption = paste0("Labels: top ", top_n,
                                       " per direction by adjusted p value; dashed cutoffs shown.")) +
    theme_classic(base_size = 11) +
    theme(legend.position = "bottom", plot.caption = element_text(size = 8))
  ggsave(file.path(FIG_DIR, filename), p, width = 8, height = 6.5)
  message(filename, ": ", nrow(df), " tested; ", nrow(labels), " labels")
}

volcano(rna, "gene", "avg_log2FC", "p_val_adj", "RNA",
        "RNA log2 fold change (OE / CTRL)", LOGFC_CUT,
        "09_volcano_RNA.pdf", forced = "LDHA")
volcano(geneact, "gene", "avg_log2FC", "p_val_adj", "geneact",
        "Gene activity log2 fold change (OE / CTRL)", LOGFC_CUT,
        "09_volcano_geneact.pdf", forced = "LDHA")
volcano(motif, "gene", "delta_deviation", "padj", "motif",
        "Mean chromVAR deviation difference (OE - CTRL)", motif_cut,
        "09_volcano_motif.pdf")

# Multiple JASPAR motifs can map to the same TF. Select the motif with the
# smallest adjusted p value, then the largest absolute deviation difference.
motif <- motif[!is.na(motif$gene) & motif$gene != "", ]
motif <- motif[order(motif$padj, -abs(motif$delta_deviation)), ]
motif <- motif[!duplicated(motif$gene), ]
shared <- merge(rna[, c("gene", "avg_log2FC", "p_val_adj")],
                motif[, c("gene", "feature", "delta_deviation", "padj")],
                by = "gene")
names(shared) <- c("gene", "rna_log2FC", "rna_padj", "motif_id",
                   "motif_delta_deviation", "motif_padj")
shared <- shared[is.finite(shared$rna_log2FC) &
                   is.finite(shared$motif_delta_deviation), ]
stopifnot(nrow(shared) > 0L, !anyDuplicated(shared$gene))
shared$association <- ifelse(shared$rna_padj < PADJ_CUT &
                               shared$motif_padj < PADJ_CUT,
                             ifelse(sign(shared$rna_log2FC) ==
                                      sign(shared$motif_delta_deviation),
                                    "Concordant", "Discordant"), "Other")
shared$labelled <- shared$association != "Other"
write.csv(shared,
          file.path(table_dir, "09_RNA_x_motif.csv"), row.names = FALSE)
p <- ggplot(shared, aes(rna_log2FC, motif_delta_deviation)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey60") +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  geom_point(aes(color = association), size = 1.6, alpha = 0.8) +
  geom_text_repel(data = shared[shared$labelled, ], aes(label = gene),
                  size = 3, max.overlaps = Inf, seed = SEED,
                  box.padding = 0.35, segment.color = "grey55") +
  scale_color_manual(values = c(Concordant = "#A23B72",
                                Discordant = "#287D74", Other = "grey75")) +
  labs(title = "RNA expression versus TF motif activity",
       subtitle = subtitle,
       x = "TF RNA log2 fold change (OE / CTRL)",
       y = "Mean chromVAR deviation difference (OE - CTRL)",
       color = "Direction",
       caption = paste0("One JASPAR motif per TF (lowest adjusted p value); ",
                        "colored points have both adjusted p values < ", PADJ_CUT,
                        ". All concordant and discordant TFs are labeled.")) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom", plot.caption = element_text(size = 8))
ggsave(file.path(FIG_DIR, "09_RNA_x_motif_logFC.pdf"), p,
       width = 11, height = 8)
message("RNA x motif: ", nrow(shared), " shared TFs; ",
        sum(shared$labelled), " labels")

figures <- sort(basename(list.files(FIG_DIR, pattern = "\\.pdf$")))
write.csv(data.frame(figure = figures), file.path(table_dir, "figure_index.csv"),
          row.names = FALSE)
writeLines(sub("[[:space:]]+$", "", capture.output(sessionInfo())),
           file.path(OUT_DIR, "sessionInfo", "09_cross_modal_volcano.txt"))
