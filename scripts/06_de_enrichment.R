#!/usr/bin/env Rscript
# One biological sample per condition: every cell-level test below is
# exploratory. No pseudobulk hypothesis test is valid with n=1 per group.
source("config/paths_config.R")
set.seed(SEED)
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})
# RNA differential expression and module scores do not depend on the
# computationally expensive phase-05 peak-to-gene links.
combined <- readRDS(file.path(OUT_DIR, "rds", "04_chromvar_prefootprint.rds"))
DefaultAssay(combined) <- "RNA"
write.csv(data.frame(condition = c(COND_CTRL, COND_OE),
                     biological_samples = c(1L, 1L),
                     inference = "exploratory cell-level comparison"),
          file.path(OUT_DIR, "tables", "06_analysis_design.csv"), row.names = FALSE)

de_presto <- presto::wilcoxauc(combined, group_by = "Condition",
                               assay = "data", seurat_assay = "RNA")
saveRDS(de_presto, file.path(OUT_DIR, "big", "06_DE_presto_all.rds"))
write.csv(subset(de_presto, group == COND_OE & padj < 0.05),
          file.path(OUT_DIR, "tables", "06_DE_presto_OE.csv"), row.names = FALSE)
de <- FindMarkers(combined, group.by = "Condition",
                  ident.1 = COND_OE, ident.2 = COND_CTRL,
                  test.use = "wilcox", min.pct = 0.05, logfc.threshold = 0.1)
de$gene <- rownames(de)
de$direction <- dplyr::case_when(
  de$avg_log2FC >= LOGFC_CUT & de$p_val_adj < PADJ_CUT ~ "Higher_OE",
  de$avg_log2FC <= -LOGFC_CUT & de$p_val_adj < PADJ_CUT ~ "Higher_CTRL",
  TRUE ~ "Other")
write.csv(de, file.path(OUT_DIR, "tables", "06_DE_condition.csv"), row.names = FALSE)
p <- ggplot(de, aes(avg_log2FC, pmin(-log10(pmax(p_val_adj, .Machine$double.xmin)), 50),
                    color = direction)) +
  geom_point(size = 0.6, alpha = 0.7) +
  scale_color_manual(values = c(Higher_OE = "#984EA3",
                                Higher_CTRL = "#4DAF4A", Other = "grey70")) +
  labs(x = "log2 fold change (OE / CTRL)",
       y = "-log10 adjusted cell-level p value (capped at 50)",
       title = "Exploratory: one biological sample per condition") +
  theme_classic()
ggsave(file.path(FIG_DIR, "06_volcano.pdf"), p, width = 6, height = 5)

for (side in c("Higher_OE", "Higher_CTRL")) {
  genes <- de$gene[de$direction == side]
  go <- if (length(genes) >= 10) {
    enrichGO(gene = genes, universe = rownames(combined[["RNA"]]),
             OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
             pAdjustMethod = "BH", readable = TRUE)
  } else NULL
  go_df <- if (is.null(go)) data.frame() else as.data.frame(go)
  write.csv(go_df, file.path(OUT_DIR, "tables", paste0("06_GO_", side, ".csv")),
            row.names = FALSE)
  if (nrow(go_df) > 0) {
    ggsave(file.path(FIG_DIR, paste0("06_GO_", side, ".pdf")),
           enrichplot::dotplot(go, showCategory = 15), width = 8, height = 6)
  }
  ids <- if (length(genes) > 0) {
    suppressMessages(bitr(genes, fromType = "SYMBOL",
                          toType = "ENTREZID", OrgDb = org.Hs.eg.db))
  } else data.frame(ENTREZID = character())
  if (is.null(ids)) ids <- data.frame(ENTREZID = character())
  kegg <- if (nrow(ids) >= 10) {
    options(timeout = 20)
    tryCatch(enrichKEGG(gene = ids$ENTREZID, organism = "hsa",
                        pAdjustMethod = "BH"), error = function(e) {
                          message("KEGG unavailable for ", side, ": ", conditionMessage(e))
                          NULL
                        })
  } else NULL
  kegg_df <- if (is.null(kegg)) data.frame() else as.data.frame(kegg)
  write.csv(kegg_df, file.path(OUT_DIR, "tables", paste0("06_KEGG_", side, ".csv")),
            row.names = FALSE)
  if (nrow(kegg_df) > 0) {
    ggsave(file.path(FIG_DIR, paste0("06_KEGG_", side, ".pdf")),
           enrichplot::dotplot(kegg, showCategory = 15), width = 8, height = 6)
  }
}

modules <- read.delim(MODULES_FILE, stringsAsFactors = FALSE)
module_scores <- list()
for (i in seq_len(nrow(modules))) {
  module <- modules$module[i]
  genes <- trimws(strsplit(modules$genes[i], ",", fixed = TRUE)[[1]])
  genes <- intersect(genes, rownames(combined[["RNA"]]))
  if (length(genes) < 5) next
  prefix <- paste0("module_", module)
  combined <- AddModuleScore(combined, features = list(genes), name = prefix,
                             assay = "RNA", seed = SEED)
  col <- paste0(prefix, "1")
  module_scores[[module]] <- data.frame(Condition = combined$Condition,
                                       module = module, score = combined[[col]][, 1])
}
scores <- do.call(rbind, module_scores)
write.csv(scores, file.path(OUT_DIR, "tables", "06_module_scores.csv"),
          row.names = FALSE)
module_summary <- scores |>
  dplyr::group_by(Condition, module) |>
  dplyr::summarise(median_score = median(score), cells = dplyr::n(), .groups = "drop")
write.csv(module_summary, file.path(OUT_DIR, "tables", "06_module_summary.csv"),
          row.names = FALSE)
p_modules <- ggplot(scores, aes(Condition, score, fill = Condition)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
  facet_wrap(~ module, scales = "free_y", ncol = 3) +
  labs(title = "Descriptive module scores; one sample per condition") +
  theme_classic() + theme(legend.position = "none")
ggsave(file.path(FIG_DIR, "06_module_violin.pdf"), p_modules,
       width = 10, height = 7)

saveRDS(combined, file.path(OUT_DIR, "rds", "06_scored.rds"))
writeLines(capture.output(sessionInfo()),
           file.path(OUT_DIR, "sessionInfo", "06_de_enrichment.txt"))
