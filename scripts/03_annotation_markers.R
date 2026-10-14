#!/home/chengennong/tools/mambaforge/envs/atac/bin/Rscript
# =====================================================================
# 03_annotation_markers.R — 簇 marker 与细胞类型注释（CODEX_PLAN Phase 03）
# 输入: results/rds/02_wnn.rds
# 输出: results/rds/03_annotated.rds
#       results/tables/03_markers_clusters.csv
#       figures/03_dotplot_markers.pdf, 03_proportion.pdf, 03_featureplot.pdf
# 修复项: 簇→细胞类型映射从 config$celltype_map 读取（禁止硬编码 levels）
# =====================================================================
source("config/paths_config.R")
set.seed(SEED)
library(Seurat); library(Signac); library(tidyverse)

combined <- readRDS(file.path(OUT_DIR, "rds", "02_wnn.rds"))
DefaultAssay(combined) <- "RNA"
cluster_col <- paste0("wsnn_res", RES_WNN)
Idents(combined) <- cluster_col

# ---- marker（presto）----
markers <- presto::wilcoxauc(combined, group_by = cluster_col,
                             assay = "data", seurat_assay = "RNA")
saveRDS(markers, file.path(OUT_DIR, "big", "03_markers_all.rds"))
top_markers <- markers |>
  dplyr::filter(padj < 0.05, logFC > 0) |>
  dplyr::group_by(group) |>
  dplyr::slice_max(order_by = auc, n = 100, with_ties = FALSE) |>
  dplyr::ungroup()
write.csv(top_markers, file.path(OUT_DIR, "tables", "03_markers_clusters.csv"), row.names = FALSE)

# ---- CD8 T 状态基因 dotplot ----
state_genes <- c("CD8A","KLF2","LEF1","TCF7","IL7R",          # naive/stem
                 "GZMB","PRF1","NKG7","TBX21","IFNG",          # effector
                 "MKI67","TOP2A","TUBB",                       # cycling
                 "CXCL13","TOX2","ENTPD1",                     # exhausted
                 "NR4A1","NR4A2","HIF1A","LDHA","HK2","SLC16A3") # activation/糖酵解
state_genes <- intersect(state_genes, rownames(combined))
state_stats <- markers |>
  dplyr::filter(feature %in% state_genes) |>
  dplyr::select(group, feature, logFC, pct_in, auc, padj)
write.csv(state_stats, file.path(OUT_DIR, "tables", "03_state_gene_markers.csv"), row.names = FALSE)
cluster_qc <- do.call(rbind, lapply(split(combined[[]], combined[[cluster_col]][, 1]),
                                   function(m) data.frame(
                                     cluster = as.character(m[[cluster_col]][1]),
                                     cells = nrow(m),
                                     fraction_OE = mean(m$Condition == COND_OE),
                                     median_RNA = median(m$nCount_RNA),
                                     median_ATAC = median(m$nCount_ATAC),
                                     median_percent_mt = median(m$percent.mt))))
write.csv(cluster_qc, file.path(OUT_DIR, "tables", "03_cluster_qc.csv"), row.names = FALSE)
p_dot <- DotPlot(combined, group.by = cluster_col, features = state_genes,
                 cols = c("lightgrey", "#CD534CFF"), dot.scale = 6) + RotatedAxis()
ggsave(file.path(FIG_DIR, "03_dotplot_markers.pdf"), p_dot, width = 12, height = 5)

# ---- 注释：映射表来自 config，人工确认后填写 ----
if (length(celltype_map) == 0) stop("请先看 marker 图，在 config/paths_config.R 填写 celltype_map 后重跑")
stopifnot(setequal(names(celltype_map), levels(Idents(combined))))
combined$celltype <- factor(as.character(combined@meta.data[[cluster_col]]),
                            levels = names(celltype_map), labels = celltype_map)
stopifnot(!any(is.na(combined$celltype)))  # 每个簇都必须有映射

# ---- 比例图 ----
pal <- scales::hue_pal()(nlevels(combined$celltype))
p_prop <- ggplot(combined@meta.data, aes(x = Condition, fill = celltype)) +
  geom_bar(position = "fill", width = 0.7, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = pal) + coord_flip() +
  labs(x = "Condition", y = "Proportion") + theme_classic()
ggsave(file.path(FIG_DIR, "03_proportion.pdf"), p_prop, width = 4, height = 3)

# ---- FeaturePlot（含 LDHA 本身）----
p_fp <- FeaturePlot(combined, features = intersect(c("KLF2","RGS1","MKI67","JUN","NFKB1","LDHA"),
                                                   rownames(combined)),
                    reduction = "wnn.umap", ncol = 3, cols = c("lightgrey", "#CD534CFF"))
ggsave(file.path(FIG_DIR, "03_featureplot.pdf"), p_fp, width = 11, height = 7)

saveRDS(combined, file.path(OUT_DIR, "rds", "03_annotated.rds"))
writeLines(capture.output(sessionInfo()),
           file.path(OUT_DIR, "sessionInfo", "03_annotation_markers.txt"))
