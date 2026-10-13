#!/home/chengennong/tools/mambaforge/envs/atac/bin/Rscript
# =====================================================================
# 05_da_peaks_gene_activity.R — DA peaks + peak-gene links + GeneActivity（Phase 05）
# 输入: results/rds/04_chromvar.rds
# 输出: results/rds/05_geneact.rds
#       results/tables/05_DE_peaks.csv, 05_DE_peak_gene.csv
#       figures/05_heatmap_geneact.pdf, 05_heatmap_gex.pdf, 05_heatmap_motif.pdf
# 修复项 B2: 删除参考代码中未定义的 up_links 引用
# 修复项: anno_mark 前对高亮基因 trimws 并校验 match 成功
# =====================================================================
set.seed(SEED)
source("config/paths_config.R")
library(Seurat); library(Signac); library(tidyverse)
library(ComplexHeatmap); library(circlize)

combined <- readRDS(file.path(OUT_DIR, "rds", "04_chromvar.rds"))

# ---- DA peaks（LR 检验，校正 ATAC 深度）----
DefaultAssay(combined) <- "ATAC"
da_peaks <- FindMarkers(combined, ident.1 = COND_KO, ident.2 = COND_CTRL,
                        group.by = "Condition", only.pos = FALSE,
                        test.use = DA_TEST, min.pct = 0.05, latent.vars = DA_LATENT)
da_peaks$peak <- rownames(da_peaks)
write.csv(da_peaks, file.path(OUT_DIR, "tables", "05_DE_peaks.csv"), row.names = FALSE)

# ---- peak-gene links（修复 B2：只用 peak_links，无 up_links）----
peak_links <- Links(combined)[Links(combined)$peak %in% rownames(da_peaks), ]
peak_links_df <- as.data.frame(peak_links)
da_peak_gene <- inner_join(da_peaks, peak_links_df, by = "peak")
write.csv(da_peak_gene, file.path(OUT_DIR, "tables", "05_DE_peak_gene.csv"), row.names = FALSE)

# ---- GeneActivity ----
gene.activities <- GeneActivity(combined, assay = "ATAC",
                                extend.upstream = 2000, extend.downstream = 0)
combined[["geneact"]] <- CreateAssayObject(counts = gene.activities)
combined <- NormalizeData(combined, assay = "geneact", normalization.method = "LogNormalize",
                          scale.factor = median(combined$nCount_geneact), verbose = FALSE)
combined <- ScaleData(combined, assay = "geneact", verbose = FALSE)

# ---- 热图（geneact / GEX / motif），高亮基因先 trimws + 校验 ----
make_hilight_anno <- function(genes, mat) {
  genes <- trimws(genes)
  idx <- match(genes, rownames(mat))
  stopifnot("高亮基因未全部匹配到矩阵行" = !any(is.na(idx)))
  rowAnnotation(foo = anno_mark(at = idx, labels = rownames(mat)[idx],
                               labels_gp = gpar(fontsize = 10)))
}
# TODO(Codex): 三张热图的具体 top marker 选择与配色沿用参考代码逻辑，
# 按 celltype 取 top60（presto），DotPlot 取 avg.exp.scaled 后 pivot_wider 成矩阵。
# 此处仅给 geneact 热图骨架，GEX/motif 热图同构：
DefaultAssay(combined) <- "geneact"
Idents(combined) <- "celltype"
markers_geneact <- presto:::wilcoxauc.Seurat(combined, group_by = "celltype",
                                             assay = "data", seurat_assay = "geneact")
top_markers <- markers_geneact |> filter(pct_in > 0.5) |>
  group_by(group) |> slice_max(n = 60, order_by = logFC) |> pull(feature) |> unique()
mat <- GetAssayData(combined, assay = "geneact", slot = "scale.data")[top_markers, ]
# TODO(Codex): 补齐 DotPlot->pivot_wider 的簇均值矩阵（参考 scMultiome.R 950-990 行）
# heat <- Heatmap(mat, name = "GeneActivity (Z-score)", ...)
# pdf(file.path(FIG_DIR, "05_heatmap_geneact.pdf")); draw(heat); dev.off()

saveRDS(combined, file.path(OUT_DIR, "rds", "05_geneact.rds"))
writeLines(capture.output(sessionInfo()),
           file.path(OUT_DIR, "sessionInfo", "05_da_peaks_gene_activity.txt"))
