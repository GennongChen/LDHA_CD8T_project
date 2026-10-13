#!/home/chengennong/tools/mambaforge/envs/atac/bin/Rscript
# =====================================================================
# 04_chromvar_footprint.R — chromVAR motif 活性 + footprint（Phase 04）
# 输入: results/rds/03_annotated.rds
# 输出: results/rds/04_chromvar.rds
#       results/tables/04_motifs_by_celltype.csv, 04_motifs_by_condition.csv
#       figures/04_footprint_condition.pdf, 04_footprint_celltype.pdf
# 修复项 B1: 落盘的 motif 差异表必须写 markers_motifs 对象（参考代码写错成 RNA）
# 修复项: topTFs 按列名取子集（参考代码按列号，presto 版本变化即失效）
# =====================================================================
set.seed(SEED)
source("config/paths_config.R")
library(Seurat); library(Signac); library(tidyverse)
library(chromVAR); library(JASPAR2020); library(TFBSTools); library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)

combined <- readRDS(file.path(OUT_DIR, "rds", "03_annotated.rds"))

# ---- chromVAR ----
DefaultAssay(combined) <- "ATAC"
pwm_set <- getMatrixSet(JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
motif.matrix <- CreateMotifMatrix(features = granges(combined), pwm = pwm_set,
                                  genome = BSgenome.Hsapiens.UCSC.hg38, use.counts = FALSE)
combined <- SetAssayData(combined, assay = "ATAC", slot = "motifs",
                         new.data = CreateMotifObject(data = motif.matrix, pwm = pwm_set))
combined <- RunChromVAR(combined, genome = BSgenome.Hsapiens.UCSC.hg38)

# ---- 差异 motif（按 celltype / Condition）----
motif_by_ct <- presto:::wilcoxauc.Seurat(combined, group_by = "celltype",
                                         assay = "data", seurat_assay = "chromvar")
motif_by_cd <- presto:::wilcoxauc.Seurat(combined, group_by = "Condition",
                                         assay = "data", seurat_assay = "chromvar")
# 修复 B1：写盘对象用 motif 结果
write.csv(motif_by_ct, file.path(OUT_DIR, "tables", "04_motifs_by_celltype.csv"), row.names = FALSE)
write.csv(motif_by_cd %>% filter(group == COND_KO),
          file.path(OUT_DIR, "tables", "04_motifs_by_condition.csv"), row.names = FALSE)

# ---- topTFs：RNA + motif 联合（按列名取子集）----
rna_by_ct <- presto:::wilcoxauc.Seurat(combined, group_by = "celltype",
                                       assay = "data", seurat_assay = "RNA")
motif_by_ct$gene <- ConvertMotifID(combined, id = motif_by_ct$feature, assay = "ATAC")

topTFs <- function(celltype, markers_rna, markers_motif, padj.cutoff = 0.05) {
  r <- dplyr::filter(markers_rna, group == celltype, padj < padj.cutoff, logFC > 0) |>
    select(group, feature, logFC, auc) |> arrange(-auc)
  m <- dplyr::filter(markers_motif, group == celltype, padj < padj.cutoff, logFC > 0) |>
    select(group, gene, logFC, auc) |> arrange(-auc)
  inner_join(r, m, by = c("feature" = "gene"), suffix = c(".rna", ".motif")) |>
    mutate(avg_auc = (auc.rna + auc.motif) / 2) |> arrange(-avg_auc)
}
for (ct in levels(combined$celltype)) {
  write.csv(topTFs(ct, rna_by_ct, motif_by_ct),
            file.path(OUT_DIR, "tables", sprintf("04_topTFs_%s.csv", ct)), row.names = FALSE)
}

# ---- footprint（内存峰值高，单独跑）----
# TODO(Codex): TF 列表按 LDHA 项目调整；建议含 HIF1A/ARNTL + NFKB/AP-1/KLF2 轴
fp_tfs <- c("NFKB1", "JUN", "KLF2", "LEF1", "HIF1A", "TBX21")
fp_tfs <- intersect(fp_tfs, rownames(combined@assays$ATAC@motifs))
combined <- Footprint(combined, motif.name = fp_tfs,
                      genome = BSgenome.Hsapiens.UCSC.hg38)
p1 <- PlotFootprint(combined, features = fp_tfs, group.by = "Condition", label = FALSE) &
  scale_color_manual(values = setNames(c("#B8DA68", "#BF8AC0"), c(COND_CTRL, COND_KO)))
ggsave(file.path(FIG_DIR, "04_footprint_condition.pdf"), p1, height = 12, width = 5)
p2 <- PlotFootprint(combined, features = fp_tfs, group.by = "celltype", label = FALSE)
ggsave(file.path(FIG_DIR, "04_footprint_celltype.pdf"), p2, height = 12, width = 5)

saveRDS(combined, file.path(OUT_DIR, "rds", "04_chromvar.rds"))
writeLines(capture.output(sessionInfo()),
           file.path(OUT_DIR, "sessionInfo", "04_chromvar_footprint.txt"))
