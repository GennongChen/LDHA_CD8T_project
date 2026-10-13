#!/home/chengennong/tools/mambaforge/envs/atac/bin/Rscript
# =====================================================================
# 02_rna_atac_wnn.R — RNA/ATAC 降维 + WNN 整合（CODEX_PLAN Phase 02）
# 输入: results/rds/01_combined_qc.rds
# 输出: results/rds/02_wnn.rds
#       figures/02_umap_res{0.1,0.15,0.2}.pdf
# 注意: CellCycleScoring 仅记录比例；LDHA 与增殖/乳酸相关，周期差异本身
#       可能是生物学信号，默认不做 vars.to.regress（TODO: 与审查者确认）
# =====================================================================
set.seed(SEED)
source("config/paths_config.R")
library(Seurat); library(Signac); library(tidyverse)

combined <- readRDS(file.path(OUT_DIR, "rds", "01_combined_qc.rds"))

# ---- RNA ----
DefaultAssay(combined) <- "RNA"
combined <- NormalizeData(combined, verbose = FALSE)
combined <- FindVariableFeatures(combined, verbose = FALSE)
combined <- ScaleData(combined, features = Features(combined), verbose = FALSE)
combined <- RunPCA(combined, verbose = FALSE)
combined <- CellCycleScoring(combined, s.features = cc.genes$s.genes,
                             g2m.features = cc.genes$g2m.genes)
print(table(combined$Condition, combined$Phase))  # 条件间周期组成，写入日志

# ---- ATAC（LSI 跳过与深度相关的第一维）----
DefaultAssay(combined) <- "ATAC"
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = "q0")
combined <- RunSVD(combined)

# ---- WNN ----
combined <- FindMultiModalNeighbors(
  combined, reduction.list = list("pca", "lsi"),
  dims.list = list(PCA_DIMS, LSI_DIMS), verbose = FALSE)
combined <- RunUMAP(combined, nn.name = "weighted.nn",
                    reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

# ---- 聚类：resolution 扫描 ----
for (r in c(0.1, 0.15, 0.2)) {
  combined <- FindClusters(combined, graph.name = "wsnn", algorithm = 3,
                           resolution = r, cluster.name = paste0("wsnn_res", r),
                           verbose = FALSE)
  p <- DimPlot(combined, reduction = "wnn.umap",
               group.by = paste0("wsnn_res", r), label = TRUE) + ggtitle(paste("res", r))
  ggsave(file.path(FIG_DIR, sprintf("02_umap_res%s.pdf", r)), p, width = 6, height = 5)
}
# TODO(Codex): 与审查者确认最终 resolution 后，在 config 中更新 RES_WNN

saveRDS(combined, file.path(OUT_DIR, "rds", "02_wnn.rds"))
writeLines(capture.output(sessionInfo()),
           file.path(OUT_DIR, "sessionInfo", "02_rna_atac_wnn.txt"))
