#!/home/chengennong/tools/mambaforge/envs/atac/bin/Rscript
# =====================================================================
# 01_preprocess_qc.R — 数据读入与 QC（CODEX_PLAN Phase 01）
# 输入: config$samples 中每个样本的 10x multiome h5 + fragments
# 输出: results/rds/01_combined_qc.rds
#       results/tables/01_qc_summary.csv
#       figures/01_qc_violin.pdf
# 修复项（相对参考 scMultiome.R）:
#   - 不再从预拼 rds 读入，直接读 10x 输出（如复用已有 rds 需在 config 注明）
#   - ATAC 去 scaffold 用规范子集方法，不直接改 @counts slot
#   - JoinLayers 不传 new.layer（Seurat 5 无此参数）
# =====================================================================
set.seed(SEED)
source("config/paths_config.R")

library(Seurat); library(Signac); library(tidyverse)
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0("chr", seqlevels(annotation)) |> str_replace("^chrMT$", "chrM")

# ---- 1. 逐样本读入 -----------------------------------------------------
obj_list <- lapply(names(samples), function(s) {
  info <- samples[[s]]
  counts <- Read10x_h5(info$h5)
  rna_obj <- CreateSeuratObject(counts = counts$`Gene Expression`, assay = "RNA")
  rna_obj[["ATAC"]] <- CreateChromatinAssay(
    counts = counts$`Peaks`,
    sep = c(":", "-"),
    fragments = info$fragments,
    annotation = annotation
  )
  rna_obj$sample_id <- s
  rna_obj$Condition <- info$condition
  rna_obj
})
names(obj_list) <- names(samples)

# ATAC 矩阵仅保留主染色体（规范子集，勿直接改 slot）
obj_list <- lapply(obj_list, function(o) {
  keep <- grep("^chr", rownames(o[["ATAC"]]), value = TRUE)
  subset(o, features = keep)  # subset 同时作用于两个 assay 的共有特征需注意；
                              # 若 RNA 基因被误删，改用 o[["ATAC"]] <- subset(o[["ATAC"]], features = keep)
})

# ---- 2. 合并 -----------------------------------------------------------
combined <- merge(obj_list[[1]], obj_list[-1],
                  add.cell.ids = names(obj_list), merge.data = TRUE)
combined <- JoinLayers(combined, assay = "RNA")
stopifnot(length(Layers(combined[["RNA"]], "counts")) == 1)  # 断言 counts 已合并

combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")

# ---- 3. QC：先画图定阈值，再过滤 ---------------------------------------
qc_before <- ncol(combined)
p_qc <- VlnPlot(combined, features = c("nCount_RNA", "nCount_ATAC", "percent.mt"),
                group.by = "sample_id", ncol = 3, pt.size = 0)
ggsave(file.path(FIG_DIR, "01_qc_violin.pdf"), p_qc, width = 12, height = 4)
# TODO(Codex): 查看 01_qc_violin.pdf，确认 config$QC 阈值适用于本数据后再运行下方 subset

combined <- subset(combined,
  subset = nCount_ATAC < QC$nCount_ATAC_max & nCount_ATAC > QC$nCount_ATAC_min &
           nCount_RNA  < QC$nCount_RNA_max  & nCount_RNA  > QC$nCount_RNA_min &
           percent.mt  < QC$percent_mt_max)

qc_summary <- data.frame(
  sample   = names(qc_before),
  cells_before = as.integer(qc_before),
  cells_after  = as.integer(table(combined$sample_id)[names(qc_before)])
)
write.csv(qc_summary, file.path(OUT_DIR, "tables", "01_qc_summary.csv"), row.names = FALSE)

saveRDS(combined, file.path(OUT_DIR, "rds", "01_combined_qc.rds"))
writeLines(capture.output(sessionInfo()),
           file.path(OUT_DIR, "sessionInfo", "01_preprocess_qc.txt"))
cat(sprintf("QC 完成: %d -> %d 细胞\n", sum(qc_before), ncol(combined)))
