#!/home/chengennong/tools/mambaforge/envs/atac/bin/Rscript
# =====================================================================
# 00_environment.R — 环境检查（CODEX_PLAN Phase 00）
# 输入: 无
# 输出: results/sessionInfo/00_environment.txt（包版本表 + 依赖可用性）
# =====================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Signac); library(tidyverse)
  library(EnsDb.Hsapiens.v86); library(BSgenome.Hsapiens.UCSC.hg38)
  library(presto); library(chromVAR); library(JASPAR2020)
  library(TFBSTools); library(motifmatchr)
  library(ggrepel); library(patchwork); library(pheatmap)
  library(ComplexHeatmap); library(circlize)
  library(clusterProfiler); library(org.Hs.eg.db); library(ggsignif)
  library(Pando); library(tidygraph); library(ggraph)
})

here <- function(...) file.path(dirname(commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))] |> sub("^--file=", ""), ...))
# 若以 source() 方式运行，直接手动 setwd 到仓库根目录
if (!dir.exists("config")) setwd("/path/to/LDHA_CD8T_project")  # TODO(Codex): 改成实际路径

source("config/paths_config.R")

# 关键依赖可用性检查
stopifnot(
  "EnsDb 不可用"        = exists("EnsDb.Hsapiens.v86"),
  "BSgenome 不可用"     = exists("BSgenome.Hsapiens.UCSC.hg38"),
  "Seurat/Signac 版本"  = packageVersion("Seurat") >= "5.0"
)
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0("chr", seqlevels(annotation)) |> str_replace("^chrMT$", "chrM")
cat("annotation ranges:", length(annotation), "\n")

out <- file.path(OUT_DIR, "sessionInfo", "00_environment.txt")
writeLines(capture.output(sessionInfo()), out)
cat("环境检查通过，版本信息已写入:", out, "\n")
