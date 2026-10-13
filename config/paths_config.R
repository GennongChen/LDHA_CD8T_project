# =====================================================================
# config/paths_config.R —— LDHA_CD8T_project 唯一需要按服务器修改的文件
# 所有脚本 source 本文件；禁止在其他脚本中出现绝对路径或硬编码阈值
# =====================================================================

# ---- 运行环境（mamba/conda）------------------------------------------
# 所有脚本用该环境下的 R 执行：
#   /home/chengennong/tools/mambaforge/envs/atac/bin/Rscript scripts/01_preprocess_qc.R
# shebang 统一写：#!/home/chengennong/tools/mambaforge/envs/atac/bin/Rscript
R_ENV_BIN <- "/home/chengennong/tools/mambaforge/envs/atac/bin"

PROJECT_DIR <- "/path/to/LDHA_CD8T_project"     # 服务器上仓库的绝对路径
DATA_DIR    <- "/path/to/raw_multiome"          # 10x multiome cellranger-arc 输出根目录

OUT_DIR <- file.path(PROJECT_DIR, "results")
FIG_DIR <- file.path(PROJECT_DIR, "figures")
dir.create(file.path(OUT_DIR, "rds"),         recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "tables"),      recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "sessionInfo"), recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- 样本与条件 -------------------------------------------------------
# condition: 对照组标签在前，LDHA 扰动在后（与差异检验 ident.1/ident.2 对应）
samples <- list(
  CTRL = list(
    h5        = file.path(DATA_DIR, "ctrl/filtered_feature_bc_matrix.h5"),
    fragments = file.path(DATA_DIR, "ctrl/alt_splice_fragments.tsv.gz"),  # 按实际文件名改
    condition = "sgNTC"
  ),
  LDHA = list(
    h5        = file.path(DATA_DIR, "ldha/filtered_feature_bc_matrix.h5"),
    fragments = file.path(DATA_DIR, "ldha/alt_splice_fragments.tsv.gz"),
    condition = "sgLDHA"
  )
)
COND_CTRL <- "sgNTC"    # ident.2
COND_KO   <- "sgLDHA"   # ident.1

# ---- 参考基因组 -------------------------------------------------------
GENOME_BUILD <- "hg38"
ENSDB_PKG    <- "EnsDb.Hsapiens.v86"
BSGENOME_PKG <- "BSgenome.Hsapiens.UCSC.hg38"

# ---- QC 阈值（Phase 01 画完 violin 后按本数据修改，勿照抄参考项目）----
QC <- list(
  nCount_ATAC_min = 5e3,
  nCount_ATAC_max = 1e5,
  nCount_RNA_min  = 1000,
  nCount_RNA_max  = 25000,
  percent_mt_max  = 20   # snRNA 核内数据可放宽
)

# ---- 聚类与注释 -------------------------------------------------------
SEED        <- 123
RES_WNN     <- 0.1        # FindClusters resolution（备选 0.1/0.15/0.2）
PCA_DIMS    <- 1:50
LSI_DIMS    <- 2:50       # 跳过与测序深度相关的第一维
# Phase 03 人工确认后填写：簇号 -> 细胞类型
celltype_map <- c(
  # "0" = "Tnaive",
  # "1" = "Teff",
  # "2" = "Tcycling"
)

# ---- 差异分析参数 -----------------------------------------------------
DA_TEST      <- "LR"
DA_LATENT    <- "nCount_ATAC"
LOGFC_CUT    <- 0.25
PADJ_CUT     <- 0.05

# ---- 模块打分基因集（每列一个通路，列内逗号/空白分隔基因）------------
MODULES_FILE <- file.path(PROJECT_DIR, "config", "gene_modules.txt")
