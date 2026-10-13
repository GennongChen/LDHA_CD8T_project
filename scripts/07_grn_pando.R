#!/home/chengennong/tools/mambaforge/envs/atac/bin/Rscript
# =====================================================================
# 07_grn_pando.R — Pando GRN 推断（Phase 07，内存峰值高）
# 输入: results/rds/06_scored.rds
# 输出: results/tables/07_Pando_combined.csv, 07_Pando_<cond>.csv
#       figures/07_pando_combined.pdf, 07_pando_<cond>.pdf
# 策略: 先下采样 ~1 万细胞跑通，确认内存可行后再决定全量（config 中 N_PANDO 控制）
# =====================================================================
set.seed(SEED)
source("config/paths_config.R")
library(Seurat); library(Signac); library(Pando)
library(BSgenome.Hsapiens.UCSC.hg38)

combined <- readRDS(file.path(OUT_DIR, "rds", "06_scored.rds"))
N_PANDO <- 10000  # TODO(Codex): 跑通后按服务器内存决定是否改为 ncol(combined)

run_pando <- function(obj, label) {
  if (ncol(obj) > N_PANDO) {
    obj <- obj[, sample(colnames(obj), N_PANDO)]
    cat(sprintf("[%s] 下采样至 %d 细胞\n", label, ncol(obj)))
  }
  muo <- initiate_grn(obj, peak_assay = "ATAC", rna_assay = "RNA")
  muo <- find_motifs(muo, pfm = getMatrixSet(JASPAR2020, opts = list(species = 9606, all_versions = FALSE)),
                     genome = BSgenome.Hsapiens.UCSC.hg38)
  muo <- infer_grn(muo)
  muo <- find_modules(muo)
  write.csv(muo@grn@networks$glm_network@modules@meta,
            file.path(OUT_DIR, "tables", sprintf("07_Pando_%s.csv", label)), row.names = FALSE)
  pdf(file.path(FIG_DIR, sprintf("07_pando_%s.pdf", label)), width = 8, height = 6)
  print(plot_network_graph(get_network_graph(muo), layout = "fr"))
  dev.off()
  muo
}

# 合并条件 GRN + 分条件 GRN（与参考设计一致）
muo_combined <- run_pando(combined, "combined")
muo_ctrl <- run_pando(subset(combined, Condition == COND_CTRL), COND_CTRL)
muo_ko   <- run_pando(subset(combined, Condition == COND_KO),   COND_KO)

# TODO(Codex): 核心 TF 子网络图（参考 get_tf_network + plot_tf_network，
# TF 列表建议 HIF1A/NFKB1/REL/JUN/FOS + KLF2/LEF1，与 Phase 04 footprint 一致）
saveRDS(list(combined = muo_combined, ctrl = muo_ctrl, ko = muo_ko),
        file.path(OUT_DIR, "rds", "07_pando_objects.rds"))
writeLines(capture.output(sessionInfo()),
           file.path(OUT_DIR, "sessionInfo", "07_grn_pando.txt"))
