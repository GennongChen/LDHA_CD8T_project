# Configuration for the 26081226 CD8 CAR-T RNA/ATAC experiment.
# Run scripts from the repository root.
PROJECT_DIR <- normalizePath(getwd(), mustWork = TRUE)
DATA_DIR <- normalizePath(file.path(PROJECT_DIR, "..", "data"), mustWork = TRUE)
OUT_DIR <- file.path(PROJECT_DIR, "results")
FIG_DIR <- file.path(PROJECT_DIR, "figures")
for (p in c(file.path(OUT_DIR, "rds"), file.path(OUT_DIR, "tables"),
            file.path(OUT_DIR, "sessionInfo"), file.path(OUT_DIR, "big"), FIG_DIR)) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

samples <- list(
  CTRL = list(
    rna = file.path(DATA_DIR, "matrix/26081226_ctrl_GE/filtered_feature_bc_matrix"),
    atac = file.path(DATA_DIR, "matrix/26081226_ctrl_arc/filtered_peaks_bc_matrix"),
    fragments = file.path(DATA_DIR, "fragments/26081226_ctrl_arc/ctrl_A_fragments.tsv.gz"),
    atac_metrics = file.path(DATA_DIR, "bed/26081226_ctrl_arc/per_barcode_metrics.csv"),
    condition = "CTRL"
  ),
  OE = list(
    rna = file.path(DATA_DIR, "matrix/26081226_oe_GE/filtered_feature_bc_matrix"),
    atac = file.path(DATA_DIR, "matrix/26081226_oe_arc/filtered_peaks_bc_matrix"),
    fragments = file.path(DATA_DIR, "fragments/26081226_oe_arc/oe_A_fragments.tsv.gz"),
    atac_metrics = file.path(DATA_DIR, "bed/26081226_oe_arc/per_barcode_metrics.csv"),
    condition = "OE"
  )
)
COND_CTRL <- "CTRL"
COND_OE <- "OE"
GENOME_BUILD <- "hg38"
SEED <- 123L
RES_WNN <- 0.1
PCA_DIMS <- 1:30
LSI_DIMS <- 2:30
# Exploratory CD8 states from phase-03 RNA markers; suffix preserves WNN cluster.
celltype_map <- c(
  "0" = "IL7R_like_0",
  "1" = "Cytotoxic_1",
  "2" = "IL7R_like_2",
  "3" = "Cytotoxic_NR4A2_3",
  "4" = "Activated_IEG_4",
  "5" = "Cytotoxic_NR4A2_5",
  "6" = "Cycling_6",
  "7" = "Cytotoxic_7"
)
# Empirical cutoffs from 01_qc_metrics_before.csv: both samples have
# median mitochondrial fractions above 20%; 30% retains ~65-69% of cells.
QC <- list(nCount_ATAC_min = 500, nCount_ATAC_max = 50000,
           nCount_RNA_min = 1000, nCount_RNA_max = 25000,
           percent_mt_max = 30)
DA_TEST <- "LR"
DA_LATENT <- "nCount_ATAC"
LOGFC_CUT <- 0.25
PADJ_CUT <- 0.05
MODULES_FILE <- file.path(PROJECT_DIR, "config", "gene_modules.txt")
