#!/home/chengennong/tools/mambaforge/envs/atac/bin/Rscript
# =====================================================================
# 06_de_enrichment.R — 条件间 DE + GO/KEGG + 模块打分（Phase 06）
# 输入: results/rds/05_geneact.rds, config/gene_modules.txt
# 输出: results/tables/06_DE_condition.csv, 06_GO_*.csv, 06_KEGG.csv
#       figures/06_volcano.pdf, 06_go_dot.pdf, 06_kegg_dot.pdf, 06_module_violin.pdf
# 统计学注意: 若每条件 n>=2 生物学重复，主结论用 pseudobulk+DESeq2（见 plan）；
#             n=1 时单细胞级检验仅作探索性证据，脚本头注明。
# =====================================================================
# 探索性分析声明：本数据每条件生物学重复数 = TODO(Codex: 填写)。
# 若 <2，所有 Condition 级差异结论仅具探索性。
set.seed(SEED)
source("config/paths_config.R")
library(Seurat); library(Signac); library(tidyverse)
library(clusterProfiler); library(org.Hs.eg.db); library(ggsignif)

combined <- readRDS(file.path(OUT_DIR, "rds", "05_geneact.rds"))

# ---- DE（presto + FindMarkers 双口径）----
de_presto <- presto:::wilcoxauc.Seurat(combined, group_by = "Condition",
                                       assay = "data", seurat_assay = "RNA")
DefaultAssay(combined) <- "RNA"
de_fm <- FindMarkers(combined, group.by = "Condition",
                     ident.1 = COND_KO, ident.2 = COND_CTRL,
                     test.use = "wilcox", min.pct = 0.05, logfc.threshold = 0.1)
de_fm$gene <- rownames(de_fm)
de_fm$regulate <- case_when(
  de_fm$avg_log2FC >=  LOGFC_CUT & de_fm$p_val_adj < PADJ_CUT ~ "Up",
  de_fm$avg_log2FC <= -LOGFC_CUT & de_fm$p_val_adj < PADJ_CUT ~ "Down",
  TRUE ~ "NS")
write.csv(de_fm, file.path(OUT_DIR, "tables", "06_DE_condition.csv"), row.names = FALSE)

# ---- 火山图 ----
p_vol <- ggplot(de_fm, aes(x = avg_log2FC, y = pmin(-log10(p_val_adj), 50))) +
  geom_point(aes(color = regulate), size = 1.5) +
  scale_color_manual(values = c(Up = "#BF8AC0", Down = "#B8DA68", NS = "grey")) +
  ggrepel::geom_text_repel(
    data = subset(de_fm, regulate != "NS" & abs(avg_log2FC) > 1 &
                    !str_detect(gene, "^MT-|^RPL|^RPS|\\.")),
    aes(label = gene), size = 3.5, max.overlaps = 10) +
  theme_bw() + labs(title = paste(COND_KO, "vs", COND_CTRL),
                    x = "log2FC", y = "-log10(adj.p)")
ggsave(file.path(FIG_DIR, "06_volcano.pdf"), p_vol, height = 4, width = 5)

# ---- GO / KEGG ----
up_genes   <- subset(de_fm, regulate == "Up")$gene
down_genes <- subset(de_fm, regulate == "Down")$gene
for (side in c("Up", "Down")) {
  genes <- get(side == "Up" ? "up_genes" : "down_genes")  # R<4.2 兼容: 用 if
  go <- enrichGO(gene = if (side == "Up") up_genes else down_genes,
                 OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "fdr")
  write.csv(go, file.path(OUT_DIR, "tables", sprintf("06_GO_%s.csv", side)), row.names = FALSE)
  ggsave(file.path(FIG_DIR, sprintf("06_go_%s_dot.pdf", side)),
         dotplot(go, showCategory = 9), height = 4, width = 5)
}
genes_id <- bitr(up_genes, "SYMBOL", "ENTREZID", org.Hs.eg.db)
kegg <- enrichKEGG(gene = genes_id$ENTREZID, organism = "hsa",
                   pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
write.csv(kegg, file.path(OUT_DIR, "tables", "06_KEGG.csv"), row.names = FALSE)
ggsave(file.path(FIG_DIR, "06_kegg_dot.pdf"), dotplot(kegg), height = 4, width = 6)

# ---- 模块打分（基因集: config/gene_modules.txt，每列一个通路）----
# TODO(Codex): 建议模块含 糖酵解/缺氧反应/T细胞耗竭/记忆/效应；
# 结构沿用参考代码 modules.txt（read.delim + 循环 AddModuleScore + ggsignif 小提琴）。
if (file.exists(MODULES_FILE)) {
  gene_sets <- read.delim(MODULES_FILE, check.names = FALSE)
  for (pw in colnames(gene_sets)) {
    genes <- na.omit(gene_sets[, pw]); genes <- unique(as.character(genes[genes != ""]))
    combined <- AddModuleScore(combined, features = list(genes), name = paste0(pw, "_score"))
  }
  # TODO(Codex): 补齐批量小提琴图（参考 scMultiome.R 760-800 行）
}

saveRDS(combined, file.path(OUT_DIR, "rds", "06_scored.rds"))
writeLines(capture.output(sessionInfo()),
           file.path(OUT_DIR, "sessionInfo", "06_de_enrichment.txt"))
