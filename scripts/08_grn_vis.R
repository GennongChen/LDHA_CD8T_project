#!/home/chengennong/tools/mambaforge/envs/atac/bin/Rscript
# =====================================================================
# 08_grn_vis.R — GRN 子网络可视化（Phase 08，基于参考 GRN_vis.R 改写）
# 输入: results/tables/07_Pando_combined.csv, 06_DE_condition.csv
# 输出: figures/08_grn_network.pdf
# 依赖: 必须在 06、07 之后运行（隐式数据依赖，见审查报告第五节）
# 修复项: 清理参考代码中 rna_df 的三次覆盖死代码；slice_head 注释修正
# =====================================================================
set.seed(SEED)
source("config/paths_config.R")
library(tidyverse); library(tidygraph); library(ggraph)

grn_df <- read.csv(file.path(OUT_DIR, "tables", "07_Pando_combined.csv"), row.names = 1)
de_df  <- read.csv(file.path(OUT_DIR, "tables", "06_DE_condition.csv"))

# TODO(Codex): 聚焦 TF 子网络——LDHA/糖酵解轴（HIF1A/ARNTL）+
# 效应/耗竭轴（TBX21/NFKB1/REL/JUN/FOS/TOX2/BATF）+ 记忆轴（KLF2/LEF1/TCF7）
gene_to_vis <- c("HIF1A", "NFKB1", "REL", "JUN", "FOS", "TBX21",
                 "KLF2", "LEF1", "TCF7", "TOX2", "BATF")
gene_to_vis <- intersect(gene_to_vis, unique(c(grn_df$tf, grn_df$target)))

grn_filtered <- grn_df |>
  filter(tf %in% gene_to_vis, padj < 0.05) |>
  arrange(padj) |>
  slice_head(n = 100) |>   # top 100 条边（按 padj）
  mutate(neg_log_padj = -log10(padj),
         regulation = ifelse(estimate > 0, "Activating", "Repressing"))

nodes <- tibble(name = unique(c(grn_filtered$tf, grn_filtered$target))) |>
  mutate(type = ifelse(name %in% grn_filtered$tf, "TF", "target"))
edges <- grn_filtered |> select(from = tf, to = target, neg_log_padj, regulation, estimate)
grn_graph <- tbl_graph(nodes, edges, directed = TRUE)

# 节点叠加条件间 DE 信息（TF 调控方向 vs 表达变化）
layout <- create_layout(grn_graph, layout = "stress") |>
  left_join(de_df |> select(gene, avg_log2FC) |> unique, by = c("name" = "gene"))

p <- ggraph(layout) +
  geom_edge_link(aes(alpha = neg_log_padj, color = regulation, linewidth = estimate)) +
  geom_node_point(aes(size = centrality_degree(mode = "in"), fill = type),
                  shape = 21, stroke = 0.3, color = "black") +
  geom_node_text(aes(label = ifelse(type == "TF", name, "")), size = 3, repel = TRUE) +
  scale_edge_color_manual(values = c(Activating = "#984EA3", Repressing = "#B8DA68")) +
  scale_edge_alpha_continuous(range = c(0.5, 1)) +
  scale_fill_manual(values = c(TF = "#984EA3", target = "#C9C0D3")) +
  theme_void() + labs(edge_color = "Regulation", fill = "Node type")
ggsave(file.path(FIG_DIR, "08_grn_network.pdf"), p, width = 6, height = 5, bg = "white")

# ---- 图表清单 ----
figs <- list.files(FIG_DIR, pattern = "\\.pdf$", full.names = FALSE)
write.csv(data.frame(figure = figs), file.path(OUT_DIR, "tables", "figure_index.csv"),
          row.names = FALSE)
writeLines(capture.output(sessionInfo()),
           file.path(OUT_DIR, "sessionInfo", "08_grn_vis.txt"))
