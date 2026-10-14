#!/home/chengennong/tools/mambaforge/envs/atac/bin/Rscript
# Visualize a focused, exploratory Pando network after phases 06 and 07.
source("config/paths_config.R")
set.seed(SEED)
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidygraph)
  library(ggraph)
})
grn <- readRDS(file.path(OUT_DIR, "big", "07_Pando_combined_all.rds"))
de <- read.csv(file.path(OUT_DIR, "tables", "06_DE_condition.csv"))
names(grn) <- tolower(names(grn))
tf_col <- intersect(c("tf", "regulator", "motif"), names(grn))[1]
target_col <- intersect(c("target", "gene"), names(grn))[1]
effect_col <- intersect(c("estimate", "coefficient", "beta"), names(grn))[1]
p_col <- intersect(c("padj", "p_val_adj", "adjusted_p_value", "pvalue", "p_val"), names(grn))[1]
if (anyNA(c(tf_col, target_col, effect_col, p_col))) {
  stop("Unrecognized Pando coefficient columns: ", paste(names(grn), collapse = ", "))
}
edges <- data.frame(tf = as.character(grn[[tf_col]]),
                    target = as.character(grn[[target_col]]),
                    estimate = as.numeric(grn[[effect_col]]),
                    p = as.numeric(grn[[p_col]]))
if (p_col %in% c("pvalue", "p_val")) edges$p <- p.adjust(edges$p, "BH")
focus <- c("HIF1A", "NFKB1", "REL", "JUN", "FOS", "TBX21",
           "KLF2", "LEF1", "TCF7", "TOX2", "BATF", "LDHA")
edges <- edges |>
  dplyr::filter(tf %in% focus, is.finite(p), is.finite(estimate), p < 0.05) |>
  dplyr::arrange(p) |>
  dplyr::slice_head(n = 100)
write.csv(edges, file.path(OUT_DIR, "tables", "08_grn_focused_edges.csv"),
          row.names = FALSE)
pdf_path <- file.path(FIG_DIR, "08_grn_network.pdf")
if (nrow(edges) == 0) {
  pdf(pdf_path, width = 7, height = 5)
  plot.new()
  text(0.5, 0.5, "No focused Pando edges pass the exploratory threshold")
  dev.off()
} else {
  plot_edges <- edges |>
    dplyr::group_by(tf, target) |>
    dplyr::slice_min(p, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
  nodes <- data.frame(name = unique(c(plot_edges$tf, plot_edges$target)))
  nodes$type <- ifelse(nodes$name %in% plot_edges$tf, "TF", "target")
  nodes <- dplyr::left_join(nodes,
                            dplyr::distinct(de[, c("gene", "avg_log2FC")]),
                            by = c("name" = "gene"))
  graph_edges <- data.frame(from = match(plot_edges$tf, nodes$name),
                            to = match(plot_edges$target, nodes$name),
                            strength = abs(plot_edges$estimate),
                            effect_sign = ifelse(plot_edges$estimate > 0, "positive", "negative"))
  stopifnot(!anyNA(graph_edges$from), !anyNA(graph_edges$to))
  graph <- tidygraph::tbl_graph(nodes = nodes, edges = graph_edges, directed = TRUE)
  p <- ggraph(graph, layout = "fr") +
    geom_edge_link(aes(edge_colour = effect_sign, edge_width = strength),
                   arrow = grid::arrow(type = "closed", length = grid::unit(2, "mm")),
                   alpha = 0.7, show.legend = TRUE) +
    geom_node_point(aes(fill = type), shape = 21, size = 4,
                    color = "black") +
    geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    scale_edge_colour_manual(values = c(positive = "#984EA3",
                                        negative = "#4DAF4A")) +
    scale_fill_manual(values = c(TF = "#C9C0D3", target = "#F1F1F1")) +
    theme_void() +
    labs(title = "Exploratory CD8 CAR-T GRN; one sample per condition")
  ggsave(pdf_path, p, width = 11, height = 8)
}
figs <- list.files(FIG_DIR, pattern = "\\.pdf$", full.names = FALSE)
write.csv(data.frame(figure = figs),
          file.path(OUT_DIR, "tables", "figure_index.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()),
           file.path(OUT_DIR, "sessionInfo", "08_grn_vis.txt"))
