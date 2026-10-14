#!/usr/bin/env Rscript
# Vendor FRiP and fragment summaries for the cells retained in phase 01.
source("config/paths_config.R")
set.seed(SEED)
suppressPackageStartupMessages(library(ggplot2))
obj <- readRDS(file.path(OUT_DIR, "rds", "01_combined_qc.rds"))
metrics <- do.call(rbind, lapply(names(samples), function(s) {
  x <- read.csv(samples[[s]]$atac_metrics)
  selected <- sub(paste0("^", s, "_"), "", colnames(obj)[obj$sample_id == s])
  idx <- match(selected, x$barcode)
  stopifnot(!anyNA(idx))
  x <- x[idx, c("barcode", "fragments_num", "fraction_frag_overlap_peak", "is_cell")]
  x$sample <- s
  x
}))
stopifnot(all(metrics$is_cell == 1))
write.csv(metrics, file.path(OUT_DIR, "tables", "01b_atac_metrics.csv"), row.names = FALSE)
summary <- do.call(rbind, lapply(split(metrics, metrics$sample), function(x) {
  data.frame(sample = x$sample[1], cells = nrow(x),
             median_fragments = median(x$fragments_num),
             median_FRiP = median(x$fraction_frag_overlap_peak),
             fraction_FRiP_below_0.1 = mean(x$fraction_frag_overlap_peak < 0.1))
}))
write.csv(summary, file.path(OUT_DIR, "tables", "01b_atac_qc_summary.csv"), row.names = FALSE)
p <- ggplot(metrics, aes(sample, fraction_frag_overlap_peak, fill = sample)) +
  geom_violin(trim = FALSE) + geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_hline(yintercept = 0.1, linetype = 2, color = "red") +
  labs(x = NULL, y = "Fraction of fragments in original peaks") + theme_classic()
ggsave(file.path(FIG_DIR, "01b_atac_FRiP.pdf"), p, width = 5, height = 4)
writeLines(capture.output(sessionInfo()),
           file.path(OUT_DIR, "sessionInfo", "01b_atac_metrics.txt"))
