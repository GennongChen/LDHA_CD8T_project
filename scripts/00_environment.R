#!/usr/bin/env Rscript
source("config/paths_config.R")
set.seed(SEED)
required <- c("Seurat", "Signac", "Matrix", "EnsDb.Hsapiens.v86",
              "BSgenome.Hsapiens.UCSC.hg38", "presto", "chromVAR",
              "JASPAR2020", "TFBSTools", "motifmatchr", "ggrepel",
              "patchwork", "pheatmap", "ComplexHeatmap", "circlize",
              "clusterProfiler", "org.Hs.eg.db", "ggsignif", "Pando",
              "tidygraph", "ggraph")
versions <- vapply(required, function(p) {
  if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else NA_character_
}, character(1))
write.csv(data.frame(package = names(versions), version = unname(versions)),
          file.path(OUT_DIR, "tables", "00_package_versions.csv"), row.names = FALSE)
if (anyNA(versions)) stop("Missing packages: ", paste(names(versions)[is.na(versions)], collapse = ", "))
stopifnot(packageVersion("Seurat") >= "5.0")
stopifnot(all(vapply(samples, function(x) all(file.exists(unlist(x[c("rna", "atac", "fragments", "atac_metrics")]))), logical(1))))
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo", "00_environment.txt"))
message("Environment and input paths verified")
