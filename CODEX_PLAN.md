# CODEX 执行计划：LDHA_CD8T_project scMultiome 分析

> 执行者：Codex（服务器端）。审查者：Biomni（通过 GitHub）。
> 参考流程：`CTSD_CD8T_project/CART_snRNAATACseq/scMultiome.R`（审查结论见 `docs/code_review_reference.md`，下称"审查报告"，其中 B1/B2 等 bug 移植时必须修复）。

## 全局规则（每个阶段都适用）

1. **分支**：所有工作在 `codex/main` 分支进行，每阶段 ≥1 次 commit，阶段结束必须 push。commit message 格式：`[phase-XX] 一句话说明`。
2. **配置唯一来源**：所有路径、阈值、样本信息只写在 `config/paths_config.R`，脚本一律 `source("config/paths_config.R")`。禁止在脚本里出现绝对路径。
3. **可复现**：每个脚本开头 `set.seed(SEED)`；结尾 `writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo", "<script>.sessionInfo.txt"))`。
4. **产物落盘**：中间 Seurat 对象存 `results/rds/`（带阶段编号），表格存 `results/tables/`，图存 `figures/`（pdf）。大文件（rds、fastq、fragments）加入 `.gitignore`，不入库。
5. **检查点**：每个阶段末尾的 "✅ 检查点" 完成后才允许进入下一阶段；push 后在 PR/commit 里 @ Biomni 审查。

---

## Phase 00 — 环境与骨架（预计 0.5 天）

- [ ] 在服务器 clone 仓库。**运行环境使用现有 mamba 环境 `/home/chengennong/tools/mambaforge/envs/atac`**：所有脚本 shebang 写 `#!/home/chengennong/tools/mambaforge/envs/atac/bin/Rscript`，用 `.../envs/atac/bin/Rscript scripts/XX.R` 执行。
- [ ] 在该环境内安装 renv 后 `renv::init()`，把该环境当前的包版本锁定进 `renv.lock`（`renv::snapshot()` 只锁项目实际用到的包）。
- [ ] 安装并锁定包版本：`Seurat`(≥5), `Signac`, `EnsDb.Hsapiens.v86`, `BSgenome.Hsapiens.UCSC.hg38`, `tidyverse`, `presto`, `chromVAR`, `JASPAR2020`, `TFBSTools`, `motifmatchr`, `ggrepel`, `patchwork`, `pheatmap`, `ComplexHeatmap`, `circlize`, `clusterProfiler`, `org.Hs.eg.db`, `ggsignif`, `Pando`, `tidygraph`, `ggraph`。
- [ ] `renv::snapshot()` 后提交 `renv.lock`。
- [ ] 写 `scripts/00_environment.R`：加载全部包、打印版本表、检查 BSgenome 与 EnsDb 可用，输出到 `results/sessionInfo/00_environment.txt`。
- [ ] 按 LDHA 数据实际路径填写 `config/paths_config.R`（样本 h5、fragments、条件标签、QC 阈值占位）。

**✅ 检查点 0**：`renv.lock` + `00_environment.R` 运行日志已 push；Biomni 确认版本组合（尤其 Seurat 5.x 与 Signac 兼容性）。

## Phase 01 — 数据读入与 QC（1 天）

- [ ] `scripts/01_preprocess_qc.R`：从 10x multiome 输出直接读入（`Read10x_h5` 取 RNA + ATAC 两个 assay；不要像参考代码那样从预先拼好的 rds 读，除非服务器上已有等价对象——若复用已有 rds，在 config 里注明其来源）。
- [ ] 每个样本：建 Seurat 对象 → `CreateChromatinAssay`（挂 fragments 与 annotation）→ 加 `Condition` 标签。
- [ ] ATAC 矩阵仅保留 `^chr` 主染色体（用规范子集方法，见审查报告第三节）。
- [ ] `merge()` + `JoinLayers`（修复审查报告中的 `new.layer` 问题）。
- [ ] QC：先画 `VlnPlot`（nCount_RNA / nCount_ATAC / percent.mt，按样本分面）到 `figures/qc_violin.pdf`，**据此确定阈值**写入 config，再 `subset()`。参考起点：nCount_ATAC 5e3–1e5、nCount_RNA 1e3–2.5e4、percent.mt<20（snRNA 可放宽）。
- [ ] 输出 `results/rds/01_combined_qc.rds`、QC 过滤前后细胞数表 `results/tables/01_qc_summary.csv`。

**✅ 检查点 1**：QC 图 + 细胞数表 push；Biomni 确认阈值合理性（双联体比例、ATAC 深度分布）。

## Phase 02 — 降维、聚类与 WNN（0.5 天）

- [ ] `scripts/02_rna_atac_wnn.R`：RNA（NormalizeData→HVG→ScaleData→PCA 1:50）与 ATAC（TFIDF→FindTopFeatures(q0)→SVD，UMAP 用 2:50，跳过与深度相关的第一维）。
- [ ] `CellCycleScoring`（s.genes/g2m.genes），记录各期比例，若条件间细胞周期组成差异大，在报告中注明（LDHA 与乳酸/增殖相关，这本身可能是生物学信号，不要盲目回归掉）。
- [ ] `FindMultiModalNeighbors`（pca 1:50 + lsi 2:50）→ `RunUMAP`（wnn.umap）→ `FindClusters`（wsnn, algorithm=3）。resolution 从 0.1 起扫 0.1/0.15/0.2，各出一张 UMAP 供选择。
- [ ] 输出 `02_wnn.rds`、`figures/02_umap_*.pdf`。

**✅ 检查点 2**：UMAP 与簇数 push，Biomni 确认无明显的按样本分裂（批次主导）再继续。

## Phase 03 — 注释与 marker（0.5–1 天）

- [ ] `scripts/03_annotation_markers.R`：presto `wilcoxauc` 按 WNN 簇做 RNA marker；dotplot 展示 CD8 T 状态基因（KLF2/LEF1/TCF7=naive/stem, GZMB/PRF1/NKG7=effector, MKI67/TOP2A=cycling, CXCL13/TOX2=exhausted/Tex, NR4A1/2, IFNG, TBX21, HIF1A 及糖酵解相关 HK2/SLC16A3 等）。
- [ ] **簇→细胞类型映射写在 config**（`celltype_map = c("0"="Teff", ...)`），由人工看 marker 后确定；禁止照抄参考代码的 levels=c(1,2,0,3,4)。
- [ ] 比例图（stacked bar + grouped bar）、FeaturePlot（KLF2/RGS1/MKI67/JUN/NFKB1 + LDHA 本身）。
- [ ] 输出 `03_annotated.rds`、`results/tables/03_markers_clusters.csv`。

**✅ 检查点 3**：注释映射表 push，Biomni 复核每个簇的身份证据。

## Phase 04 — chromVAR 与 footprint（1 天）

- [ ] `scripts/04_chromvar_footprint.R`：JASPAR2020（species=9606, all_versions=FALSE）→ `CreateMotifMatrix`（基于 `granges(combined_seurat)` 当前 peaks）→ `RunChromVAR`。
- [ ] 按 celltype 和 Condition 分别做 motif 的 wilcoxauc（**修复 B1：写盘对象用 motif 结果**）。
- [ ] `topTFs()` 重写为按列名取子集（审查报告第三节）。
- [ ] `AddMotifs` + `Footprint`（TF 列表按 LDHA 项目调整，建议含 HIF1A/ARNTL 糖酵解轴 + NFKB/AP-1/KLF2 轴）；`PlotFootprint` 按 Condition 与 celltype 各出一张。
- [ ] 输出 `04_chromvar.rds`、`results/tables/04_motifs_*.csv`、`figures/04_footprint_*.pdf`。

**✅ 检查点 4**：差异 motif 表 push；Biomni 检查 HIF1A/糖酵解相关 motif 是否在 LDHA 扰动侧富集（这是本项目的核心预期之一）。

## Phase 05 — DA peaks 与 GeneActivity（1 天）

- [ ] `scripts/05_da_peaks_gene_activity.R`：`FindMarkers(test.use="LR", latent.vars="nCount_ATAC")` 按 Condition；删除参考代码中未定义的 `up_links`（B2）。
- [ ] peak-gene link：`Links()` 关联 DA peaks 与基因，输出 `results/tables/05_DE_peak_gene.csv`。
- [ ] `GeneActivity()`（extend.upstream=2000, downstream=0 起步，参考代码 3000/100 可在 config 里保留为参数）→ NormalizeData/ScaleData → `geneact` assay。
- [ ] 三张热图（geneact / GEX / motif，ComplexHeatmap，anno_mark 前先 `trimws` 并校验 match 成功）。
- [ ] 输出 `05_geneact.rds` 与对应表格。

**✅ 检查点 5**：DA peaks 表 + 热图 push。

## Phase 06 — 条件间 DE 与富集（0.5–1 天）

- [ ] `scripts/06_de_enrichment.R`：按 Condition 的 RNA DE（presto + Seurat FindMarkers 双口径），火山图。
- [ ] **若 LDHA 项目有 ≥2 生物学重复**：增加 pseudobulk 分析（按 sample×condition 聚合 counts → DESeq2），作为主结论；单细胞级检验放补充。n=1 则在脚本头注释明确"探索性，无生物学重复"。
- [ ] GO/KEGG（clusterProfiler，SYMBOL→ENTREZID 用 `bitr`），up/down 分侧。
- [ ] AddModuleScore：模块基因集放 `config/gene_modules.txt`（参考代码的 modules.txt 结构：每列一个通路），小提琴 + ggsignif；LDHA 项目建议加入糖酵解、缺氧反应、T 细胞耗竭、记忆/效应模块。
- [ ] 输出 `results/tables/06_DE_*.csv`、`figures/06_*.pdf`。

**✅ 检查点 6**：DE 表 + 富集图 push；Biomni 解读 LDHA 扰动的转录表型。

## Phase 07 — Pando GRN（1–2 天，内存峰值高）

- [ ] `scripts/07_grn_pando.R`：先下采样 ~1 万细胞跑通全流程（`set.seed` 后 `sample(colnames(...))`），确认内存可行后再决定全量。
- [ ] 合并条件 + 分条件（LDHA / 对照）各建一次 GRN；`find_motifs` → `infer_grn` → `find_modules` → 网络图。
- [ ] 输出 `results/tables/07_Pando*.csv`、`figures/07_pando*.pdf`。

**✅ 检查点 7**：GRN 网络表 push；Biomni 检查 HIF1A/NFKB/AP-1 等核心 TF 的调控关系。

## Phase 08 — GRN 可视化（0.5 天）

- [ ] `scripts/08_grn_vis.R`：基于参考 `GRN_vis.R` 改写（清理死代码、修正注释），聚焦 LDHA/糖酵解轴 + 效应/耗竭 TF 的子网络。
- [ ] 最终把全部图表清单写入 `results/tables/figure_index.csv`。

**✅ 检查点 8（终审）**：Biomni 通读全部结果，出具结果解读报告。

---

## 风险与回退

| 风险 | 应对 |
|---|---|
| Pando/BSgenome 内存超限 | 下采样；Footprint 分 TF 批次跑 |
| Seurat 5 与 Signac 版本不兼容 | 以 Phase 00 的版本表为准，锁死在 renv.lock |
| 簇数与参考项目不同导致注释映射失效 | 映射表放 config，Phase 03 人工确认 |
| 无生物学重复 | 所有 Condition 级差异结论标注"探索性"，主推 motif/peak 一致性证据 |
