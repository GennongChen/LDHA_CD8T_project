# 参考代码审查报告：CTSD_CD8T_project/CART_snRNAATACseq

审查对象：`scMultiome.R`（1313 行）、`GRN_vis.R`（174 行）。
总体评价：流程设计合理（Signac 标准路线：QC → RNA/ATAC 各自降维 → WNN → 注释 → chromVAR → DA peaks → GeneActivity → 差异与富集 → Pando GRN），可以直接作为 LDHA 项目的骨架；但存在 **2 个会直接出错的 bug**、若干可复现性缺陷和大量硬编码路径，移植时必须修正。

## 一、必须修复的 bug

### B1. 写错对象：motif 结果表实际写入的是 RNA 结果
```r
markers_motifs_ct <- presto:::wilcoxauc.Seurat(..., seurat_assay = 'chromvar')
write.csv(markers_rna_ct, ".../DEGs_celltype_nMotif.csv")   # ← 写的是 markers_rna_ct
```
`DEGs_celltype_nMotif.csv` 里存的是 RNA 的 DEG，不是 motif。下游 `topTFs()` 用的是内存里的正确对象所以图没错，但落盘的表是错的。**移植时改为 `write.csv(markers_motifs_ct, ...)`。**

### B2. 引用未定义变量 `up_links`
```r
peak_to_plot <- up_links$peak[1]
up_links$gene
```
`up_links` 从未定义（上文只有 `peak_links`），这两行执行会直接报错。删除或改为 `peak_links`。

## 二、可复现性缺陷（移植时必须补上）

1. **无 `set.seed()`**：PCA/UMAP/clustering 均有随机性，两次运行结果会有微小差异。所有脚本开头统一 `set.seed(SEED)`，SEED 放 config。
2. **无 `sessionInfo()`**：每个脚本结尾输出 `sessionInfo()` 到 `results/sessionInfo/`，版本不一致时第一时间可查。
3. **环境未锁定**：shebang 指向个人 conda env（`/home/chengennong/tools/mambaforge/envs/atac`），包版本完全不受控。用 `renv` 锁定。
4. **中间产物无版本标记**：`combined_seurat.rds` 覆盖式保存。建议产物带阶段编号（如 `02_combined_qc.rds`）。

## 三、脆弱写法（换数据集就会坏，移植时重构）

| 位置 | 问题 | 修正 |
|---|---|---|
| celltype 注释 | `factor(wsnn_res.0.1, levels=c(1,2,0,3,4), labels=c("Tnaive",...))` 硬编码簇号→细胞类型映射，簇数或顺序一变即错 | 先跑 marker 基因（CD8A/GZMB/MKI67/KLF2/LEF1/TUBB...），人工确认每个簇的身份后再在 config 里写映射表 |
| top5 marker | `mutate(group = factor(group, levels = c(0:9)))` 假设恰好 10 个簇 | 改为 `levels = as.character(sort(unique(cluster_id)))` |
| `si.hilight` | 含尾部空格（`"LDLRAD4 "`、`"AOAH "`），`match()` 静默返回 NA，`anno_mark` 会错位或报错 | `trimws()` 后再 match，并 `stopifnot(!any(is.na(...)))` 校验 |
| topTFs 列索引 | `ctmarkers_rna[, c(2, 11, 6, 7)]` 按列号取子集，presto 版本一变列序就变 | 按列名取：`select(group, feature, logFC, auc)` 等 |
| 注释 `top_n <- 60` | 注释写"每个 cluster 选 5 个"，实际取 60 | 无功能影响，改注释 |
| `p_prop1` | `scale_fill_manual()` 调用了两次，第二次覆盖第一次 | 删一个 |
| `JoinLayers(combined_seurat, assays="RNA", layers=..., new.layer="counts")` | Seurat 5 的 `JoinLayers` 没有 `new.layer` 参数，取决于版本可能报 unused argument | 直接 `combined_seurat <- JoinLayers(combined_seurat, assay="RNA")`，并断言 `length(Layers(combined_seurat[["RNA"]], "counts")) == 1` |
| ATAC 矩阵去 scaffold | `combined_seurat@assays$ATAC@counts <- ATAC_counts[keep_rows, ]` 直接改 slot | 用 `subset(combined_seurat, features = keep_rows)` 或在 CreateChromatinAssay 前过滤矩阵 |

## 四、方法学上需要讨论的点（不阻塞，但要有意识）

1. **每个条件 n=1（n116 vs n216 各一只/一个池）**：所有按 Condition 的差异检验（wilcoxauc、FindMarkers）把细胞当独立样本，实际只有 1 个生物学重复，p 值会严重偏乐观。若 LDHA 项目同样无重复，结论表述要限定为"探索性"；若有 ≥2 重复/供者，**强烈建议改 pseudobulk（按 sample 聚合 counts 后用 DESeq2/edgeR）做 Condition 层面的 DE**，单细胞级检验仅作辅助。
2. **QC 阈值硬编码**（nCount_ATAC 5e3–7e4、nCount_RNA 1e3–2.5e4、percent.mt<20）：这是核内 RNA（sn），20% 线粒体可接受，但上下限必须先画 violin 图按本数据重新定，不能照抄。
3. **未做批次/条件整合**：两个样本直接 merge 后 WNN。对纯 CD8 T 细胞、同平台数据可接受；若 LDHA 项目样本间批次明显，考虑 Harmony（RNA）或 Seurat integration 后再 WNN，但要记录这一选择。
4. **`FindMarkers(test.use="LR")` 的 DA peaks 用 `latent.vars='nCount_ATAC'` 校正深度**：正确，保留。
5. **Pando GRN**：`infer_grn()` 默认参数在 2 万+ 细胞上很慢且内存大；建议先 `sample` 下采样到 ~1 万细胞跑通，再决定是否全量。合并条件与分条件各建一次 GRN 的设计保留。
6. **`AddMotifs()` + `Footprint()`**：需要 BSgenome 全基因组索引，内存峰值高（>32 GB 常见），服务器上单独跑并监控内存。

## 五、GRN_vis.R 备注

- 逻辑可用；`set.seed(123)` 已固定布局，好。
- `slice_head(n = 100)` 注释写 top200，改注释。
- `rna_df` 被连续覆盖赋值三次（csv → txt），前两次是死代码，清理。
- `layout2` 依赖 `rna_motifs.txt` 的列名（`motif.logFC`），该文件由 scMultiome.R 生成，两个脚本间有隐式数据依赖——在 README/脚本头注释里写明运行顺序。
