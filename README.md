# LDHA_CD8T_project

LDHA 过表达（OE）与对照（CTRL）CD8 CAR-T 细胞的单细胞 RNA + ATAC 分析。

本批次每组只有一个生物学样本。条件与样本完全重合；条件间差异检验及富集结果只用于探索，不能据此推断可重复的 LDHA 因果效应。

分析流程改编自 [CTSD_CD8T_project/CART_snRNAATACseq](https://github.com/GennongChen/CTSD_CD8T_project/tree/main/CART_snRNAATACseq)（Signac + Seurat WNN + chromVAR + Pando GRN）。

## 目录结构

```
LDHA_CD8T_project/
├── README.md
├── CODEX_PLAN.md              # Codex 执行计划（分阶段，含检查点）
├── docs/
│   ├── code_review_reference.md   # 对 CTSD 参考代码的审查结论
│   ├── annotation_notes.md        # CD8 状态注释依据
│   └── results_report.md          # 本批次结果与限制
├── config/
│   ├── paths_config.R         # 路径、样本与分析参数
│   └── gene_modules.txt       # RNA 模块基因集
├── scripts/                   # 按阶段编号的分析脚本
│   ├── 00_environment.R
│   ├── 01_preprocess_qc.R
│   ├── 01b_atac_metrics.R
│   ├── 02_rna_atac_wnn.R
│   ├── 03_annotation_markers.R
│   ├── 04_chromvar_footprint.R
│   ├── 05_da_peaks_gene_activity.R
│   ├── 06_de_enrichment.R
│   ├── 07_grn_pando.R
│   ├── 08_grn_vis.R
│   ├── 09_cross_modal_volcano.R
│   └── 10_condition_coverage.R
├── results/                   # 表格及本地 RDS 缓存；大 RDS 不入库
└── figures/                   # PDF 图表
```

## 工作分工与分支

- **Codex**：负责代码审核、执行、结果核查及 Git pull/push；已在 `codex/main` 分阶段提交。
- **Biomni**：仓库原计划中的独立复核者；本次尚未执行该复核。
- 通过 PR 将 `codex/main` 合入 `main`；`renv.lock` 记录本次使用的 R 包版本。

## 快速开始（Codex）

```bash
# 在项目根目录运行；输入为 ../data 中独立的 RNA/ATAC 矩阵与 fragments。
# 数据路径与参数见 config/paths_config.R。
/home/chengennong/tools/mambaforge/envs/atac/bin/Rscript scripts/00_environment.R
/home/chengennong/tools/mambaforge/envs/atac/bin/Rscript scripts/01_preprocess_qc.R
/home/chengennong/tools/mambaforge/envs/atac/bin/Rscript scripts/01b_atac_metrics.R
/home/chengennong/tools/mambaforge/envs/atac/bin/Rscript scripts/02_rna_atac_wnn.R
/home/chengennong/tools/mambaforge/envs/atac/bin/Rscript scripts/03_annotation_markers.R
/home/chengennong/tools/mambaforge/envs/atac/bin/Rscript scripts/04_chromvar_footprint.R
/home/chengennong/tools/mambaforge/envs/atac/bin/Rscript scripts/05_da_peaks_gene_activity.R
/home/chengennong/tools/mambaforge/envs/atac/bin/Rscript scripts/06_de_enrichment.R
/home/chengennong/tools/mambaforge/envs/atac/bin/Rscript scripts/07_grn_pando.R
/home/chengennong/tools/mambaforge/envs/atac/bin/Rscript scripts/08_grn_vis.R
/home/chengennong/tools/mambaforge/envs/atac/bin/Rscript scripts/09_cross_modal_volcano.R
/home/chengennong/tools/mambaforge/envs/atac/bin/Rscript scripts/10_condition_coverage.R
```

截至 2026-10-14，Codex 已完成 00–10 阶段。第 09 阶段绘制 RNA、gene activity、motif 火山图，以及 TF RNA 与 motif 变化的散点图；motif 的横轴为 chromVAR 偏离分数均值差，并非 log2 fold change。第 10 阶段绘制 marker 和候选 TF 靶基因的 ATAC coverage 与 RNA 表达轨道。结果摘要见 [docs/results_report.md](docs/results_report.md)，全部图表见 [figure_index.csv](results/tables/figure_index.csv)。中间 RDS 只保存在分析服务器；若输入数据或参数变化，需清除相应阶段的缓存 RDS 后重跑。独立 Biomni 复核与合并到 `main` 尚未进行。
