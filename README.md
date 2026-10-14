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
│   └── code_review_reference.md   # 对 CTSD 参考代码的审查结论
├── config/
│   └── paths_config.R         # 唯一需要按服务器实际情况修改的配置文件
├── scripts/                   # 按阶段编号的分析脚本（Codex 填充）
│   ├── 00_environment.R
│   ├── 01_preprocess_qc.R
│   ├── 02_rna_atac_wnn.R
│   ├── 03_annotation_markers.R
│   ├── 04_chromvar_footprint.R
│   ├── 05_da_peaks_gene_activity.R
│   ├── 06_de_enrichment.R
│   ├── 07_grn_pando.R
│   └── 08_grn_vis.R
├── results/                   # 运行产物（rds、表格），不入库大文件
└── figures/                   # 图表输出（pdf/svg）
```

## 双 agent 协作约定

- **Codex（服务器）**：负责装环境、跑数据、按 `CODEX_PLAN.md` 分阶段实现脚本。
- **Biomni（沙盒）**：负责代码审查、方法学把关、结果解读；通过 GitHub 读取已推送的代码与结果表。
- 每个 agent 只在自己的分支工作（`codex/...` / `biomni/...`），通过 PR 合入 `main`；小步提交、每阶段结束必须 push。
- 环境由 `renv` 锁定（`renv.lock` 入库），任何依赖变更走 PR，另一方 `renv::restore()` 同步。

## 快速开始（Codex）

```bash
# 1. 在项目根目录，使用 mamba 环境 atac
# 2. 输入为 ../data/matrix 下每组独立的 RNA/ATAC Matrix Market 目录，
#    以及 ../data/fragments 下对应的 fragments 文件；见 config/paths_config.R
# 3. 按顺序运行 scripts/00 -> 08，每阶段核查结果后继续
/home/chengennong/tools/mambaforge/envs/atac/bin/Rscript scripts/00_environment.R
/home/chengennong/tools/mambaforge/envs/atac/bin/Rscript scripts/01_preprocess_qc.R
```
