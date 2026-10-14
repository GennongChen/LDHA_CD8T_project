# LDHA 过表达 CD8 CAR-T 单细胞 RNA + ATAC 分析报告

## 实验设计与解释范围

CTRL 和 LDHA 过表达（OE）各只有一个生物学样本，因此条件与样本完全重合。下述 RNA、ATAC、motif 和网络的组间差异均为探索性描述。以细胞为单位计算的校正后 p 值不能代表生物学重复间的可重复性，也不能证明 LDHA 的因果效应。

## 数据与质控

按 barcode 配对独立的 RNA 和 ATAC Matrix Market 矩阵，并在共同的 52,189 个 peak 上重新计数 fragments。CTRL 输入 4,887 个细胞，保留 3,097 个；OE 输入 5,054 个，保留 3,202 个。供应商给出的 FRiP 中位数分别为 0.822 和 0.860；ATAC fragments 中位数分别为 1,690 和 7,065.5。约四倍的深度差可能影响 ATAC 结构，即使差异可及性模型校正了单细胞 ATAC 计数也应谨慎解释。当前环境没有可用的双细胞检测包；质控设置了计数上限，但仍可能残留双细胞。

## CD8 状态与样本分离

RNA PCA、ATAC LSI 和 WNN 得到八个探索性的 CD8 CAR-T 状态，注释依据见 [annotation_notes.md](annotation_notes.md)。WNN 簇与样本高度相关：簇 0 的 OE 比例为 2.7%，簇 5 为 99.7%。因此细胞状态比例变化无法归因于 LDHA。G2/M 比例为 CTRL 332/3,097（10.7%）、OE 594/3,202（18.6%）。LDHA 的平均 log 标准化 RNA 表达为 CTRL 0.861、OE 0.954，检出比例为 64.8% 和 65.8%。仅凭这些 RNA 信号不能验证过表达效率。

## chromVAR 与 footprint

对 633 个 JASPAR 2020 motif 计算 chromVAR。OE 侧数值较高的 motif 包括 NRF1、CTCF 和 ELK3/ELK1/ELF2 等 ETS 家族；BATF、FOS 和 FOS::JUN 等 AP-1 motif 较低。HIF1A 和 NFKB1 在 presto 结果表中的 `logFC` 列分别为 0.250 和 0.069。这些数值受细胞状态组成和 ATAC 深度差异影响；motif 偏离分数不能直接证明 TF 蛋白活性或靶基因调控。六个 TF 的 footprint 已绘图并目视检查。

## 差异可及性与 peak–gene 关联

校正单细胞 ATAC 计数的探索性逻辑回归检验了 21,568 个 peak；细胞级校正后 p < 0.05 的有 380 个，其中 OE 较高 369 个、CTRL 较高 11 个。这些数目不是经生物学重复验证的差异 peak。对选定 RNA 可变基因计算的 peak–gene 关联产生 1,015 行记录，涉及 959 个受检 peak 和 384 个基因；其中 44 行对应达到上述细胞级阈值的 peak。GeneActivity、RNA 和 chromVAR 热图概括了八个状态。关联分数不能证明直接调控。

## RNA 差异、通路与模块

Seurat 单细胞检验在表达与 fold-change 过滤后评估 5,345 个基因。按校正后 p < 0.05 且 |log2FC| ≥ 0.25，OE 较高 219 个、CTRL 较高 721 个。LDHA 的 log2FC 为 0.220，低于该效应阈值。OE 较高基因的 GO 结果包含细胞杀伤和白细胞介导的细胞毒作用；KEGG 在线重试得到 OE 侧 41 条、CTRL 侧 12 条通路，OE 侧首位为抗原加工与呈递。OE 的细胞毒、增殖、糖酵解和缺氧模块评分中位数均略高于 CTRL，具体数值见 `results/tables/06_module_summary.csv`。这些都是本批次内部描述，不能分离条件与样本效应。

## Pando 调控网络

Pando 作为有界试运行，每套网络最多使用 2,000 个细胞，请求 511 个靶基因，其中 371 个进入模型拟合。合并、CTRL、OE 网络分别报告 80、46、88 个 TF 模块。合并网络的聚焦结果包含 33 行达到模型阈值的记录，涉及 FOS、JUN、NFKB1、REL、BATF、KLF2 和 TOX2；所选阈值下没有 HIF1A 边。网络图对每个 TF–靶基因对显示最显著的一条边，[完整聚焦边表](../results/tables/08_grn_focused_edges.csv)保留全部 33 行。条件网络间的差别可能由样本、下采样和模型选择造成，不能归因于 LDHA。

## RNA、gene activity 与 motif 火山图

第 09 阶段以 `min.pct = 0.05`、`logfc.threshold = 0` 重新计算完整的 RNA 和 gene activity 单细胞 Wilcoxon 差异表，分别覆盖 8,160 和 13,126 个基因。motif 图覆盖 633 个 JASPAR motif，使用 chromVAR deviation Z 分数的 OE−CTRL 均值差；它与原 presto 表中的 `logFC` 数值一致，但**不是 log2 fold change**，因为 deviation 分数可以为负。三图均以校正后 p < 0.05 和效应阈值着色（RNA 与 gene activity 的 |log2FC| ≥ 0.25；motif 的 |分数差| ≥ 0.25），每个方向标注校正后 p 最小的五个不同基因或 TF，并额外标注 RNA 和 gene activity 中的 LDHA。图上 −log10 校正后 p 值封顶于 50，避免数值下溢和极小 p 值挤占画面。

RNA 图中 OE 较高 219 个、CTRL 较高 721 个；gene activity 图中分别为 4,455 和 1,388 个；motif 图中分别为 148 和 32 个。LDHA RNA log2FC 为 0.220，gene activity log2FC 为 0.147，均低于着色效应阈值。gene activity 的大范围差异尤其可能受到两组 ATAC 测序深度差异影响。相关图见 [RNA 火山图](../figures/09_volcano_RNA.pdf)、[gene activity 火山图](../figures/09_volcano_geneact.pdf)、[motif 火山图](../figures/09_volcano_motif.pdf)；对应全量表和标注记录以 `09_` 为前缀保存在 `results/tables/`。

[RNA × motif 图](../figures/09_RNA_x_motif_logFC.pdf)比较 TF 基因的 RNA log2FC 与相同 TF motif 的 chromVAR 分数差。精确基因名匹配后保留 182 个 TF；一个 TF 对应多个 JASPAR motif 时，选校正后 p 最小者，再按绝对效应排序。两种细胞级校正后 p 均 < 0.05 的 TF 中，方向一致 12 个、相反 16 个；各组按两种绝对效应的乘积标注前六个。详细结果与所选 motif ID 见 [RNA × motif 表](../results/tables/09_RNA_x_motif.csv)。RNA 与 motif 的方向一致性只是相关描述，不能证明 TF 活性或 LDHA 调控作用。

## 复现与后续验证

分析代码位于 `scripts/`，参数和数据路径位于 `config/paths_config.R`，R 版本与包版本见 `renv.lock`。结果表和 PDF 图已提交到仓库；较大的 RDS 保存在分析服务器的 `results/rds/` 和 `results/big/`。若需验证 LDHA 效应，应补充独立生物学重复，并用每样本汇总的 pseudobulk 或相应样本级模型分析；还需用独立实验确认 LDHA 过表达效率。
