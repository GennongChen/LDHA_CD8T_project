# Phase 03 exploratory CD8 CAR-T state labels

The WNN clusters are strongly associated with sample. CTRL and OE are each one
biological sample, so these labels describe expression patterns only. They do
not establish LDHA-driven state changes.

| WNN cluster | Label | Supporting RNA markers | Sample composition |
|---|---|---|---|
| 0 | IL7R_like_0 | IL7R | 2.7% OE |
| 1 | Cytotoxic_1 | PRF1, GZMB, NKG7, CD8A | 97.8% OE |
| 2 | IL7R_like_2 | IL7R | 17.9% OE |
| 3 | Cytotoxic_NR4A2_3 | NR4A2, PRF1, GZMB | 31.9% OE |
| 4 | Activated_IEG_4 | JUN, RGS1, NR4A1/2, IFNG | 72.5% OE |
| 5 | Cytotoxic_NR4A2_5 | NR4A2, PRF1, CD8A | 99.7% OE |
| 6 | Cycling_6 | MKI67, TOP2A | 60.7% OE |
| 7 | Cytotoxic_7 | GZMB, PRF1, NKG7 | 98.3% OE |

Cluster 6 has higher RNA depth; some cell-cycle marker enrichment may partly
reflect its higher UMI count. The complete marker scores and cluster QC are in
`results/tables/03_state_gene_markers.csv` and
`results/tables/03_cluster_qc.csv`.

LDHA RNA signal is only modestly higher in OE in the retained cells: mean
log-normalized expression 0.954 versus 0.861 in CTRL, with detection in
65.8% versus 64.8% of cells. This alone does not verify overexpression
efficiency; see `results/tables/03_LDHA_expression.csv`.
