# Step 03 — Mutual-information SNP diagnostics

Adds an information-theoretic SNP ranking/diagnostic layer to the Step 02 panel resources.

This step is supportive and does not replace the primary Step 02 diagnostic SNP ranking used for hierarchical origin assignment.

## Inputs

```text
results/01_qc/metadata_clean.tsv
results/02_snp_panels/panels_index.tsv
results/02_snp_panels/panels/<panel_id>/snp_matrix.rds
```

## Main outputs

```text
results/03_mi_diagnostics/panels_mi_index.tsv
results/03_mi_diagnostics/panels_mi_index.rds
results/03_mi_diagnostics/step03_run_info.rds
```

For compatibility with downstream diagnostic-panel development, per-panel MI tables are also written next to the Step 02 panel resources:

```text
results/02_snp_panels/panels/<panel_id>/snp_mi.tsv
results/02_snp_panels/panels/<panel_id>/snp_mi.rds
```

## Run

```bash
Rscript steps/03_mutual_information_snp_diagnostics/run.R
```
