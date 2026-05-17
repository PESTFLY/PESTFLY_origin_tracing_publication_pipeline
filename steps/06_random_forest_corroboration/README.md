# Step 06 — Random Forest corroboration

Runs Random Forest classification as an independent corroborative layer for Step 04 origin assignments.

Random Forest results are not used to override Step 04. They are interpreted as support, weak support, inconclusive evidence, or caution depending on agreement and confidence.

## Inputs

```text
results/01_qc/metadata_clean.tsv
results/02_snp_panels/
results/04_origin_assignment/final_assignment.tsv
```

## Main outputs

```text
results/06_rf_corroboration/
```

## Run

```bash
Rscript steps/06_random_forest_corroboration/run.R
```
