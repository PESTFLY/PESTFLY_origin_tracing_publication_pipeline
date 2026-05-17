# Step 05 — Full-panel leave-one-out validation

Validates the full-panel Step 04 assignment logic using leave-one-out cross-validation on reference samples.

This step evaluates whether the hierarchical assignment logic correctly recovers known reference labels under the same multi-K stability logic used for query/intercept samples.

## Inputs

```text
results/01_qc/metadata_clean.tsv
results/02_snp_panels/
```

## Main outputs

```text
results/05_loo_validation/
```

## Run

```bash
Rscript steps/05_full_panel_loo_validation/run.R
```
