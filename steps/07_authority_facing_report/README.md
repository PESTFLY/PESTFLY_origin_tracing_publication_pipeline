# Step 07 — Authority-facing report

Integrates Step 04 hierarchical assignments, Step 05 validation summaries, and Step 06 Random Forest corroboration into an authority-facing report.

This is the endpoint of the public manuscript workflow.

## Inputs

```text
results/04_origin_assignment/final_assignment.tsv
results/05_loo_validation/
results/06_rf_corroboration/
```

## Main outputs

```text
results/07_authority_report/FINAL_origin_tracing_authority_simplified.tsv
results/07_authority_report/FINAL_origin_tracing_authority_extended.tsv
results/07_authority_report/FINAL_origin_tracing_authority_REPORT.xlsx
```

## Run

```bash
Rscript steps/07_authority_facing_report/run.R
```
