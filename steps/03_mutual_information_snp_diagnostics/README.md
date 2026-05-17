# Step 03 — Mutual-information SNP diagnostics

Adds an information-theoretic SNP ranking/diagnostic layer to the Step 02 panel resources.

This step is supportive and does not replace the primary Step 02 diagnostic SNP ranking.

## Inputs

```text
results/01_qc/metadata_clean.tsv
results/02_snp_panels/
```

## Main outputs

MI-ranking files written into the Step 02/03 panel resources under:

```text
results/02_snp_panels/
```

## Run

```bash
Rscript steps/03_mutual_information_snp_diagnostics/run.R
```
