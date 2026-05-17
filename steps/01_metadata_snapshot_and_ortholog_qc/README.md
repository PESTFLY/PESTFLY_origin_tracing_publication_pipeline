# Step 01 — Metadata snapshot and ortholog QC

Cleans and harmonizes sample metadata, checks sample labels, summarizes ortholog/alignment availability, and creates the metadata snapshot used by downstream steps.

## Inputs

```text
data/000_input_data/metadata.xlsx
results/00_fasta/
data/000_input_data/phy/
```

## Main outputs

```text
results/01_qc/metadata_clean.tsv
results/01_qc/
```

## Run

```bash
Rscript steps/01_metadata_snapshot_and_ortholog_qc/run.R
```
