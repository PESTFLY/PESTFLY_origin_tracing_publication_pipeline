# Results directory

Each step writes to a matching numbered result folder:

```text
results/00_fasta/             Step 00 output
results/01_qc/                Step 01 output
results/02_snp_panels/        Step 02 output
results/03_mi_diagnostics/    Step 03 output
results/04_origin_assignment/ Step 04 output
results/05_loo_validation/    Step 05 output
results/06_rf_corroboration/  Step 06 output
results/07_authority_report/  Step 07 output
```

This repository includes lightweight and medium-sized derived result tables required to inspect the manuscript workflow and verify the reported origin-tracing conclusions.

## Tracked result folders

The repository is intended to track the main derived result folders from Steps 01 to 07, except for individual files that are too large for normal GitHub upload.

These folders contain the QC summaries, SNP-panel tables, mutual-information diagnostics, hierarchical assignments, leave-one-out validation summaries, Random Forest corroboration summaries, and final authority-facing report.

## Excluded result folders and files

The folder:

```text
results/00_fasta/
```

is not tracked in GitHub because it contains generated FASTA files produced by Step 00 from the PHYLIP ortholog alignment inputs in:

```text
data/000_input_data/phy/
```

These FASTA files are reproducible and are therefore excluded from version control.

The file:

```text
results/01_qc/sample_missingness.tsv
```

is also excluded because it is a large generated QC matrix. Summary QC outputs are retained in the repository where possible.

## Recreating excluded outputs

To rerun the full workflow from Step 00, users need access to the controlled data package containing:

```text
data/000_input_data/phy/
```

Once that folder is restored, run:

```bash
Rscript run_public_pipeline.R
```

This will regenerate `results/00_fasta/` and downstream outputs.

## Final report outputs

The main final report files are:

```text
results/07_authority_report/FINAL_origin_tracing_authority_REPORT.xlsx
results/07_authority_report/FINAL_origin_tracing_authority_simplified.tsv
results/07_authority_report/FINAL_origin_tracing_authority_extended.tsv
```

These files summarize the final conservative origin-tracing interpretation, including macroregion assignment, conditional subregion assignment, leave-one-out validation context, Random Forest corroboration, and authority-facing reporting logic.
