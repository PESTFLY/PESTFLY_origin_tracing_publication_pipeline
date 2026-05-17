# PESTFLY origin-tracing pipeline for Molecular Ecology Resources

This repository contains the publication-ready, genome-wide nuclear SNP origin-tracing workflow used for the manuscript. It is intentionally limited to the manuscript/authority-reporting analysis and stops at Step 07.

## Pipeline steps and result folders

| Step | Folder | Purpose | Default output / key input for next step |
|---|---|---|---|
| 00 | `steps/00_convert_phylip_to_fasta/` | Convert OMA/read2tree PHYLIP ortholog alignments to FASTA | `results/00_fasta/` |
| 01 | `steps/01_metadata_snapshot_and_ortholog_qc/` | Clean metadata, harmonize sample labels, summarize ortholog QC | `results/01_qc/` |
| 02 | `steps/02_snp_panel_discovery_and_ranking/` | Extract and rank diagnostic nuclear SNPs for origin-tracing panels | `results/02_snp_panels/` |
| 03 | `steps/03_mutual_information_snp_diagnostics/` | Add mutual-information SNP diagnostics to the Step 02 panel resources | `results/02_snp_panels/` |
| 04 | `steps/04_hierarchical_origin_assignment/` | Hierarchical macroregion-first origin assignment | `results/04_origin_assignment/` |
| 05 | `steps/05_full_panel_loo_validation/` | Leave-one-out validation of the full-panel assignment logic | `results/05_loo_validation/` |
| 06 | `steps/06_random_forest_corroboration/` | Random Forest corroboration of Step 04 assignments | `results/06_rf_corroboration/` |
| 07 | `steps/07_authority_facing_report/` | Integrate results into the authority-facing report | `results/07_authority_report/` |

Step 03 deliberately writes into `results/02_snp_panels/` because it annotates the SNP-panel resources created by Step 02 rather than creating an independent downstream dataset.

## Run

From the repository root:

```bash
Rscript run_publication_pipeline.R
```

Or run individual steps:

```bash
Rscript steps/00_convert_phylip_to_fasta/run.R
Rscript steps/01_metadata_snapshot_and_ortholog_qc/run.R
Rscript steps/02_snp_panel_discovery_and_ranking/run.R
Rscript steps/03_mutual_information_snp_diagnostics/run.R
Rscript steps/04_hierarchical_origin_assignment/run.R
Rscript steps/05_full_panel_loo_validation/run.R
Rscript steps/06_random_forest_corroboration/run.R
Rscript steps/07_authority_facing_report/run.R
```

## Main manuscript outputs

The main authority-reporting files are:

```text
results/07_authority_report/FINAL_origin_tracing_authority_simplified.tsv
results/07_authority_report/FINAL_origin_tracing_authority_extended.tsv
results/07_authority_report/FINAL_origin_tracing_authority_REPORT.xlsx
```

The primary origin-assignment table used by downstream steps is:

```text
results/04_origin_assignment/final_assignment.tsv
```

## Scope

Reduced diagnostic SNP-panel development is not part of this publication release. It is handled separately in the internal diagnostic SNP-development pipeline, where the published Steps 00-07 are repeated unchanged and the workflow continues with Steps 08-09.
