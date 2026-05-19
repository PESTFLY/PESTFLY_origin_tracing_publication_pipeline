# PESTFLY origin-tracing pipeline for *Bactrocera dorsalis* interceptions

This repository contains the publication-ready workflow used for genome-wide nuclear SNP-based origin tracing of intercepted/trapped *Bactrocera dorsalis* samples in the Mol Ecol Res manuscript.

The workflow stops at Step 07, which produces the authority-facing report. Reduced diagnostic SNP-panel development is **not** part of this public manuscript pipeline and is maintained as a separate internal extension.

## Pipeline overview

| Step | Folder | Purpose | Main result folder |
|---:|---|---|---|
| 00 | `steps/00_convert_phylip_to_fasta/` | Convert OMA/read2tree PHYLIP ortholog alignments to FASTA | `results/00_fasta/` |
| 01 | `steps/01_metadata_snapshot_and_ortholog_qc/` | Clean metadata, snapshot sample labels, and perform ortholog QC | `results/01_qc/` |
| 02 | `steps/02_snp_panel_discovery_and_ranking/` | Extract and rank diagnostic nuclear SNP panels | `results/02_snp_panels/` |
| 03 | `steps/03_mutual_information_snp_diagnostics/` | Add supportive mutual-information SNP diagnostics | `results/03_mi_diagnostics/` |
| 04 | `steps/04_hierarchical_origin_assignment/` | Perform hierarchical macroregion-first origin assignment | `results/04_origin_assignment/` |
| 05 | `steps/05_full_panel_loo_validation/` | Validate the full-panel assignment logic by leave-one-out | `results/05_loo_validation/` |
| 06 | `steps/06_random_forest_corroboration/` | Run Random Forest as independent corroboration | `results/06_rf_corroboration/` |
| 07 | `steps/07_authority_facing_report/` | Produce the final authority-facing report | `results/07_authority_report/` |

## Quick start

Run the full workflow from the repository root:

```bash
Rscript run_public_pipeline.R
```

Or run individual steps manually:

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

## Inputs

The expected public metadata file is:

```text
data/000_input_data/metadata.xlsx
```

The raw PHYLIP alignment folder is expected at:

```text
data/000_input_data/phy/
```

The `.gitignore` is configured to exclude the `phy/` folder and generated Step 00 FASTA outputs by default. If large alignments are deposited separately, provide the corresponding DOI or download instructions in `data/README_data.md`.

## Main final output

```text
results/07_authority_report/FINAL_origin_tracing_authority_REPORT.xlsx
results/07_authority_report/FINAL_origin_tracing_authority_simplified.tsv
results/07_authority_report/FINAL_origin_tracing_authority_extended.tsv
```

## Notes on Step 03

Step 03 has its own result folder (`results/03_mi_diagnostics/`) to keep the public pipeline numbering one-to-one. It also writes per-panel `snp_mi.tsv` files into `results/02_snp_panels/panels/<panel_id>/` for compatibility with downstream diagnostic-panel development.
