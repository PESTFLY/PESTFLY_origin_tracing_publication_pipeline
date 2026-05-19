#!/usr/bin/env Rscript

# Run the complete publication pipeline from Step 00 to Step 07.
# Execute from the repository root:
#
#   Rscript run_public_pipeline.R

steps <- c(
  "steps/00_convert_phylip_to_fasta/run.R",
  "steps/01_metadata_snapshot_and_ortholog_qc/run.R",
  "steps/02_snp_panel_discovery_and_ranking/run.R",
  "steps/03_mutual_information_snp_diagnostics/run.R",
  "steps/04_hierarchical_origin_assignment/run.R",
  "steps/05_full_panel_loo_validation/run.R",
  "steps/06_random_forest_corroboration/run.R",
  "steps/07_authority_facing_report/run.R"
)

for (s in steps) {
  message("\n============================================================")
  message("Running: ", s)
  message("============================================================")

  status <- system2("Rscript", s)

  if (!identical(status, 0L)) {
    stop("Pipeline stopped because this step failed: ", s, call. = FALSE)
  }
}

message("\nPublication pipeline completed successfully.")
