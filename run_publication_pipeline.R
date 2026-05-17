    #!/usr/bin/env Rscript

    # Run the public Molecular Ecology Resources origin-tracing pipeline.
    # Execute from the repository root:
    #
    #   Rscript run_publication_pipeline.R

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

    for (script in steps) {
      message("
============================================================")
      message("Running: ", script)
      message("============================================================")

      status <- system2("Rscript", script)

      if (!identical(status, 0L)) {
        stop("Pipeline stopped because this step failed: ", script)
      }
    }

    message("
Public origin-tracing pipeline completed successfully.")
    message("Final report folder: results/07_authority_report")
