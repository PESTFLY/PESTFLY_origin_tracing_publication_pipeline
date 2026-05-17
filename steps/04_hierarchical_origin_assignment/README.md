# Step 04 — Hierarchical origin assignment

Performs the primary hierarchical origin assignment using Step 02 diagnostic SNP panels.

Assignment hierarchy:

1. `P1`: macroregion assignment, Africa vs Asia.
2. `P2`: conditional African subregion assignment, only if P1 supports Africa.
3. `P3`: conditional Asian subregion assignment, only if P1 supports Asia.

The step evaluates multiple top-K SNP subsets and accepts an assignment only when posterior support, posterior gap, global K agreement, and large-K tail stability pass the reporting thresholds.

## Inputs

```text
results/01_qc/metadata_clean.tsv
results/02_snp_panels/panels_index.tsv
results/02_snp_panels/panels/<panel_id>/snp_map.tsv
results/02_snp_panels/panels/<panel_id>/snp_matrix.rds
```

## Main outputs

```text
results/04_origin_assignment/final_assignment.tsv
results/04_origin_assignment/P1_macroregion_raw.tsv
results/04_origin_assignment/P1_macroregion_stability.tsv
results/04_origin_assignment/P2_africa_subregion_raw.tsv
results/04_origin_assignment/P3_asia_subregion_raw.tsv
```

## Publication patch

The included script explicitly restricts references by panel:

- P1: Africa and Asia references;
- P2: African references only;
- P3: Asian references only.

## Run

```bash
Rscript steps/04_hierarchical_origin_assignment/run.R
```
