# Step 02 — Diagnostic SNP-panel discovery and ranking

Extracts biallelic nuclear SNPs from ortholog alignments and constructs ranked diagnostic panels for macroregion and subregion assignment.

Default panels:

- `P1_macroregion_africa_vs_asia`
- `P2_subregion_within_africa`
- `P3_subregion_within_asia`

## Inputs

```text
results/00_fasta/
results/01_qc/metadata_clean.tsv
```

## Main outputs

```text
results/02_snp_panels/panels_index.tsv
results/02_snp_panels/panels/<panel_id>/snp_map.tsv
results/02_snp_panels/panels/<panel_id>/snp_matrix.rds
```

## Run

```bash
Rscript steps/02_snp_panel_discovery_and_ranking/run.R
```
