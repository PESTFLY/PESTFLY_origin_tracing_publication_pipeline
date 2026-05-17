# Data directory

Expected input structure:

```text
data/000_input_data/
├── metadata.xlsx
└── phy/
    └── *.phy
```

The workflow expects:

- a metadata spreadsheet at `data/000_input_data/metadata.xlsx`;
- OMA/read2tree PHYLIP ortholog alignments under `data/000_input_data/phy/`.

Large raw data files are intentionally not included in the public code package unless explicitly permitted for redistribution.

The Step 00 script writes converted FASTA files to `results/00_fasta/` by default. Subsequent steps use these validated default paths.
