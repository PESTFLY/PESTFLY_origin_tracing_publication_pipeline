# Manuscript pipeline description

This public workflow reproduces the genome-wide SNP-based origin-tracing analyses for intercepted/trapped *Bactrocera dorsalis* samples.

## Conceptual design

The workflow uses nuclear SNPs extracted from ortholog alignments. SNPs are ranked into diagnostic panels and used in a hierarchical origin-tracing framework:

1. macroregion assignment: Africa vs Asia;
2. conditional subregion assignment within the macroregion supported by the first step;
3. full-panel leave-one-out validation;
4. independent Random Forest corroboration;
5. authority-facing reporting.

## Multi-K stability

Step 04 evaluates top-K subsets of ranked SNPs. Small K values are retained in raw output, but only K values with enough usable/non-missing SNPs contribute to stability decisions.

A call is accepted only when posterior support, posterior gap, sufficient usable K coverage, global K agreement, and large-K tail agreement all pass the thresholds defined in the script parameters.

## Conservative reporting

Step 07 reports subregion only when Step 04 supports a High/Moderate subregion call. Otherwise it falls back to macroregion if macroregion is High/Moderate. If macroregion is uncertain, origin is reported as uncertain.

Random Forest support is corroborative only and never overrides Step 04.

## Downstream diagnostic SNP-panel development

The reduced SNP-panel/lab-test development workflow is intentionally excluded from this manuscript release. It uses the ranked SNP resources generated here but answers a different question: how few common, assayable SNPs are needed for simplified laboratory testing.
