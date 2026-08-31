# Longitudinal Pseudobulk Design

## Status

Implementation and inferential regeneration are technically complete. The
bounded Docker qualification, longitudinal metadata audit, repeated-measures
model, one-sample-per-subject sensitivity, diagnostic PDF, and technical output
validator pass. Human gene- and pathway-level interpretation remains pending.
This document supersedes the earlier Tube-as-independent-replicate design.

## Unit Contract

- Aggregate raw counts by `sample_unit_id x cell_type`.
- Retain `subject_id` on each profile to model repeated samples.
- Retain `technical_library_id` provenance and confirm exact count conservation.
- Never use cells or sample units as independent people.

## Primary Model

For each pre-specified cell type:

```r
~ sex + batch + age_decade + (1 | subject_id)
```

The implementation uses edgeR for count containers, design-aware filtering and
TMM normalization, followed by `voomWithDreamWeights` and `dream`. The age
coefficient is reported as log2 fold change per 10 years.

A cell-type model is excluded rather than simplified when the fixed design is
rank deficient, age is absent, age support is insufficient, residual degrees of
freedom are inadequate, or no subjects have repeated profiles. There is no
automatic age-only fallback.

## Pre-specified Sensitivity

Select one sample per subject by earliest age, breaking ties by lexical
`sample_unit_id`. Fit the same fixed effects without the random intercept using
voom/limma. Store this result separately as
`one_sample_per_subject_results.csv.gz`.

The sensitivity discards observations and depends on its selection rule. It is
a robustness diagnostic, not a substitute for the primary model or a source of
automatic biological validation.

## Multiplicity And Outputs

- Report BH FDR within each cell type and globally across all gene-cell-type
  tests, separately for primary and sensitivity analyses.
- Preserve profile metadata, genes, sparse counts, eligibility decisions,
  design diagnostics, runtime versions, session information, and logs.
- Do not interpret or promote genes/pathways during technical validation.

## Acceptance Gates

- exact sparse count conservation;
- unique sample-unit-by-cell-type profiles with one subject per sample;
- configured minimum support in both qualifying sample units and distinct
  subjects, with subject count governing independent biological support;
- full-rank fixed designs and the configured residual degrees of freedom;
- completed primary and sensitivity models for every accepted cell type;
- pinned R/Bioconductor/edgeR/variancePartition versions;
- independently reviewed diagnostic plots and technical report.

Thresholds in the tracked profile are study decisions, not universal defaults.
Changing them requires a new decision record and output namespace.
