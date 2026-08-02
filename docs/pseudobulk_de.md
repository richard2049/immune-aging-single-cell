# Donor-Aware Pseudobulk Differential Expression

## Scope

This targeted stage tests age-associated expression within approved immune
cell populations while treating `biological_replicate_id` as the experimental
unit. It does not treat cells or technical libraries as independent
replicates.

The implementation is technically validated on both a bounded fixture and the
full GSE164378 corrected-replicate checkpoint. The full gene-level results have
not yet undergone biological interpretation, so this document defines an
execution and technical acceptance contract rather than a biological
conclusion.

## Approved Contract

- Sum raw counts within each `biological_replicate_id x cell_type` profile.
- Require at least 50 cells per profile, 12 qualifying replicates, 12 years of
  age coverage, and five residual model degrees of freedom.
- Analyze 12 well-supported populations as primary and NK cells as
  exploratory.
- Model continuous age in decades with `~ sex + batch + age_decade`.
- Use an age-only sensitivity model only when the adjusted design is not
  estimable; never drop covariates silently.
- Use edgeR quasi-likelihood, TMM normalization, and design-aware
  `filterByExpr`.
- Report BH FDR within each cell type and globally across all tested
  gene-cell-type combinations.

## Runtime

The Dockerfile pins the Bioconductor 3.23 base-image digest and verifies:

- R 4.6.x;
- Bioconductor 3.23;
- edgeR 4.10.1.

Build the image from PowerShell:

```powershell
docker --context desktop-linux build --progress plain `
  -t immune-aging-edger:bioc-3.23-edger-4.10.1 `
  -f containers/edger/Dockerfile .
```

Run the optional bounded integration test:

```powershell
$env:RUN_EDGER_DOCKER_TESTS = "1"
python -m unittest `
  tests.test_pseudobulk_de.PseudobulkDifferentialExpressionTests.test_pinned_edger_runtime_on_bounded_fixture `
  -v
```

## Full-Data Target

The stage is intentionally absent from the default workflow target. Run it
explicitly:

```powershell
python -m snakemake `
  -s workflows/Snakefile.replicate_corrected `
  -c 1 pseudobulk_de_corrected `
  --configfile config/config.real.full.replicate_corrected.yml `
  --printshellcmds
```

Validate existing outputs without rerunning aggregation or edgeR:

```powershell
python -u -m src.validate_pseudobulk_de `
  --config config/config.real.full.replicate_corrected.yml `
  --outdir results/gse164378_full_replicate_corrected/pseudobulk_de `
  --report results/gse164378_full_replicate_corrected/pseudobulk_de/technical_validation.json
```

The aggregation command prints cell progress. The edgeR wrapper streams output
to the terminal, writes `edgeR.log`, records the working directory and command,
and enforces the configured timeout.

## Output Contract

Outputs are written under
`results/gse164378_full_replicate_corrected/pseudobulk_de/`:

- `pseudobulk_counts.mtx.gz`: sparse profile-by-gene raw-count matrix;
- `pseudobulk_profiles.csv`: profile metadata and library sizes;
- `pseudobulk_genes.csv`: ordered gene identifiers;
- `celltype_eligibility.csv`: inclusion status and explicit reasons;
- `aggregation_audit.json`: count-conservation and provenance checks;
- `gene_level_results.csv.gz`: combined gene-level results;
- `by_cell_type/*.csv.gz`: cell-type-specific result tables;
- `celltype_result_manifest.csv`: per-cell-type result index;
- `design_diagnostics.csv`: support, rank, residual df, and model status;
- `edgeR_diagnostic_plots.pdf`: technical mean-difference plots;
- `runtime_versions.csv`, `session_info.txt`, and `edgeR.log`: runtime
  provenance.
- `technical_validation.json`: machine-readable acceptance checks and a compact
  run inventory, explicitly separated from biological interpretation.

These outputs remain analytical results pending review. Pathway enrichment,
public gene-level figures, and biological claims require a later D-class
interpretation stage.

## Full-Data Execution Status

The full corrected-profile run completed on 2026-07-31. It produced 2,944
pseudobulk profiles from 317 biological replicates and 107,948 tests across 13
cell types. All populations used the approved adjusted model, and the technical
validator passed. See
`docs/validation/pseudobulk_de_full_run_acceptance.md` for the acceptance record
and remaining manual review requirements.

## Robustness Audit And Evidence Preparation

Decision card 005 defines a separate C-class stage that reuses the validated
pseudobulk counts and primary result table. It evaluates a higher profile-cell
floor, covariate-omission diagnostics, and leave-one-batch-out models. These
models do not replace the adjusted primary analysis and are not independent
replications.

The full audit completed on 2026-08-02. It preserved 7,970 global-FDR
candidate rows, recorded 212 completed and three explicitly non-estimable
cell-type-by-scenario models, and produced a bounded 249-row manual-review
queue. No model failed. See `docs/pseudobulk_robustness_audit.md` for the
technical record and interpretation boundary.

Run the sensitivity models explicitly:

```powershell
python -u -m src.run_pseudobulk_robustness `
  --config config/config.real.full.replicate_corrected.yml `
  --matrix results/gse164378_full_replicate_corrected/pseudobulk_de/pseudobulk_counts.mtx.gz `
  --profiles results/gse164378_full_replicate_corrected/pseudobulk_de/pseudobulk_profiles.csv `
  --genes results/gse164378_full_replicate_corrected/pseudobulk_de/pseudobulk_genes.csv `
  --primary-results results/gse164378_full_replicate_corrected/pseudobulk_de/gene_level_results.csv.gz `
  --sensitivity-out results/gse164378_full_replicate_corrected/pseudobulk_de/robustness/sensitivity_gene_results.csv.gz `
  --diagnostics-out results/gse164378_full_replicate_corrected/pseudobulk_de/robustness/scenario_diagnostics.csv `
  --manifest-out results/gse164378_full_replicate_corrected/pseudobulk_de/robustness/scenario_manifest.csv `
  --log-out results/gse164378_full_replicate_corrected/pseudobulk_de/robustness/edgeR_robustness.log
```

After recording the diagnostic-plot review, prepare the evidence tables:

```powershell
python -m src.prepare_pseudobulk_evidence `
  --config config/config.real.full.replicate_corrected.yml `
  --pseudobulk-dir results/gse164378_full_replicate_corrected/pseudobulk_de `
  --report-out docs/pseudobulk_robustness_audit.md
```

The compact `review_queue.csv` is a workload-management artifact. Selection
for that table does not approve a gene, pathway, or biological claim. Any such
disposition requires a later D-class human review.
