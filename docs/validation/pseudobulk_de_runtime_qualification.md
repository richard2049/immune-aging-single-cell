# Pseudobulk DE Runtime Qualification

> [!IMPORTANT]
> This record qualifies the superseded fixed-effect edgeR runtime only. The
> maintained longitudinal workflow uses the separately pinned `dream` image,
> which still requires bounded qualification.

## Status

Passed on 2026-07-30 for bounded technical validation. The full GSE164378
analysis was not run.

## Runtime

- Base image:
  `bioconductor/bioconductor_docker:RELEASE_3_23`
- Base-image digest:
  `sha256:1d871e1ca9cca76b220eb16e22677e728f4352f81a9ee91aaf29e24aea43e624`
- Local image:
  `immune-aging-edger:bioc-3.23-edger-4.10.1`
- Built image ID:
  `sha256:dd1c123f48bd346fd025c13bbc23ac8cd88ae8ee33b55682b42a3c3bacf120a6`
- R: 4.6.1
- Bioconductor: 3.23
- edgeR: 4.10.1

The Docker build asserts the approved versions and fails if the runtime drifts.

## Bounded Fixture

- One synthetic cell type.
- 12 independent pseudobulk profiles.
- 40 genes with integer counts.
- Sex and three-level batch covariates.
- Continuous age spanning 55 years.
- Primary design: `~ sex + batch + age_decade`.

The test completed design-aware filtering, TMM normalization, dispersion
estimation, quasi-likelihood fitting, the age contrast, within-cell-type FDR,
global FDR, result-table generation, runtime reporting, and a diagnostic PDF.
It returned 40 gene-level tests under `adjusted_primary`.

## Additional Validation

- The ordinary suite completed 19 tests: 18 passed and the Docker integration
  test was skipped unless explicitly enabled.
- The explicit Docker integration test passed.
- Sparse aggregation tests confirmed exact profile and total-count
  conservation.
- Dense and fractional input matrices were rejected.
- Matrix Market output was read back and checked against profile metadata.

A Snakemake dry-run timed out during DAG construction on Windows. Process and
lock checks found no orphaned Python/Snakemake process and no remaining lock.
The targeted rules are therefore implemented, but workflow-engine DAG
validation remains an execution-environment limitation to recheck before the
full run.

## Interpretation Boundary

This qualification validates software behavior and runtime compatibility. It
does not validate any immune-aging gene association, pathway, causal claim, or
clinical interpretation.
