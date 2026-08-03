# Pseudobulk DE One-Million-Cell Run Technical Acceptance

## Scope

This record covers execution and technical acceptance of the approved
donor-aware, cell-type-specific pseudobulk differential-expression stage on the
corrected GSE164378 profile. It does not approve genes, pathways, biological
interpretations, or public-facing claims.

The source H5AD contains 1,916,367 cells. This record applies only to the
accepted 1,000,000-cell checkpoint and must not be interpreted as technical
acceptance of an uncapped source-data run.

## Execution

- Date: 2026-07-31
- Configuration: `config/gse164378_corrected.yaml`
- Input checkpoint: `results/gse164378_full/06_annotated.h5ad` (read-only)
- Replicate mapping:
  `results/gse164378_full_replicate_corrected/tables/cell_to_biological_replicate.csv`
- Runtime: R 4.6.1, Bioconductor 3.23, edgeR 4.10.1 in
  `immune-aging-edger:bioc-3.23-edger-4.10.1`
- Output root:
  `results/gse164378_full_replicate_corrected/pseudobulk_de/`

The heavy steps were run as narrow direct module commands because Snakemake DAG
construction repeatedly stalled on this Windows environment. The aggregation
and edgeR logs record the commands, working directory, timestamps, and exit
status. The workflow still exposes the corresponding targeted rules.

## Observed Run Inventory

| Item | Observed |
| --- | ---: |
| Input cells | 1,000,000 |
| Included cells | 949,518 |
| Excluded cells | 50,482 |
| Biological replicates | 317 |
| Pseudobulk profiles | 2,944 |
| Input genes | 36,601 |
| Aggregated counts | 2,861,289,479 |
| Primary cell types | 12 |
| Exploratory cell types | 1 (NK cells) |
| Gene-by-cell-type tests | 107,948 |

The sparse aggregation audit reported exact selected-count conservation and no
dense cell-by-gene conversion. All 13 approved populations completed with
full-rank designs and at least the configured residual degrees of freedom. All
used `~ sex + batch + age_decade`; no age-only fallback was used.

## Automated Acceptance

Command:

```powershell
python -u -m src.validate_pseudobulk_de `
  --config config/gse164378_corrected.yaml `
  --outdir results/gse164378_full_replicate_corrected/pseudobulk_de `
  --report results/gse164378_full_replicate_corrected/pseudobulk_de/technical_validation.json
```

Result: passed. The validator checked required artifacts, aggregation and
profile contracts, approved analysis tiers, eligibility, unique result keys,
finite statistics, model formulas and diagnostics, per-cell-type tables,
portable manifest paths, pinned runtime versions, absence of incomplete files,
UTF-8 execution logs with successful completion markers, and independent
recalculation of within-cell-type and global BH FDR.

For review triage, the machine-readable report counts 8,408 tests below 0.05
within cell type and 7,970 below 0.05 under the global correction. These counts
are not evidence that every row is biologically credible, independent, or
suitable for presentation. They require effect-size review, diagnostic review,
annotation checks, and a separate human scientific interpretation stage.

## Corrections During Acceptance

The initial cell-type manifest stored Docker-internal `/work/...` paths. The R
writer was corrected to emit `by_cell_type/...` paths, and the existing
manifest was normalized without rerunning edgeR. The validator confirms that
all 13 relative paths resolve to non-empty tables with matching cell types and
row counts.

## Subsequent Robustness Review

- All 13 pages of `edgeR_diagnostic_plots.pdf` were rendered and inspected on
  2026-08-02. The review found no empty, clipped, or malformed page and retained
  follow-up flags for low-abundance effects in MAIT, Memory B, Naive B, and the
  exploratory NK population.
- The pre-specified robustness audit completed 212 sensitivity models and
  recorded three non-estimable scenarios without substituting another model.
  See `docs/pseudobulk_robustness_audit.md`.

## Remaining Review

- Review candidate effect sizes, profile support, gene annotation, and known
  technical artifacts before interpreting statistical significance.
- Define and approve a human-review protocol before ranking genes, running
  pathway enrichment, producing public figures, or changing the README results
  narrative.
- Investigate the persistent Snakemake DAG-construction delay separately; a
  60-second target-listing attempt timed out and left no orphaned process.
