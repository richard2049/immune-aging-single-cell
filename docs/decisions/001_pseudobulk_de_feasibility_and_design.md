# Decision Card: Donor-Aware Pseudobulk DE Feasibility And Design

## Metadata

- Stage: Donor-aware pseudobulk DE design and implementation
- Priority: High
- Status: Implemented - full-data execution and technical validation complete;
  biological interpretation pending
- Risk classification: C - scientifically important

## Context

Observed:

- `IMMUNE-AGING_NOTES.md`, `README.md`, and
  `docs/phase2_next_steps_checklist.md` identify donor-aware,
  cell-type-specific differential expression as the main remaining scientific
  gap.
- The workflow already produces an annotated checkpoint at
  `results/gse164378_full/06_annotated.h5ad`.
- The initial donor-composition table contains 166 donor labels, 84 cell-type
  labels, and 14 batches. The feasibility audit established that those donor
  labels are not globally unique biological replicates.
- In the full profile, QC is disabled and scVI reads `adata.X` directly rather
  than a dedicated `counts` layer.

This stage is scientifically important because the experimental unit, count
source, cell-type eligibility rules, model design, and multiple-testing policy
can materially change the biological conclusions. Human approval is required
before implementation.

## Decision

The bounded feasibility stage is implemented in
`src/pseudobulk_de_feasibility.py`. On 2026-07-30, the user approved the
following implementation contract:

1. Sum raw counts over all technical libraries within each
   `biological_replicate_id x cell_type` group. Never treat cells or technical
   libraries as independent age observations.
2. Use continuous age in decades. Report the primary age coefficient as log2
   fold change per 10 years.
3. Retain the feasibility thresholds of at least 50 cells per pseudobulk, 12
   qualifying biological replicates, 12 years of age coverage, and five
   residual model degrees of freedom.
4. Treat the 12 well-supported populations as primary. Analyze NK cells using
   the same method, but label them exploratory because only 18 replicates and
   eight residual degrees of freedom qualify.
5. Use the primary design `~ sex + batch + age_decade`. Require a full-rank
   design and adequate residual degrees of freedom. If it is not estimable,
   exclude that population from primary adjusted inference and report an
   age-only sensitivity result rather than silently dropping covariates.
6. Use `edgeR` quasi-likelihood with TMM normalization and design-aware
   `filterByExpr` in a dedicated, version-pinned R/Bioconductor environment.
7. Report BH FDR within each cell type and global BH FDR across all tested
   gene-cell-type combinations. Reserve cross-cell-type or headline findings
   for globally corrected results.
8. Limit this implementation stage to pseudobulk counts, metadata, eligibility
   records, design diagnostics, gene-level result tables, and technical
   diagnostic plots. Pathway interpretation and public gene-level claims
   require a later D-class review.

## Feasibility Outcome

- The complete 1,000,000 x 36,601 CSR matrix passed raw-count integrity checks:
  no non-finite, negative, fractional, or out-of-range stored values.
- `donor_id` contains 166 reused labels and is not a valid global replicate
  key. It produces age conflicts for 102 labels.
- `sample_id` is a multiplexed technical-library identifier.
- Source `Tube_id` resolves 317 biological replicates with complete coverage
  and no age, sex, batch, or donor-label conflicts.
- Thirteen cell types pass the candidate support and design rules. Twelve have
  at least 105 qualifying replicates; NK cells are borderline with 18
  qualifying replicates and eight residual degrees of freedom.
- Sparse aggregation is feasible without a dense cell-by-gene conversion. The
  projected dense int64 pseudobulk matrix would be about 862 MB, but chunked
  sparse aggregation remains the required strategy.
- `Rscript` and `edgeR` are not available in the current environment.

Detailed results and the proposed contract are in
`docs/pseudobulk_de_design.md`. Generated audit outputs are under
`results/gse164378_full/pseudobulk_de_feasibility/`.

The replicate remediation and controlled downstream regeneration are complete.
The corrected analyses expose source `Tube_id` as
`biological_replicate_id`; technical validation found 1,000,000 mapped cells,
317 biological replicates, no metadata conflicts, and no grouped
cross-validation leakage. Public-result review was completed separately before
promotion.

The preferred statistical direction for review is summed raw counts per
biological replicate x cell type followed by a replicate-aware count model.
`edgeR` quasi-likelihood was approved because it directly supports
pseudobulk count analysis, filtering, normalization, design matrices, and
dispersion estimation. It will use a dedicated R/Bioconductor environment so
that the Python/scVI environment remains isolated from the count-model
runtime.

## Approval Gates

The design gates below are approved:

- the raw-count source and integrity checks pass (completed);
- the biological replicate and aggregation unit are unambiguous and propagated
  through the corrected workflow (completed);
- eligible cell types and exclusion reasons are approved (completed);
- the age model and estimable covariates are fixed (completed);
- the gene-filtering and multiple-testing policies are fixed (completed);
- the engine and environment strategy are approved;
- the full-scale memory plan does not require dense matrix conversion
  (completed).

The dedicated `edgeR` runtime passed a bounded qualification test on
2026-07-30 using R 4.6.1, Bioconductor 3.23, and edgeR 4.10.1. The full
GSE164378 analysis then completed on 2026-07-31 with the same pinned runtime.

No gene-level biological claim should be made during this design stage.

## Implementation Outcome

- `src/pseudobulk_aggregate.py` performs sparse, chunked aggregation and checks
  exact count conservation.
- `src/pseudobulk_edger.R` implements the approved adjusted model, age-only
  fallback, TMM normalization, design-aware filtering, quasi-likelihood test,
  within-cell-type FDR, and global FDR.
- `src/run_edger.py` runs the model through a time-bounded, logged Docker or
  native-R interface.
- `containers/edger/Dockerfile` pins the Bioconductor base-image digest and
  asserts the approved R, Bioconductor, and edgeR versions.
- `workflows/Snakefile.replicate_corrected` exposes targeted aggregation and DE
  rules without adding the heavy stage to the default target.
- `src/validate_pseudobulk_de.py` verifies the approved analysis contract,
  output schemas, model diagnostics, portable manifest paths, package versions,
  and independent recalculation of both BH corrections.
- `tests/test_pseudobulk_de.py` covers count conservation, input rejection,
  file interchange, runtime command construction, and an optional Docker
  integration test.

The full run aggregated 949,518 eligible cells into 2,944 profiles from 317
biological replicates and retained 36,601 input genes without a dense
cell-by-gene conversion. edgeR completed 107,948 gene-by-cell-type tests across
12 primary populations and exploratory NK cells. All 13 designs were full rank,
met the residual-df threshold, and used the approved adjusted formula; no
age-only fallback was used. The generated
`technical_validation.json` passed all automated checks.

These outcomes establish technical completion, not biological validity. The
gene-level tables, FDR counts, pathway analysis, and any public-facing claim
remain subject to a separate D-class interpretation and review stage.

## Rationale And Tradeoffs

Pseudobulk aggregation preserves the donor as the independent biological
replicate and avoids treating thousands of cells from one donor as independent
age observations. This closes a scientific gap that predefined signature
scores cannot address.

The tradeoff is that the stage adds method and environment complexity. Rare
cell types will lose power or be excluded, adjustment for sex and batch may be
limited by the observed design, and an R/Bioconductor engine would expand the
runtime contract. Those costs are preferable to producing anti-conservative
per-cell differential-expression results.

Alternatives considered:

- Naive per-cell differential expression: rejected because it does not respect
  donor-level replication.
- Continue with signature scores only: lower effort, but it does not provide
  genome-wide discovery.
- Begin external Phase 2 validation first: premature while the primary
  within-dataset transcriptional question remains incomplete.
- Productize data acquisition first: valuable engineering work, but lower
  immediate scientific priority than resolving the documented DE gap.

Supporting methodological sources:

- [Squair et al., Confronting false discoveries in single-cell differential
  expression](https://www.nature.com/articles/s41467-021-25960-2)
- [edgeR User's Guide, single-cell differential expression with
  pseudo-bulking](https://master.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)

## Interview Question

- Question: Why should age-associated expression in this single-cell dataset
  be tested with donor-level pseudobulk profiles rather than by treating cells
  as independent replicates?
- Skill tested: Experimental-unit identification, pseudoreplication,
  pseudobulk aggregation, and count-based differential-expression design.
- Expected answer should mention: Age is measured per donor; cells from one
  donor are correlated observations rather than independent replicates; raw
  counts should be summed within an eligible cell type and biological
  replicate; the model should estimate between-donor variability, check
  covariate identifiability, and control multiple testing.
- Status: Ready

## Follow-up Biological Question

- Question: Which genes and pathways show donor-aware age associations within
  sufficiently supported immune cell populations, beyond composition shifts
  and predefined signature scores?
- Biological theme: Cell-type-specific transcriptional immune aging.
- Why it matters: It addresses genome-wide within-cell-type changes that the
  current composition and signature analyses cannot establish.
- Possible next analysis: After the technical implementation and result review,
  perform direction-aware pathway enrichment on ranked gene-level effects.
- Status: Open

## No-AI Survival Exercise

- Estimated time: 20-30 minutes
- Skill tested: Identifying the biological replicate and evaluating whether a
  cell type has adequate donor-level support for pseudobulk analysis.
- Allowed resources: This decision card,
  `results/gse164378_full/pseudobulk_de_feasibility/replicate_celltype_support.csv`,
  a spreadsheet or short local script, and package documentation; no
  generative AI.
- Expected output: A compact table showing, for each cell type, the number of
  biological replicates with at least 50 cells and the observed age span, plus
  a short explanation of why `Tube_id`, rather than `donor_id` or `sample_id`,
  is the proposed experimental unit.
- Pass criteria: Correctly treats donors rather than cells as independent
  replicates, identifies unsupported cell types without interpreting absence
  as biology, and states that summed raw counts are required for the planned
  count model.
- Follow-up biological question: Do the eligible cell types cover the immune
  populations needed to distinguish broad aging effects from lineage-specific
  transcriptional changes?
- Status: Not started
