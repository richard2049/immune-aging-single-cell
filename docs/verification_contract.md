# Immune Aging Verification Contract

## Biological Units

- Biological replicates are the independent units for replicate-level age
  associations and age prediction.
- Cells from one biological replicate are not independent replicates.
- `biological_replicate_id` must come from an authoritative metadata field.
  For GSE164378, source `Tube_id` is canonical; reused `donor_id` labels are
  retained only as source metadata.
- Cell-level operations may be used for representation learning, QC,
  clustering, or annotation, but inferential claims must respect replicate
  structure.

## Metadata

- Donor, biological replicate, age, sex, batch, cohort, sample, tissue, and
  cell-type fields must have documented provenance.
- Joins must validate cardinality, preserve cell count and order, and report
  unmatched identifiers.
- Exact duplicate metadata rows may be collapsed. A join key with conflicting
  payloads is invalid and must stop execution.
- Parsing sample names is not sufficient evidence when authoritative metadata
  exist.
- Missing values must not be silently imputed. Configured model covariates must
  not be silently removed.
- Explicitly configured adjustment covariates must be present and complete, and
  non-estimable nuisance designs must stop inferential analysis rather than
  trigger an unadjusted fallback.

## Composition

- Eligible biological replicates must be crossed with the observed cell-type
  universe before fractions are calculated. Absence of a population is a zero
  count, not a missing replicate.
- Report total independent replicates separately from replicates with a
  positive count for the population.
- Fractions over the complete cell-type grid must sum to one within each
  eligible biological replicate, subject only to numerical tolerance.

## Cross-Validation And Prediction

- No biological replicate may occur in more than one fold.
- Model candidates must use identical grouped folds for comparison.
- Feature selection, scaling, imputation, prediction-specific dimensionality
  reduction, and hyperparameter selection must use training data only.
- Requested model dependencies are part of the configuration contract; a
  missing library must fail explicitly rather than alter the candidate set.
- Report replicate-level predictions, fold variability, uncertainty, and a
  simple baseline.
- The current regressor uses grouped cross-validation on `X_scVI`, but
  `X_scVI` was learned once from the full cohort. Its metrics are internal,
  transductive estimates, not fully end-to-end inductive performance for unseen
  donors.

## Single-Cell Constraints

- Large sparse matrices must not be densified without a documented and bounded
  memory justification.
- Differential or association analyses must account for biological
  replicates; cell counts must not inflate the apparent sample size.
- Outputs generated from demo fixtures are not biological evidence.
- Stochastic stages must use the configured seed and record it with outputs.

## External Resources And Environments

- The runtime must contain every dependency requested by tracked configs.
- The CellTypist output must record the package version, model identifier,
  resolved model path, model metadata when available, and SHA-256 checksum.
- A configured model checksum, when provided, must be verified before
  annotation.
- Environment changes require a clean environment solve and import smoke test
  before broad workflow execution.

## Checkpoint And Failure Contracts

- Maintained real-data profiles must fail on missing analytical inputs,
  dependencies, invalid designs, or incomplete provenance. Placeholder outputs
  are allowed only when a profile opts in explicitly for demo use.
- Reused checkpoints must record the resolved configuration, source-file
  identity, required data representations, and stage provenance. A successful
  prior run does not establish compatibility with current code or config.
- QC reports must account for stage-wise cell retention and expose metadata and
  annotation uncertainty by relevant strata without imposing universal
  biological thresholds.

## Interpretation

- Association with chronological age is not evidence of causal aging
  mechanisms.
- Cell-type composition, gene-expression signatures, pseudobulk expression,
  and latent embeddings are distinct evidence layers.
- Automated cell-type annotations are not ground truth and require biological
  review.
- Technical acceptance is separate from biological interpretation and public
  promotion.
