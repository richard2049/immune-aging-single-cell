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

## Interpretation

- Association with chronological age is not evidence of causal aging
  mechanisms.
- Cell-type composition, gene-expression signatures, pseudobulk expression,
  and latent embeddings are distinct evidence layers.
- Automated cell-type annotations are not ground truth and require biological
  review.
- Technical acceptance is separate from biological interpretation and public
  promotion.
