# Immune Aging Verification Contract

## Biological Units

- `subject_id` is the dependency and cross-validation group. For the Blood Age
  Atlas it is derived from source `Donor_id`.
- `sample_unit_id` is the longitudinal observation and aggregation unit. For
  this study it is derived from source `Tube_id`.
- `technical_library_id` records the technical input and is derived from source
  `File_name`; it is not an independent biological replicate.
- Cells and repeated sample units from one subject are not independent subjects.
- Every sample unit must map to exactly one subject, age, sex, and configured
  batch. A subject may map to several sample units and ages.
- Cell-level operations may be used for representation learning, QC,
  clustering, or annotation, but inferential claims must respect replicate
  structure.

## Metadata

- Subject, sample unit, technical library, age, sex, batch, cohort, tissue, and
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

- The primary Blood Age Atlas composition denominator retains every QC-passed
  cell. Cells without an approved primary analysis label contribute only to the
  derived `Other/unresolved` denominator category; raw and exploratory labels
  remain preserved and `Other/unresolved` is not tested as a biological
  population.
- Eligible sample units must be crossed with the observed cell-type
  universe before fractions are calculated. Absence of a population is a zero
  count, not a missing sample.
- Report sample units and subjects separately, including the number with a
  positive count for the population.
- Fractions over the complete cell-type grid must sum to one within each
  eligible sample unit, subject only to numerical tolerance.
- Primary uncertainty must account for repeated samples within subjects.
- The deterministic one-sample-per-subject sensitivity uses earliest age and
  then lexical sample ID; it does not replace the primary repeated-measures model.
- Signature associations attempt exchangeable Gaussian GEE first. Only the
  specific GEE nonconvergence condition may trigger the configured independence
  working-correlation fallback; data, covariates, subject groups, family, and
  population-average estimand remain unchanged.
- The signature fallback retains robust subject-clustered covariance and must
  record its model, covariance, trigger, standard error, and manual-review
  status. Missing data, rank deficiency, non-finite inference, invalid
  configuration, or a second nonconvergence must stop rather than trigger an
  additional substitute model.
- Fallback estimates remain subject to the prespecified one-sample-per-subject
  comparison and cannot independently authorize a biological claim.

## Cross-Validation And Prediction

- No subject may occur in more than one fold. All sample units and cell-type
  rows from a subject must receive the same outer-fold assignment.
- Model candidates must use identical grouped folds for comparison.
- Feature selection, scaling, imputation, prediction-specific dimensionality
  reduction, and hyperparameter selection must use training data only.
- Requested model dependencies are part of the configuration contract; a
  missing library must fail explicitly rather than alter the candidate set.
- Report sample-unit predictions, fold variability, subject-clustered
  uncertainty, and a
  simple baseline.
- The current regressor uses grouped cross-validation on `X_scVI`, but
  `X_scVI` was learned once from the full cohort. Its metrics are internal,
  transductive estimates, not fully end-to-end inductive performance for unseen
  subjects.

## Pseudobulk

- Raw counts are summed within `sample_unit_id x cell_type`; technical-library
  contributions are retained in metadata and count conservation is exact.
- The primary model uses a subject random intercept and must pass fixed-effect
  rank, age support, residual-degree-of-freedom, and repeated-subject checks.
- The one-sample-per-subject sensitivity uses the pre-specified deterministic
  rule and a fixed-effect model without a subject random intercept.
- Primary and sensitivity results remain separate. Concordance is a robustness
  diagnostic, not automatic biological validation.

## scVI Checkpoint Reconstruction

- Architecture and training arguments must be explicit in the resolved profile;
  package defaults are not an acceptable durable execution contract.
- The reconstructed checkpoint must preserve the exact input cell and gene
  identifiers, order, and shape.
- `X_scVI` must have the configured latent dimension and only finite values.
- Stored provenance must match the resolved seed, layer, covariates, model
  arguments, training arguments, accelerator, devices, and config fingerprint.
- A non-empty persisted scVI model and a passing machine-readable checkpoint
  report are required before clustering.

## Single-Cell Constraints

- Large sparse matrices must not be densified without a documented and bounded
  memory justification.
- Differential or association analyses must account for biological
  replicates; cell counts must not inflate the apparent sample size.
- Outputs generated from demo fixtures are not biological evidence.
- Stochastic stages must use the configured seed and record it with outputs.

## Clustered Checkpoint Reconstruction

- Neighbours, UMAP, and Leiden parameters, backend, worker limit, resolution,
  iteration count, and seed must be explicit and preserved in provenance.
- The clustered checkpoint must preserve the exact ordered cell and gene
  identities of the accepted scVI checkpoint.
- UMAP coordinates and sparse graph values must be finite, graphs must match
  the cell-by-cell shape, and every cell must receive a Leiden label.
- A passing machine-readable cluster report is required before CellTypist
  annotation. Similarity to a historical partition is diagnostic evidence, not
  an automatic pass/fail threshold or biological validation.

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
