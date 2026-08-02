# Decision Card: Scientific Reproducibility Remediation

## Metadata

- Stage: Metadata and model reproducibility remediation
- Priority: High
- Status: Implemented; local environment and Git recovery remain external actions
- Risk classification: C - Scientifically important
- Verification level: L3 - Claim-critical

## Biological Question

Can the immune-aging workflow preserve the intended biological replicate,
model design, and evaluation scope consistently across machines and reruns?

## Why This Stage Matters

The repository audit found failure modes that could silently alter sample
identity or model specification: ambiguous metadata keys were deduplicated,
one metadata helper aggregated a reused donor label, configured scVI
covariates could be omitted, and optional prediction dependencies could change
the candidate model set. Random seeds and external annotation-model provenance
were also incomplete. These are reproducibility defects even when existing
result files remain technically readable.

This remediation strengthens the workflow contract. It does not reinterpret
existing biological results or approve new public claims.

## Required Inputs And Metadata

- Global and repository `AGENTS.md` instructions.
- `IMMUNE-AGING_NOTES.md` and `docs/verification_contract.md`.
- Metadata join keys and authoritative source fields.
- Configured biological-replicate, age, batch, and model covariates.
- Existing demo and corrected full-data profiles.
- Existing test fixtures and validated corrected-result manifests.

## Expected Outputs

- Metadata joins that reject conflicting duplicate keys and validate
  many-to-one cardinality.
- A metadata-table helper that uses an explicit configurable biological
  replicate key and writes only to derived-data locations.
- Explicit failure when required scVI covariates or requested prediction-model
  dependencies are unavailable.
- Configured random seeds for scVI, SOLO, neighbors, UMAP, and Leiden steps.
- Recorded CellTypist model identity and package/runtime provenance in AnnData.
- Environment and documentation contracts consistent with configured models.
- Focused tests covering failure boundaries and reproducibility settings.

## Key Scientific And Statistical Decisions

- Exact duplicate metadata rows may collapse to one row; duplicate join keys
  with conflicting payloads are invalid and must stop the workflow.
- Every metadata merge must preserve cell count and cell order and validate a
  many-to-one relationship from cells to metadata.
- Derived replicate metadata must use an explicitly named key. For GSE164378,
  that key is `Tube_id` exposed as `biological_replicate_id`; reused
  `donor_id` labels are not independent biological units.
- Covariates listed in the scVI configuration are required by default. Their
  silent removal would constitute a different model.
- A requested prediction model is part of the analysis contract. Missing
  libraries must cause an actionable error rather than silently changing the
  comparison set.
- A single configured seed is propagated to stochastic workflow stages and
  recorded with their outputs.
- Current age-prediction estimates use donor-grouped cross-validation of the
  regressor but a latent representation trained once on the complete cohort.
  They are therefore internal, transductive estimates rather than a fully
  end-to-end inductive evaluation for unseen donors.

## Risks And Possible Artifacts

- Stricter validation can expose previously hidden metadata conflicts and stop
  runs that used to continue.
- Re-running stochastic upstream stages after adding explicit seeds can change
  embeddings, clusters, annotations, and downstream results.
- Rebuilding the environment can change numerical output unless dependency
  versions are controlled and the regenerated outputs are reviewed.
- CellTypist model files obtained by model name can differ if an external cache
  changes; model provenance must therefore be recorded at execution time.
- Full fold-specific retraining of scVI would be substantially more expensive
  and is outside this remediation. Existing prediction claims must retain the
  internal/transductive limitation.

## Minimal Validation

- Positive join fixture with unique keys preserves row count and order.
- Exact duplicate metadata rows are accepted; conflicting duplicate keys fail.
- Missing configured scVI covariates fail before model training.
- Missing requested optional prediction dependencies fail with installation
  guidance.
- Seed propagation is verified without running a heavy model.
- Replicate aggregation tests demonstrate that reused donor labels remain
  separate when the configured biological-replicate key differs.
- Config files parse and the Python source compiles.
- Existing independent biological-replicate and pseudobulk acceptance tests
  continue to pass without weakened assertions.

## Classification

**C - Scientifically important / L3 - claim-critical.** The changes affect
metadata cardinality, experimental-unit handling, model specification, and the
scope of predictive claims. No biological result will be promoted as part of
this implementation.

## Decision Card

- Required: Yes
- Rationale: The remediation changes failure behavior at scientific validity
  boundaries and establishes durable reproducibility contracts.
- Approval status: Approved by the user's request to correct the audited
  deficiencies and verify each correction.

## Decision

Implement strict, explicit, config-driven validation while preserving existing
validated outputs as immutable checkpoints. Document the transductive
prediction boundary instead of presenting the existing analysis as an
end-to-end unseen-donor evaluation.

## Rationale And Tradeoffs

Failing early is preferable to silently changing biological units or model
design. The stricter behavior may require users to correct metadata or install
declared dependencies, but it makes runs comparable and auditable. A complete
fold-specific scVI retraining experiment remains a separate, computationally
expensive C-stage analysis rather than being mixed into this repair.

## Interview Question

- Question: Which workflow checks prevent metadata and environment differences
  from silently changing a single-cell analysis across machines?
- Skill tested: Metadata cardinality, experimental-unit validation, dependency
  contracts, stochastic reproducibility, and interpretation boundaries.
- Expected answer should mention: many-to-one joins; rejection of conflicting
  duplicate keys; explicit biological-replicate identifiers; required
  covariates and dependencies; controlled seeds; external-model provenance;
  and the distinction between grouped regressor CV and end-to-end inductive
  representation learning.
- Status: Ready

## Follow-up Biological Question

- Question: Do the retained immune-aging findings remain stable when the full
  representation-learning and prediction procedure is evaluated on unseen
  donors?
- Biological theme: External validity of latent immune-aging representations.
- Why it matters: A donor-grouped regressor cannot by itself establish that a
  representation trained on the full cohort generalizes to new donors.
- Possible next analysis: Compare the current internal estimate with a bounded
  fold-specific scVI retraining experiment or an independent-cohort analysis.
- Status: Open

## No-AI Survival Exercise

- Estimated time: 20-30 minutes
- Skill tested: Auditing a scientific workflow for silent design changes.
- Allowed resources: This decision card, `environment.yml`, one real-data
  config, and the metadata/scVI/age-prediction modules; no generative AI.
- Expected output: A short table listing one metadata, one dependency, one
  stochastic, and one validation-scope failure mode, with the required safe
  behavior for each.
- Pass criteria: Identifies conflicting join keys, reused replicate labels,
  missing configured dependencies or covariates, seed propagation, and the
  transductive limitation of a full-cohort latent representation.
- Follow-up biological question: Which reported associations or prediction
  metrics are most sensitive to end-to-end unseen-donor evaluation?
- Status: Not started

## Implementation Evidence

- Strict metadata join and replicate-helper boundaries are covered by focused
  positive and negative fixtures in `tests/test_reproducibility_contracts.py`.
- Configured covariates and prediction-model dependencies now fail explicitly
  when unavailable.
- `run.seed` is explicit in each primary profile and is propagated to scVI,
  SOLO, neighbors, UMAP, and Leiden; model-stage provenance is stored in
  AnnData.
- CellTypist provenance includes the resolved model SHA-256 and enforces an
  optional configured checksum.
- The tracked environment declares version-constrained scVI, CellTypist, and
  XGBoost dependencies and passed a `mamba` dry-run solve.
- The full unittest suite ran 33 tests: 32 passed, and the opt-in Docker edgeR
  integration test was skipped without changing its contract.
- Detailed verification and manual recovery instructions are recorded in
  `docs/validation/scientific_reproducibility_remediation.md`.
- Git history could not be recovered because the local `.git` directory is
  empty and no authoritative remote URL or backup is recorded in the working
  tree.
