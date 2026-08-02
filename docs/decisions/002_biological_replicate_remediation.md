# Decision Card: Biological Replicate Remediation And Result Regeneration

## Metadata

- Stage: Biological replicate remediation and donor-level result regeneration
- Priority: High
- Status: Implemented - corrected interpretation review completed
- Risk classification: C - Scientifically important

## Biological Question

Which immune-aging composition, signature, and age-prediction findings remain
supported when every analysis uses the correct independent biological
replicate?

## Why This Stage Matters

The feasibility audit established that `donor_id` contains 166 labels reused
across source pools, while `Tube_id` identifies 317 conflict-free biological
replicates. Existing donor-level analyses therefore merge distinct individuals.
Correcting this unit is a prerequisite for valid age associations,
donor-grouped cross-validation, sensitivity analysis, and pseudobulk
differential expression.

This stage restores validity; it should precede any new biological analysis.

## Required Inputs And Metadata

- `results/gse164378_full/06_annotated.h5ad` as the validated expression,
  embedding, clustering, and annotation checkpoint.
- `data/raw/raw_counts_h5ad/all_pbmcs/all_pbmcs_metadata.csv`.
- Cell-level join key corresponding to AnnData `obs_names`.
- Source `Tube_id`, `donor_id`, `Age`, `Sex`, `Batch`, and `File_name`.
- `config/config.real.full.yml`.
- Feasibility outputs under
  `results/gse164378_full/pseudobulk_de_feasibility/`.
- Current provisional donor-level tables and figures for before/after
  comparison.

## Expected Outputs

- A validated canonical `biological_replicate_id` derived from `Tube_id`, with
  source and conflict fields.
- A machine-readable cell-to-replicate mapping and replicate-level metadata
  audit.
- Versioned, corrected composition, signature, age-prediction, sensitivity,
  and supplementary tables and figures.
- A comparison report describing which provisional results remain stable,
  change direction, lose support, or become newly supported.
- Updated provenance and public documentation after scientific review.
- No overwrite of provisional outputs until the corrected run is validated.

## Key Scientific And Statistical Decisions

- Treat `Tube_id` as the independent biological replicate and retain
  `donor_id` only as a non-unique source label.
- Treat `sample_id` / `File_name` as multiplexed technical-library metadata,
  not as independent donors.
- Aggregate or group every age analysis by
  `biological_replicate_id`.
- Group all age-prediction cross-validation folds by
  `biological_replicate_id` to prevent leakage.
- Retain sex and batch adjustment only when the corrected design is estimable.
- Reuse the validated scVI embedding, clustering, and CellTypist annotation if
  verification confirms that none of those stages used `donor_id` as a model
  input. Avoid unnecessary retraining that would confound metadata correction
  with a new representation.
- Write corrected outputs to a distinct versioned location and compare them
  with provisional outputs before replacing curated assets or README results.
- Do not interpret changed biological findings during implementation; that
  review is a separate D-class task.

## Risks And Possible Artifacts

- An incomplete or order-dependent cell-level join could assign the wrong
  replicate.
- Technical libraries may contribute unequal numbers of cells, creating
  library-size and batch imbalance.
- Correcting the unit from 166 reused labels to 317 replicates can materially
  change effect estimates, confidence intervals, FDR values, and prediction
  error.
- Rare cell types may gain or lose eligibility after regrouping.
- Failure to update every donor-level module could leave mixed-key outputs.
- Reusing stale curated figures or README text could present invalidated
  results.
- Retraining scVI unnecessarily could make before/after differences
  uninterpretable because both metadata and the learned representation changed.

## Minimal Validation

- Map exactly 1,000,000 of 1,000,000 cells to a non-missing
  `biological_replicate_id`.
- Confirm 317 unique biological replicates and zero within-replicate conflicts
  for donor label, age, sex, and batch.
- Confirm cell identifiers are unique and the join is key-based, never
  positional.
- Confirm `adata.X`, `X_scVI`, cell-type labels, and cell order are unchanged
  when only metadata is remediated.
- Confirm all donor-level configs and output provenance name
  `biological_replicate_id`.
- Confirm no biological replicate appears in both training and test data in
  any age-prediction fold.
- Run only the affected downstream targets and verify their expected tables,
  figures, row counts, and exclusion reports.
- Produce a structured before/after comparison and keep all corrected
  interpretations provisional until human review.

## Classification

**C - Scientifically important.** The stage changes the experimental unit and
therefore the validity of every donor-level estimate. The engineering work is
delegable with strict tests, but the identity definition, before/after
assessment, and approval of changed conclusions require human review.

Any subsequent decision to publish or biologically interpret changed findings
is **D - High-risk interpretation or claim**.

## Decision Card

- Required: Yes
- Rationale: The stage changes the independent biological replicate, affects
  existing scientific outputs, and gates pseudobulk DE.
- Approval status: Approved and implemented

## Decision

Implement the remediation and controlled regeneration using `Tube_id` as
`biological_replicate_id`, while keeping all resulting biological
interpretations provisional until separate human review.

## Rationale And Tradeoffs

Correcting the replicate key before adding new analyses prevents invalid
grouping and cross-validation leakage from propagating into pseudobulk DE.
Reusing the validated representation avoids an unnecessary full scVI retrain,
but only if code/config inspection confirms that upstream representation stages
did not depend on the reused label. Versioned outputs require additional disk
space but preserve an auditable before/after comparison.

## Interview Question

- Question: How would you correct a reused donor identifier in an established
  single-cell workflow without unnecessarily retraining upstream models?
- Skill tested: Experimental-unit definition, metadata provenance, checkpoint
  reuse, leakage prevention, and reproducible result migration.
- Expected answer should mention: Validate a stable source identifier such as
  `Tube_id`; join by cell ID rather than row order; expose a canonical
  biological-replicate field; update every donor-level grouping and CV split;
  verify upstream models did not use the old identifier; version corrected
  outputs; compare old and new results before changing claims.
- Status: Ready

## Follow-up Biological Question

- Question: Which provisional immune-aging findings remain stable after
  regrouping cells by the correct biological replicate?
- Biological theme: Robustness of immune-aging associations to replicate
  definition.
- Why it matters: Stability after correction determines which composition,
  signature, and predictive signals can support further biological analysis.
- Possible next analysis: Compare corrected and provisional effect sizes,
  confidence intervals, FDR values, eligible cell types, and donor-grouped
  prediction errors.
- Status: Open

## No-AI Survival Exercise

- Estimated time: 20-30 minutes
- Skill tested: Detecting identifier reuse and defining an auditable biological
  replicate.
- Allowed resources: This decision card,
  `results/gse164378_full/pseudobulk_de_feasibility/replicate_metadata_audit.csv`,
  `config/config.real.full.yml`, and a spreadsheet or short local script; no
  generative AI.
- Expected output: A one-page replicate-key validation stating why
  `donor_id`, `sample_id`, and `Tube_id` have different roles, plus the checks
  required before accepting `biological_replicate_id`.
- Pass criteria: Identifies `Tube_id` as the independent replicate, recognizes
  multiplexed technical libraries, requires key-based joins and conflict
  checks, and explains why cross-validation must group by the corrected key.
- Follow-up biological question: Which previous age associations are most
  sensitive to the correction in replicate definition?
- Status: Not started

## Implementation Evidence

- Corrected profile:
  `config/config.real.full.replicate_corrected.yml`.
- Targeted workflow:
  `workflows/Snakefile.replicate_corrected`.
- Corrected output root:
  `results/gse164378_full_replicate_corrected/`.
- Mapping audit: 1,000,000 of 1,000,000 cells mapped to 317 biological
  replicates, with zero within-replicate conflicts for donor label, age, sex,
  and batch.
- Upstream checkpoint reuse: the existing annotated H5AD was opened read-only;
  no expression, embedding, clustering, or annotation checkpoint was rewritten.
- Corrected analyses: composition, signature associations, grouped age
  prediction, five sensitivity scenarios, supplementary figures, and a
  descriptive before/after comparison.
- CV validation: no biological replicate spans more than one fold within a
  model/evaluation level.
- Workflow validation: the targeted Snakemake dry-run reports all requested
  corrected outputs present and up to date.
- Interpretation boundary: comparison outputs remain descriptive. The separate
  D-class review is complete in decision card 003, and only its explicitly
  retained corrected results were promoted.
