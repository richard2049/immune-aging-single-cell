# Annotation Mapping Human-Review Guide

## Purpose

This guide defines the outcome-blind human review required before the Blood
Age Atlas annotation mapping can authorize longitudinal inference. The review
approves how each CellTypist label is handled analytically; it does not claim
that automated annotation establishes biological identity.

Do not inspect age-association estimates, effect directions, p-values, FDR,
prediction metrics, or differential-expression results during this review.

## Files And Ownership

Generated evidence under `results/` is read-only:

- `annotation_mapping_review_table.csv`: one-row-per-raw-label review index;
- `annotation_raw_label_marker_summary.csv`: raw-count marker evidence by raw label;
- `annotation_cluster_crosstab.csv`: raw and proposed labels by Leiden cluster;
- `annotation_design_eligibility.csv`: method-specific support and estimability;
- `annotation_design_qualification.json`: fingerprints and qualification status;
- `annotation_design_review.md`: compact qualification summary.

Human conclusions are written only to:

- `docs/reviews/annotation_mapping_review.xlsx` as the compact working sheet;
- `docs/reviews/annotation_mapping_review.yml` during review;
- `config/annotation_mappings/blood_age_atlas_celltypist_v1.csv` after the
  review record is complete and internally consistent;
- decision card 012 when final approval or a revised mapping is recorded.

Never manually edit generated evidence tables or the qualification report.
Start in the workbook rather than the YAML. Its first 17 rows are the primary
proposals, evidence columns are fixed, and yellow columns are the human-review
fields. The YAML remains the canonical final record and should be updated from
the completed workbook before approval.

## Assessment Vocabulary

Use these values in each `human_assessment` field:

- `coherent` or `adequate`: evidence supports the proposed treatment;
- `inconclusive`: available evidence cannot resolve the question;
- `incoherent` or `inadequate`: evidence contradicts or does not support it;
- `not_assessable`: the configured evidence does not address the question;
- `not_applicable`: the criterion does not apply;
- `not_started`: review has not occurred.

Use these `final_decision.status` values:

- `accept_as_proposed`: retain the proposed analysis label and disposition;
- `modify`: change the analysis label, disposition, or both;
- `reject_pending_resolution`: reject the proposal but do not yet have a
  defensible replacement;
- `pending_review`: review is incomplete.

`primary` means that cells enter the stated broad analysis population.
`exploratory` means that the raw label is retained but excluded from primary
inference. `unresolved` means that no defensible analysis label is assigned.

### Method-Scoped Primary Use

`primary` is not blanket authorization for every downstream endpoint. A
population may be primary only for methods whose executable eligibility
contract it passes. Reports, tables, and figure captions must identify any
method restriction that materially affects interpretation.

Under the current review, the DC2-derived `Dendritic cells` population is
eligible for primary composition analysis only. It must not be used for
primary signature, age-prediction, pseudobulk, gene-level, pathway-level, or
DC2-versus-DC3 conclusions. A composition result would support only a cautious
statement about the relative abundance of the approved broad dendritic-cell
mapping, subject to the longitudinal model, denominator decision, annotation
uncertainty, and human result review. It would not establish dendritic-cell
function, transcriptional ageing, causality, or a specific cDC2 mechanism.

The `Regulatory T cells` population denotes a transcriptomic Treg annotation,
not experimentally demonstrated suppressive function. It is eligible for
composition and pseudobulk under the current contracts. Signature estimates
may be reported only with their limited subject support and corresponding
uncertainty; age prediction is not eligible. Do not infer thymic or peripheral
origin, functional suppression, or a finer Treg state without independent
marker or functional evidence.

The `pDC` population is eligible for primary composition analysis only. Its
strong annotation evidence does not compensate for insufficient per-profile
support for signatures, age prediction, or pseudobulk. A reported composition
association may describe relative pDC abundance under the approved model, but
must not be presented as evidence of interferon activity, molecular pathway
change, differential expression, predictive value, or pDC function.

If later evidence supports a wider use, revise and reapprove the mapping and
method-eligibility record before regenerating the affected outputs.

### Required NK Mapping Sensitivity

The low-confidence raw `NK cells` label and the less certain CD16-negative
assignment may contribute to the aggregated primary `NK cells` population,
but neither may be analysed as a standalone primary population. Every
reported broad-NK association must be compared across three prespecified
mappings: all approved NK assignments; exclusion of the raw generic `NK cells`
label; and a restrictive bound retaining only the strongly supported
CD16-positive assignment. The last variant is a robustness bound, not a claim
that CD16-positive cells represent the complete NK compartment.

Report changes in effect direction, magnitude, uncertainty, and method
eligibility. If the result changes materially, prioritize or explicitly
qualify the restricted estimate before making a biological claim. Do not
define materiality after inspecting whichever version gives the preferred
result.

## Review One Raw Label

For the current raw label, open its row in
`annotation_mapping_review_table.csv` and the corresponding block in the YAML.

### 1. Marker Coherence

Filter `annotation_raw_label_marker_summary.csv` by `raw_label`. Inspect
`expected_marker`, `gene_present`, `detection_fraction`, and `mean_raw_count`.

- Confirm expression of multiple expected lineage markers where a panel exists.
- Inspect markers from competing lineages; one highly expressed generic gene is
  not sufficient evidence.
- Treat raw counts as descriptive cell-level evidence, not biological replicates.
- Use `not_assessable` when `marker_panel_label` is empty. Promotion of such a
  label requires an independently reviewed marker panel and regenerated evidence.

Do not apply a universal detection-fraction threshold.

### 2. Cluster Coherence

Use `dominant_leiden`, `dominant_leiden_fraction`, and `n_leiden_clusters` from
the review table. Then inspect all matching raw rows in
`annotation_cluster_crosstab.csv`.

- Concentration in one or a few biologically compatible clusters supports coherence.
- Broad dispersion across unrelated clusters is a warning, not an automatic failure.
- `dominant_leiden_fraction` is the fraction of the label in its dominant cluster.
- `cluster_fraction` is the fraction of a cluster assigned to the label; the two
  quantities answer different questions.

### 3. Sample And Subject Support

Inspect `n_cells`, `n_sample_units_detected`, `n_subjects_detected`, and the
method-specific eligible sample-unit and subject fields.

- Cells do not count as independent biological replicates.
- Broad subject and sample support strengthens generalizability.
- Sparse labels may remain biologically plausible while being unsuitable for
  primary inference.
- Do not override a failed executable method contract by changing the review text.

### 4. Technical Distribution

Inspect `n_batches_detected`, `dominant_batch`, and
`dominant_batch_fraction`, together with method-specific design rank and
residual degrees of freedom.

- Strong concentration in one batch requires explanation or demotion.
- Full global design rank does not guarantee estimability for a rare label.
- Batch concentration is a warning requiring judgement; no universal cutoff is used.

### 5. PBMC Plausibility

Decide whether the proposed population is plausible in peripheral blood and
whether its broad parent is taxonomically coherent. Tissue-resident,
non-haematopoietic, developmental, cycling, and hybrid labels require extra
caution. Plausibility alone cannot override contradictory marker or technical evidence.

### 6. Record The Decision

Complete all five `human_assessment` fields. Then set the final decision,
analysis label, disposition, and a concise evidence-based rationale. Do not
change the executable CSV yet.

## Recommended Review Order

Review primary proposals first:

1. B cells: generic, memory, and naive B labels;
2. CD4 and CD8 T-cell compartments;
3. classical and non-classical monocytes;
4. NK, MAIT, and regulatory T populations;
5. DC2 and pDC.

Then review exploratory labels by proposed parent. For every fine label, ask
separately whether the subtype is independently supportable and whether its
cells can safely enter the broad parent. A subtype can be unsuitable for
standalone inference but still have a defensible broad lineage.

Review unresolved labels last. Approving an unresolved row means approving its
exclusion from labelled primary populations, not validating the classifier label.

## Composition Denominator Decision

The approved `composition_denominator_decision` is recorded in the review YAML.

Recommended option:

- `all_qc_passed_cells_with_other_unresolved`: retain every QC-passed cell in
  the denominator and represent cells without a primary analysis label as
  `Other/unresolved`.

`Other/unresolved` is a derived category used only in the primary composition
denominator. It does not overwrite `cell_type`, remove the raw CellTypist
classification, or discard exploratory labels and subtype-validation records.
Those fields remain available for audit, sensitivity analyses, and subsequent
outcome-blind validation. A later promotion changes the versioned analysis
mapping and requires regeneration of affected downstream outputs; it does not
rewrite the accepted annotation checkpoint.

Alternative:

- `primary_mapped_cells_only_with_explicit_label`: calculate fractions only
  among primary-mapped cells and label the estimand explicitly. This requires
  documenting sample-level exclusion fractions and assessing differential loss.

The implementation now uses the approved first behavior and records
`Other/unresolved` as a denominator-only category. It is excluded from
cell-type trend fitting and must not be interpreted as a biological population.

## Method Eligibility Decision

Review `annotation_design_eligibility.csv` and complete
`method_eligibility_review`. The recommended decision is
`accept_current_method_specific_contracts`. Passing one method does not qualify
the population for every method, and failed rows must remain outside that
method unless a separate scientific decision changes the contract.

## Apply The Completed Review

This procedure was completed on 2026-08-26. For any future mapping revision,
repeat the following steps:

1. Confirm `outcome_blind_confirmation: true`, reviewer, and review date.
2. Apply each accepted final label, disposition, and rationale to the mapping CSV.
3. Set each resolved mapping row to `review_status: approved`; leave unresolved
   review work as `pending_review` and rejected proposals as `rejected`.
4. Implement the approved composition-denominator decision.
5. Regenerate only `annotation_design_qualification`.
6. Confirm mapping and evidence fingerprints, 84 unique labels, reconciled cell
   totals, and absence of any approval sentinel before final review.
7. Record the final mapping hash and reviewer confirmation in the YAML.
8. Run `annotation_design_approval`; only then may downstream targets execute.

Any mapping change after approval invalidates the qualification report,
approval sentinel, and all affected downstream outputs.

## Optional Secondary Subtype Validation

After the 17 primary proposals have been reviewed, consider a separate,
outcome-blind validation of B- and T-cell subtypes that may be supportable from
the existing dataset. This is a scientific extension, not a prerequisite for
regenerating results that use only the approved broad populations.

The minimum useful scope is within-lineage reclustering, prespecified subtype
marker panels, concordance with raw CellTypist labels, stability across
subjects, sample units and batches, and orthogonal mapping to a suitable PBMC
reference. Keep finer labels exploratory until this evidence is reviewed.
Distinctions defined mainly by surface protein or tissue localization may need
multimodal or tissue data; circulating PBMC transcriptomes cannot establish
true tissue residency.

If a future primary figure or claim depends on a fine subtype, complete this
validation before generating that output. Otherwise, first freeze the broad
mapping, regenerate the corrected broad-label results, and treat subtype
validation as the next scientific extension.
