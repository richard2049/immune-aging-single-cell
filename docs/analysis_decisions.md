# Analysis Decisions

This summary records the current scientific scope and should be read together
with the [verification contract](verification_contract.md).

## Study And Units

The analysed cohort is the Terekhova et al. Blood Age Atlas (`syn49637038`),
not GEO accession `GSE164378`. Source `Donor_id` identifies the subject,
`Tube_id` identifies a longitudinal sample unit, and `File_name` identifies a
technical library. Cells and repeated sample units from one subject are not
independent people.

Historical analyses promoted `Tube_id` to an independent biological replicate.
Those numerical results and presentation decisions are superseded and remain
non-promotable until the current workflow is regenerated and reviewed.

## Composition And Signatures

Composition and signature scores are aggregated by sample unit and cell type
while retaining the subject identifier. Primary associations use
subject-clustered GEE and require support in distinct subjects. Composition
restores structural zero counts before calculating fractions.

Each table also records a deterministic one-sample-per-subject sensitivity,
selecting earliest observed age and then lexical sample ID. This sensitivity
does not replace the repeated-measures model, and agreement is not automatic
biological validation.

## Age Prediction

Outer and inner cross-validation folds are grouped by `subject_id`; all visits
and cell-type rows from one person remain in one fold. Training weights balance
subjects, and bootstrap intervals resample subjects. The scVI representation is
still learned from the full cohort, so performance is internal and transductive,
not fully inductive or externally validated.

## Pseudobulk

Raw counts are summed by `sample_unit_id x cell_type` with exact count
conservation. The primary model uses `voomWithDreamWeights`/`dream` with fixed
effects for sex, batch, and age per decade plus `(1 | subject_id)`. A
non-identifiable fixed design, insufficient age support, insufficient residual
degrees of freedom, or absence of repeated subjects stops that cell-type model.

The pre-specified sensitivity selects one sample per subject and fits the same
fixed effects without a random intercept. Primary and sensitivity outputs are
separate and require technical validation before gene-level review.

## Checkpoints And Promotion

New checkpoints and outputs use `blood_age_atlas` namespaces. Historical
`results/gse164378*` directories are immutable provenance. The longitudinal
workflow no longer invokes historical comparison, disposition, or public-asset
generators.

Technical acceptance does not authorize a biological claim. Cell-population,
gene, pathway, predictive, causal, or clinical interpretation requires a later
D-stage human review.

## Interpretation Boundary

The primary estimand is a population-level association with chronological age
while accounting for repeated observations. It is not a within-person ageing
effect. The data and models do not establish causality, biological age,
diagnostic performance, or clinical utility.
