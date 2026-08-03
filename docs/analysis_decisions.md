# Analysis Decisions

This document summarizes the main analysis choices that determine the
scientific scope of the public workflow. It records approved methods and
interpretation boundaries without replacing the linked design and validation
evidence.

## Biological replicate definition

The source `donor_id` values were reused across pools and therefore could not
serve as independent donor-level units. The corrected workflow uses `Tube_id`
as `biological_replicate_id`. This definition produced 317 biological
replicates across 1,000,000 cells with no detected metadata conflicts, and it
is used for donor-level summaries and cross-validation grouping. The upstream
annotated AnnData checkpoint is reused read-only because its learned cell-level
representation did not depend on the invalid donor label.

Evidence: [corrected-results review](corrected_results_review.md) and
[scientific verification contract](verification_contract.md).

## Donor-aware pseudobulk inference

Raw counts are aggregated by biological replicate and cell type. The primary
edgeR quasi-likelihood model tests continuous age per decade while adjusting
for sex and batch. A simpler age-only model is allowed only when the adjusted
design is not estimable and must be reported explicitly. Profiles require at
least 50 cells, 12 biological replicates, 12 years of age coverage, and five
residual degrees of freedom. TMM normalization and design-aware gene filtering
are applied before inference.

Multiple testing is reported both within cell type and globally across
gene-cell-type tests. Twelve predefined populations form the primary scope;
NK cells remain exploratory because of annotation granularity. In the accepted
full run, 949,518 eligible cells yielded 2,944 profiles, 317 replicates, 36,601
genes, and 107,948 tests. All 13 fitted adjusted designs were full rank and no
fallback model was used.

Evidence: [analysis design](pseudobulk_de_design.md), [workflow details](pseudobulk_de.md),
and [one-million-cell run acceptance](validation/pseudobulk_de_1m_run_acceptance.md).

## Corrected-result promotion

Automated evidence screening prioritizes results for review but does not
decide which findings are suitable for presentation. Human review of 109
evidence rows assigned the 13 priority associations to six retained, six
exploratory, and one excluded disposition. Three figures were approved for
public presentation; the signature heatmap was not selected. Signature-score
associations remain exploratory and are not direct measurements of pathway
activity.

Evidence: [corrected-results review](corrected_results_review.md),
[result dispositions](reviews/corrected_results_dispositions.yml), and the
[public-promotion record](reviews/corrected_public_promotion_record.md).

## Predictive and reproducibility scope

Metadata joins and replicate definitions are validated explicitly, configured
covariates and model dependencies fail with actionable errors, and stochastic
steps receive configured seeds. External annotation resources are recorded
with version and checksum information where applicable. The age-prediction
analysis uses replicate-grouped nested cross-validation, but its scVI
representation was learned from the full cohort. It is therefore an internal,
transductive evaluation rather than an end-to-end assessment in unseen donors,
and it has no established clinical validity.

Evidence: [scientific verification contract](verification_contract.md) and
[corrected-results review](corrected_results_review.md).

## Robustness before gene-level interpretation

Predefined pseudobulk sensitivity analyses raise the minimum cell count to 100,
omit sex, omit batch, and leave one batch out when the design remains
estimable. These analyses diagnose dependence on modeling choices and do not
replace the primary model. The candidate universe is defined by global FDR
below 0.05; candidates are preserved without an additional discovery filter,
and no automated label declares a result biologically validated.

The completed audit fitted 212 models, recorded three non-estimable designs and
no execution failures, and reduced 7,970 statistical candidates to a bounded
249-row review queue. All 13 diagnostic plot sets were visually reviewed.
Gene- and pathway-level interpretation remains pending and requires separate
human scientific review.

Evidence: [robustness audit](pseudobulk_robustness_audit.md) and
[diagnostic plot review](reviews/pseudobulk_diagnostic_plot_review.yml).

## Interpretation boundary

The study is cross-sectional. Associations with chronological age do not
estimate within-person change and do not establish causal effects, biological
age, diagnostic performance, or clinical utility. Technical validation makes
the outputs reviewable; it does not by itself authorize biological claims.
