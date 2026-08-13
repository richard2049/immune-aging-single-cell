# Corrected Results Review

> [!IMPORTANT]
> This record describes a historical review. Its public authorization is
> suspended while the annotated checkpoint and structural-zero composition
> analysis are rebuilt under the current verification contract. Results require
> renewed validation and human review before promotion.

## Review Status

**Corrected public wording and assets approved.**

This document records the pre-specified evidence screen, human dispositions, and final approval. Public interpretation is limited to the exact retained signals, wording, and assets listed in the promotion record.

- Corrected output root: `results\gse164378_full_replicate_corrected`
- Required grouping field: `biological_replicate_id`
- Evidence rows: 109
- Reviewer: `Repository owner`
- Review date: `2026-07-30`
- Public scope: `conservative_core`
- Candidate dispositions complete: **Yes**
- Public changes authorized: **Yes**

## Screening Criteria

Association rows are candidates for human review only when all configured
conditions pass:

- corrected primary FDR < 0.05;
- corrected effect interval excludes zero;
- at least 20 biological replicates;
- testable in all 5 sensitivity scenarios;
- effect direction agrees in all scenarios;
- primary FDR remains below the threshold in all scenarios;
- the main and sensitivity tables use `biological_replicate_id`.

These criteria are deliberately conservative. Failure means that a row is not
prioritized or remains sensitivity-limited; it does not prove the null.

## Automated Screen

- `not_prioritized`: 86
- `candidate_for_human_review`: 12
- `exploratory_sensitivity_limited`: 10
- `candidate_internal_performance_only`: 1

### Association Candidates

| analysis | cell_type | signature | n_replicates | effect_per_10y | effect_ci_low | effect_ci_high | primary_fdr | sensitivity_fdr_supported_count | sensitivity_expected_scenarios | human_disposition | presentation_tier |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| composition | Tcm/Naive cytotoxic T cells |  | 316 | -0.01768 | -0.01961 | -0.01578 | 5.152e-53 | 5 | 5 | retain | core |
| composition | MAIT cells |  | 315 | -0.003036 | -0.003777 | -0.002325 | 2.706e-15 | 5 | 5 | retain | secondary |
| composition | CD16+ NK cells |  | 316 | 0.007552 | 0.004178 | 0.01081 | 6.432e-05 | 5 | 5 | retain | core |
| composition | Regulatory T cells |  | 316 | 0.001245 | 0.0007713 | 0.001728 | 6.432e-05 | 5 | 5 | exploratory | exploratory |
| composition | Tem/Effector helper T cells |  | 316 | 0.003818 | 0.002135 | 0.005555 | 0.0007953 | 5 | 5 | retain | secondary |
| composition | Non-classical monocytes |  | 316 | 0.00273 | 0.001439 | 0.004033 | 0.005374 | 5 | 5 | retain | secondary |
| composition | CD8a/a |  | 86 | -7.984e-05 | -0.0001178 | -4.286e-05 | 0.009501 | 5 | 5 | exclude | excluded |
| signature | CD16+ NK cells | proteostasis_upr | 286 | -0.005851 | -0.008595 | -0.003107 | 0.001856 | 5 | 5 | exploratory | exploratory |
| signature | CD16+ NK cells | inflammatory_nfkb | 286 | -0.01893 | -0.02901 | -0.008656 | 0.002717 | 5 | 5 | exploratory | exploratory |
| signature | Classical monocytes | inflammatory_nfkb | 302 | -0.04252 | -0.06706 | -0.01817 | 0.008607 | 5 | 5 | exploratory | exploratory |
| signature | Classical monocytes | sasp_proxy | 302 | -0.01856 | -0.02996 | -0.006868 | 0.009101 | 5 | 5 | exploratory | exploratory |
| signature | Tem/Effector helper T cells | inflammatory_nfkb | 268 | -0.01868 | -0.03104 | -0.006296 | 0.0267 | 5 | 5 | exploratory | exploratory |

`retain` means eligible for cautious drafting, not approved public wording.
`exploratory` results must not be used as headline findings. The excluded
`CD8a/a` result had only 144 cells in total and a median of one cell per
contributing replicate despite passing the automated association screen.

## Internal Age-Prediction Evidence

The selected `nested_selected` model used grouped
cross-validation over 316 biological replicates
and 5 outer folds. Its donor-level MAE was
12.39 years (95% interval
11.63-13.15), versus a
fold-specific mean-age baseline MAE of
15.13 years; the internal improvement was
2.74 years. R-squared was
0.288 (95% interval
0.224-0.344).

The regressor was cross-validated by biological replicate, but the upstream
`X_scVI` representation was learned once from the full cohort. This therefore
supports only an internal, transductive predictive signal. It is not an
end-to-end unseen-donor evaluation or external validation and must not be
described as a clinical biomarker or validated aging clock. Its approved presentation tier is
`technical_only`.

## Interpretation Boundaries

- Composition effects are cross-sectional fraction changes per 10 years, not
  causal effects of aging.
- Signature scores are configured gene-set proxies, not direct measurements of
  pathway activation or suppression.
- FDR, confidence intervals, support, and sensitivity are screening evidence;
  biological plausibility and confounding still require human assessment.
- Provisional-versus-corrected comparison status is context only and is not a
  selection criterion.
- Age, batch, sex, annotation uncertainty, sparse populations, and unequal
  cells per replicate remain possible limitations.

## Human Review Checklist

- [x] Reconcile candidate numbers with the source CSVs.
- [x] Review replicate, age-range, and cell-count support per candidate.
- [x] Assign each prioritized row `retain`, `exploratory`, `exclude`, or
      `defer`.
- [x] Approve exact claim wording and captions without causal or clinical language.
- [x] Approve each proposed asset replacement.
- [x] Record reviewer, date, and rationale in the reviewed disposition record.

## Promotion State

The machine-readable promotion manifest records the exact authorization state.
The approved wording and captions are preserved in
`docs/reviews/corrected_public_promotion_record.md`; generated analytical
outputs remain separate from curated public assets.
