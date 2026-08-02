# Corrected Public-Promotion Record

## Status

**Approved source text and captions for corrected public promotion.**

- Public scope: conservative core
- Signature heatmap: not selected for curated public assets
- Public claims authorized: Yes
- Public asset replacement authorized: Yes

## Proposed README Text

### Replicate-corrected GSE164378 analysis

Review of the source metadata showed that the original `donor_id` field was
reused across pools. Donor-level analyses were therefore repeated using
`Tube_id` as `biological_replicate_id`. The source metadata contained 317
biological replicates, of which 316 met the requirements for the main
composition and prediction analyses.

After adjustment for sex and batch, a 10-year difference in age was associated
with a 1.77 percentage-point lower Tcm/Naive cytotoxic T-cell fraction (95% CI
-1.96 to -1.58; FDR 5.15e-53), a 0.30 percentage-point lower MAIT-cell
fraction (95% CI -0.38 to -0.23; FDR 2.71e-15), and a 0.76 percentage-point
higher CD16+ NK-cell fraction (95% CI 0.42 to 1.08; FDR 6.43e-5). These are
cross-sectional associations between participants, not estimates of
within-person change or causal effects of aging. Tcm/Naive cytotoxic T cells
form a combined annotation category rather than a single resolved subtype.

The scVI-derived features also retained an age-related predictive signal. In
nested cross-validation grouped by biological replicate, the selected model
had a donor-level mean absolute error of 12.39 years (95% CI 11.63 to 13.15),
compared with 15.13 years for a fold-specific mean-age baseline. This is modest
internal predictive performance and has not been validated as an aging clock
or clinical biomarker.

Several predefined gene-set scores were associated with age, but these results
remain exploratory. Signature scores are proxies for the configured gene sets
and do not directly measure pathway activation or suppression.

## Proposed Figure Captions

### Selected Composition Trends

Cell-type fractions per biological replicate plotted against age for the three
associations selected for the main figure. Lines are unadjusted linear
summaries included for visualization; Spearman rho and FDR are from analyses
adjusted for sex and batch. Tcm/Naive cytotoxic T cells form a combined
annotation category. The associations are cross-sectional and do not establish
causality.

Approved source: `results\gse164378_full_replicate_corrected\review\figures\composition_core_trends_draft.png`

### Internal Age Prediction

Observed donor age and cross-validated predicted age for the selected model.
Outer folds were grouped by biological replicate. The identity line represents
perfect prediction, while the fitted line shows the relationship observed in
the held-out predictions. Reported metrics describe internal cross-validation
and do not establish external validity or clinical utility.

Approved source: `results\gse164378_full_replicate_corrected\review\figures\age_prediction_internal_cv_draft.png`

### Retained Composition Effects

Sex- and batch-adjusted differences in cell-type fraction per 10-year
difference in age for all retained composition associations. Points show
effect estimates in percentage points, and bars show 95% bootstrap confidence
intervals.
Filled markers identify core presentation results; open markers identify
secondary retained results.

Approved source: `results\gse164378_full_replicate_corrected\review\figures\composition_effect_forest_draft.png`

## Final Approval

- [x] README wording approved with professional tone revision.
- [x] Three captions approved.
- [x] Three draft images approved after visual inspection.
- [x] Provisional curated assets approved for removal.
- [x] Approved corrected assets authorized for `docs/assets/`.
