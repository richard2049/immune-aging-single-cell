# Decision Card: Corrected Results Interpretation And Promotion

## Metadata

- Stage: Human review of corrected immune-aging results and public promotion
- Priority: High
- Status: Completed - corrected public promotion approved and applied
- Risk classification: D - High-risk interpretation or claim

## Biological Question

Which corrected composition, signature, and age-prediction findings are
sufficiently supported to present as cautious immune-aging results, and which
should remain exploratory or be excluded?

## Why This Stage Matters

The replicate-remediation stage passed technical validation for 1,000,000
cells, 317 biological replicates, and zero grouped-CV leakage. However, the
descriptive comparison reports direction changes, lost or newly supported
associations, and changes in test eligibility after replacing reused
`donor_id` labels with `biological_replicate_id`.

At the start of this stage, the README and curated assets described the
original provisional analysis. The structured review documented here was
required before corrected figures and numerical results replaced them and
before future pseudobulk findings could be integrated into the same narrative.

## Required Inputs And Metadata

- Corrected validation report:
  `results/gse164378_full_replicate_corrected/validation/replicate_correction_validation.json`.
- Corrected composition, signature, and age-prediction tables under
  `results/gse164378_full_replicate_corrected/tables/`.
- Corrected figures under
  `results/gse164378_full_replicate_corrected/figures/`.
- Five corrected sensitivity scenarios and their manifest under
  `results/gse164378_full_replicate_corrected/sensitivity_age/`.
- Descriptive before/after tables under
  `results/gse164378_full_replicate_corrected/comparison/`.
- Biological-replicate audit, cell counts, age range, sex, batch, and
  cell-type support information.
- Analysis settings from
  `config/config.real.full.replicate_corrected.yml`.
- Current README statements and curated assets for explicit
  retain/replace/remove decisions.

## Expected Outputs

- A machine-readable evidence table with one row per reviewed composition,
  signature, or prediction result.
- For each result: corrected effect estimate, uncertainty interval, FDR,
  replicate support, sensitivity status, comparison status, limitations, and
  review disposition.
- A human-readable review report separating observations, statistical
  evidence, biological interpretation, and unsupported speculation.
- A promotion manifest listing which corrected figures and numerical
  statements may replace provisional README assets.
- Explicit records for results retained as exploratory, excluded for weak
  support, or deferred to external validation.
- Updated README and curated assets only after human approval of the review.

## Key Scientific And Statistical Decisions

- Define the primary inferential statistic and multiple-testing scope for each
  analysis before selecting results for presentation.
- Require the corrected biological-replicate analysis as the sole basis for
  current donor-level claims; use provisional results only as sensitivity to
  the identifier correction.
- Review effect magnitude, confidence intervals, FDR, number of contributing
  replicates, age coverage, and sensitivity consistency together rather than
  selecting by FDR rank alone.
- Predefine how much direction or support instability across sensitivity
  scenarios is acceptable for a result to be called robust.
- Treat cell-type composition and predefined signature scores as
  cross-sectional associations, not causal effects of aging.
- Treat signature scores as proxies for configured gene sets, not direct proof
  of pathway activation or suppression.
- Treat age prediction as internal grouped cross-validation, not an externally
  validated aging clock or clinical biomarker.
- Examine sparse cell populations and age-batch imbalance separately before
  promoting their associations.
- Require direct human approval for every public-facing biological statement
  and figure replacement.

## Risks And Possible Artifacts

- Selecting only the largest or most significant corrected effects can create
  post hoc confirmation bias.
- FDR significance can coexist with a small, unstable, or poorly interpretable
  effect size.
- Unequal cells per replicate and uneven cell-type support can affect
  precision despite donor-level grouping.
- Age, batch, sex, and technical-library structure may leave residual
  confounding in this cross-sectional dataset.
- CellTypist labels and broad cell-state categories may conceal annotation
  uncertainty or subtype heterogeneity.
- Signature results depend on gene-set definition, gene coverage, sampling,
  and normalization.
- Prediction performance may reflect cohort structure and regression to the
  mean; it has no demonstrated external transportability.
- Replacing figures without updating captions, paths, and limitations could
  leave mixed provisional and corrected claims in public documentation.

## Minimal Validation

- Confirm every reviewed row names `biological_replicate_id` as its grouping
  field and traces to the validated corrected output root.
- Reconcile every reported number against its source CSV and corrected figure.
- Verify that promoted associations meet the predeclared support, uncertainty,
  FDR, and sensitivity criteria.
- Inspect age, sex, batch, replicate count, and cell-count support for every
  promoted cell type.
- Confirm that no claim uses causal, clinical, rejuvenation, or biomarker
  language unsupported by this cross-sectional internal analysis.
- Confirm that all promoted age-prediction metrics come from
  biological-replicate-grouped folds and are compared with the configured
  baseline.
- Remove or explicitly archive all provisional README statements and assets
  when corrected replacements are approved.
- Require a final human sign-off on the evidence table, captions, and README
  wording.

## Classification

**D - High-risk interpretation or claim.** This stage does not merely compute
new outputs. It decides which statistical results receive biological meaning
and public visibility. Those decisions cannot be delegated to automated
ranking and require direct human scientific review.

## Decision Card

- Required: Yes
- Rationale: The stage approves or rejects biological interpretations and
  public-facing claims after a material correction to the experimental unit.
- Approval status: Approved and applied on 2026-07-30

## Decision

Implement a pre-specified evidence assessment of the corrected results before
replacing any curated assets, changing public claims, or integrating these
findings with future pseudobulk results. The automated assessment may
prioritize rows for review but cannot assign a final disposition or approve a
biological claim. No corrected biological claim is approved at the current
review stage.

On 2026-07-30, the repository owner approved a conservative public scope, the
recommended dispositions for all 13 prioritized evidence rows, professionally
edited README wording, three corrected figures and captions, and removal of
the four provisional curated assets.

## Rationale And Tradeoffs

A structured review delays public promotion but prevents statistically ranked
outputs from being converted directly into biological conclusions. Reviewing
all major evidence dimensions in one table reduces cherry-picking and creates
an auditable boundary between analysis and interpretation. External validation
would provide stronger evidence but is a later stage and should not be implied
by this internal review.

## Interview Question

- Question: How would you decide whether corrected single-cell age
  associations are ready for public presentation after changing the
  experimental unit?
- Skill tested: Scientific result triage, effect-size interpretation,
  multiple-testing awareness, sensitivity analysis, confounding assessment,
  and cautious communication.
- Expected answer should mention: Use only the corrected replicate-level
  analysis; predefine evidence criteria; review effect size, confidence
  intervals, FDR, replicate support, age coverage, and sensitivity together;
  inspect confounding and annotation uncertainty; distinguish association from
  causation; require traceable human sign-off before replacing claims.
- Status: Ready

## Follow-up Biological Question

- Question: Which corrected immune-aging associations reproduce in an
  independent donor cohort with comparable cell populations and metadata?
- Biological theme: External validity and reproducibility of immune-aging
  associations.
- Why it matters: Internal significance and sensitivity do not establish that
  an association generalizes beyond the GSE164378 analysis cohort.
- Possible next analysis: Apply the frozen replicate-aware analysis contract
  to an independent PBMC aging cohort and compare effect directions,
  uncertainty intervals, supported cell populations, and prediction error.
- Status: Open

## No-AI Survival Exercise

- Estimated time: 20-30 minutes
- Skill tested: Evidence-based triage of corrected computational-biology
  results.
- Allowed resources: This decision card, the corrected composition and
  signature association tables, age-prediction metrics, sensitivity summary,
  and comparison table; no generative AI.
- Expected output: A one-page evidence table for three selected results,
  including effect, interval, FDR, replicate support, sensitivity, main
  limitation, and a retain/exploratory/exclude disposition.
- Pass criteria: Uses corrected results only for the disposition, does not rank
  by FDR alone, identifies at least one confounding or support limitation,
  distinguishes statistical observation from biological interpretation, and
  avoids causal or clinical language.
- Follow-up biological question: What additional cohort or assay evidence
  would be needed to strengthen each retained interpretation?
- Status: Not started

## Implementation Evidence

- Review module: `src/review_corrected_results.py`.
- Configured screening contract:
  `config/config.real.full.replicate_corrected.yml`.
- Workflow target: `review_corrected_results` in
  `workflows/Snakefile.replicate_corrected`.
- Machine-readable evidence:
  `results/gse164378_full_replicate_corrected/review/evidence_table.csv`.
- Review summary and promotion manifest:
  `results/gse164378_full_replicate_corrected/review/`.
- Human-readable report: `docs/corrected_results_review.md`.
- Automated screen result: 109 evidence rows, including 12 association
  candidates for human review, 10 sensitivity-limited associations, 86
  non-prioritized associations, and one internal-prediction evidence row.
- Final state: 13 prioritized dispositions are recorded; six retained signals
  are approved for cautious public presentation, six remain exploratory, and
  one is excluded. Three corrected assets are approved and the signature
  heatmap is not selected for public curation.
- Validation: focused tests passed; all 108 association rows and
  their sensitivity records use `biological_replicate_id`; no comparison
  status is missing. The targeted Snakemake dry-run completed with all review
  outputs present and up to date.

The minimum-replicate threshold is an inclusion floor for review, not evidence
that a cell population is biologically robust. Cell-count distribution,
annotation quality, age/batch/sex structure, plausibility, and wording still
require manual assessment for each candidate.

## Human Disposition Record

- Reviewer: Repository owner
- Review date: 2026-07-30
- Public scope: Conservative core
- Tracked source:
  `docs/reviews/corrected_results_dispositions.yml`
- Retain for cautious drafting:
  - Tcm/Naive cytotoxic T-cell composition, with a combined-label caveat.
  - CD16+ NK-cell composition as a core result.
  - MAIT-cell, Tem/Effector helper T-cell, and non-classical monocyte
    composition as secondary results.
  - Age prediction as technical internal-CV evidence only.
- Exploratory:
  - Regulatory T-cell composition.
  - Five prioritized signature-score associations.
- Exclude:
  - `CD8a/a` composition because only 144 cells contributed in total and the
    median support was one cell per replicate.
- Public claims approved: Yes, for the six retained signals and exact approved
  wording only.
- Asset replacements approved: Yes, for the three selected corrected figures.

## Draft Asset Decisions

- Decision date: 2026-07-30
- Composition figure: Three panels showing Tcm/Naive cytotoxic T cells,
  CD16+ NK cells, and MAIT cells.
- Signature heatmap: Do not select for curated public assets; retain signature
  associations as exploratory generated results.
- Age-prediction figure: Generate a clearer internal-CV figure using
  descriptive model names and no redundant rankings.
- Forest plot: Show all five retained composition associations and no
  exploratory signatures.
- Draft outputs:
  `results/gse164378_full_replicate_corrected/review/figures/`.
- Proposed wording and captions:
  `docs/reviews/corrected_public_promotion_record.md`.
- Public asset replacements approved: Yes
- README wording approved: Yes, after revision for natural professional
  scientific language.
- Provisional asset action: Removed from `docs/assets/`; generated provisional
  results remain available in their original local result root.
- Final curated assets:
  - `docs/assets/gse164378_corrected_composition_core_trends.png`
  - `docs/assets/gse164378_corrected_composition_effect_forest.png`
  - `docs/assets/gse164378_corrected_age_prediction_internal_cv.png`
