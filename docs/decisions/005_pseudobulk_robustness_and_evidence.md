# Decision Card: Pseudobulk Robustness Audit And Evidence Preparation

## Metadata

- Stage: Pseudobulk robustness audit and evidence preparation
- Priority: High
- Status: Implemented and technically validated - D-class review pending
- Risk classification (A/B/C/D): C - scientifically important

## Context

The corrected GSE164378 pseudobulk analysis passed its approved technical
contract: raw counts were summed by biological replicate and cell type, all 13
models were full rank, and edgeR produced 107,948 gene-by-cell-type tests.
Technical acceptance does not establish that individual effects are stable to
sample-support choices, dependent on a single batch, or ready for biological
interpretation.

The existing outputs contain many globally significant rows. Ranking those
rows by FDR alone would be vulnerable to post hoc selection and would not
address effect stability, uneven profile support, the exploratory NK-cell
tier, or visual model diagnostics. This stage therefore prepares traceable
evidence for later human review without approving genes, pathways, or public
claims.

## Decision

Implement a post-model robustness audit against the immutable validated
pseudobulk inputs and primary result table.

1. Keep `~ sex + batch + age_decade` as the sole primary model. Sensitivity
   results must never replace it silently.
2. Refit the same edgeR quasi-likelihood model after requiring at least 100
   cells per pseudobulk profile. This tests dependence on profiles near the
   approved 50-cell inclusion floor.
3. Fit covariate-omission sensitivities without sex and without batch. These
   are confounding diagnostics, not preferred alternative models.
4. Fit leave-one-batch-out sensitivities when the remaining design still meets
   the approved replicate, age-span, rank, and residual-df requirements.
5. Limit gene-level sensitivity output to candidates that passed global BH FDR
   below 0.05 in the primary analysis. Preserve all candidate rows; do not use
   sensitivity significance as a second discovery filter.
6. Report direction concordance, effect ranges, maximum absolute effect
   change, and the number of estimable sensitivity scenarios. Do not assign an
   automated `robust` or `validated` label.
7. Produce a compact review queue by taking the configured number of strongest
   primary candidates in each direction and cell type. Selection is for review
   workload only and does not approve interpretation.
8. Keep the exploratory NK-cell results explicitly separated from the 12
   primary populations.
9. Record profile support, design diagnostics, runtime provenance, and visual
   review of every edgeR mean-difference diagnostic page.
10. Defer gene annotation review, pathway enrichment, biological disposition,
    and public wording to a later D-class human-review stage.

Expected generated outputs are a scenario manifest, scenario diagnostics,
candidate-level sensitivity results, a cell-type robustness summary, a full
candidate evidence table, a compact review queue, and a machine-readable audit
summary under the corrected pseudobulk result root. A tracked report may
describe the audit contract and observed technical limitations but must not
make gene-level biological claims.

## Rationale And Tradeoffs

The selected sensitivities target the main observable threats in this design:
profiles close to the cell-count floor, adjustment dependence, and influence
of individual processing batches. Reusing the frozen pseudobulk counts avoids
changing cell annotation or aggregation while retaining the same count-model
engine and runtime.

Leave-one-replicate-out refitting for every gene would be substantially more
expensive and would encourage unstable per-gene influence thresholds. The
leave-one-batch-out analysis is a more focused diagnostic for the documented
age-batch concern. Sensitivity p-values are not treated as independent
confirmation because the scenarios reuse the same cohort and candidates were
selected from the primary analysis.

The audit can identify unstable or weakly supported candidates, but it cannot
establish biological validity, external replication, causal aging effects, or
clinical relevance.

## Interview Question

- Question: How would you audit robustness of donor-level pseudobulk differential-expression results before interpreting the top genes?
- Skill tested: Sensitivity design, confounding assessment, experimental-unit preservation, multiple-testing awareness, and evidence triage.
- Expected answer should mention: Preserve the primary donor-level count model; vary one support or adjustment assumption at a time; check estimability; assess effect direction and magnitude rather than demanding repeated significance; inspect batch influence and model diagnostics; keep exploratory populations separate; and require human review before biological claims.
- Status: Ready

## Follow-up Biological Question

- Question: Which cell-type-specific age-associated genes retain comparable effect direction and magnitude across support and batch-sensitivity analyses?
- Biological theme: Robustness of cell-type-specific transcriptional immune-aging associations.
- Why it matters: Stable statistical evidence is a prerequisite for deciding which genes merit annotation review and pathway-level follow-up.
- Possible next analysis: Human-review the compact evidence queue, verify gene annotation and known technical artifacts, then run direction-aware enrichment under a separately approved interpretation protocol.
- Status: Open

## No-AI Survival Exercise

- Estimated time: 20-30 minutes
- Skill tested: Reading a pseudobulk sensitivity table and distinguishing statistical stability from biological validation.
- Allowed resources: This decision card, the generated candidate evidence table and scenario diagnostics, a spreadsheet or short local script, and package documentation; no generative AI.
- Expected output: A one-page comparison of three candidate genes across the primary, high-cell-support, and leave-one-batch-out scenarios, including effect direction, effect range, estimability, population tier, and one limitation per candidate.
- Pass criteria: Keeps the biological replicate as the unit, does not require repeated FDR significance as proof, identifies direction or effect instability correctly, treats NK cells as exploratory, and avoids causal or pathway claims.
- Follow-up biological question: What external cohort or orthogonal assay would be most informative for each candidate after internal robustness review?
- Status: Not started

## Implementation Evidence

- Completed on 2026-08-02 against the immutable, technically validated
  corrected-replicate pseudobulk inputs.
- `src/pseudobulk_robustness.R` fitted the approved high-cell-support,
  covariate-omission, and leave-one-batch-out scenarios in the pinned
  R 4.6.1/Bioconductor 3.23/edgeR 4.10.1 runtime.
- The full audit recorded 215 cell-type-by-scenario models: 212 completed,
  three were explicitly non-estimable, and none failed during model fitting.
- The evidence stage preserved all 7,970 global-FDR candidate rows from the
  primary analysis: 7,961 primary and nine exploratory NK-cell rows.
- The compact queue contains 249 rows and remains bounded at ten candidates
  per cell type and effect direction. Every candidate retains
  `requires_human_review=true` and `automated_disposition=not_assigned`.
- All 13 primary edgeR diagnostic pages were independently rendered and
  inspected. Follow-up flags remain for MAIT, Memory B, Naive B, and
  exploratory NK populations; these flags do not assign gene-level meaning.
- The Windows Docker client remained open after the container had exited and
  written complete outputs. The two identified host processes were stopped;
  output acceptance used independent schema, key, count, scenario, and model
  status checks. The wrapper now decodes subprocess output explicitly as
  UTF-8, and this behavior has a regression test.
- Generated evidence is under
  `results/gse164378_full_replicate_corrected/pseudobulk_de/robustness/`.
  The tracked technical report is `docs/pseudobulk_robustness_audit.md`.
- Gene annotation, enrichment, biological disposition, and public claims were
  not performed and still require a separately approved D-class stage.
