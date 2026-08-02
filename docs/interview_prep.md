# Interview Prep Notes

This is a companion to `docs/repo_retrospective.md`, especially
`## 8. Parts To Understand For Interviews`. It does not replace the
retrospective or duplicate its full audit. Instead, it expands the interview
topics into concise explanations, Q&A prompts, decision rationale, tradeoffs,
and points to review before discussing this project.

## Source Notes

No formal decision-card or run-log files were found in the maintained
repository. The decision context below is drawn from:

- `docs/repo_retrospective.md`
- `IMMUNE-AGING_NOTES.md`
- `README.md`
- `docs/phase2_next_steps_checklist.md` (local roadmap; ignored by `.gitignore`)
- `workflows/Snakefile`
- `config/config.demo.yaml`, `config/config.real.yml`,
  `config/config.real.full.yml`
- Source modules under `src/`

When a point is inferred from code/configuration rather than stated directly in
project documentation, it is marked as inferred.

## Core Project Pitch

Short version:

> This is a reproducible Snakemake pipeline for PBMC single-cell RNA-seq
> immune-aging analysis. It takes demo or real `.h5ad` inputs through metadata
> integration, QC, optional doublet filtering, scVI modeling, clustering,
> CellTypist annotation, donor-level age-composition analysis, signature-age
> analysis, age prediction from scVI embeddings, and sensitivity/supplementary
> reporting.

What to emphasize:

- The project is not just a notebook; it is a modular workflow.
- The biological focus is immune aging in PBMCs.
- The main statistical guardrail is donor-level analysis for donor-level
  covariates such as age.
- The outputs are tables and figures that can be reviewed, reproduced, and
  discussed.
- The current scope is Phase 1; donor-aware pseudobulk differential expression
  remains the main scientific extension.

## What I Should Be Able To Explain

### 1. What Problem The Workflow Solves

Good answer:

> The workflow turns a PBMC single-cell dataset into interpretable immune-aging
> outputs. It standardizes or ingests data, integrates donor metadata, builds a
> scVI representation, annotates immune cells, and then asks whether cell-type
> composition, predefined gene-program signatures, or latent embeddings are
> associated with donor age.

Key details to know:

- Input is either demo PBMC3K or a custom `.h5ad`.
- Real-data profiles use local raw data paths under `data/raw/`.
- Generated outputs go to `results/` or profile-specific subdirectories.
- Curated public figures are copied to `docs/assets/`.

Common interview follow-up:

> Why not just use a notebook?

Answer:

> A notebook is useful for exploration, but this project is structured as a
> reproducible workflow. Snakemake records dependencies between steps, YAML
> configs separate parameters from code, and each analysis stage writes explicit
> outputs. That makes the project easier to rerun, review, and scale.

### 2. Why Donor-Level Inference Matters

Good answer:

> Age is measured at the donor level, not independently for every cell. If I
> treated cells as independent biological replicates, I would inflate the sample
> size and overstate evidence. The pipeline therefore aggregates or tests at the
> donor or donor-cell-type level for age-related analyses.

Where this appears:

- `src/composition_age.py`: donor-level cell-type fractions.
- `src/signature_age.py`: donor-cell-type signature summaries.
- `src/age_prediction.py`: donor-cell-type latent aggregation and
  donor-grouped cross-validation.

Tradeoff:

- Donor-level analysis is statistically more appropriate for age effects.
- It reduces apparent sample size compared with per-cell tests.
- Sparse donor/cell-type combinations need thresholds to avoid unstable
  estimates.

### 3. How Metadata Enters The Pipeline

Good answer:

> Metadata is integrated in the workflow, not as a manual side step. The
> metadata stage can parse sample, donor, batch, and barcode information from
> `obs_names`, then merge an external metadata table using configured join keys.
> It prefers table-derived donor and batch values when available and records
> source/conflict fields.

Where this appears:

- `src/metadata_integrate.py`
- `config/config.real.yml`
- `config/config.real.full.yml`
- `docs/real_data.md`

Decision rationale:

- Observed: `IMMUNE-AGING_NOTES.md` says real-data metadata integration happens
  inside the workflow.
- Inferred: keeping metadata integration in Snakemake reduces hidden manual
  state and makes age/sex/batch assumptions easier to audit.

Tradeoff:

- More upfront code and config complexity.
- Better provenance and fewer undocumented spreadsheet edits.

What to be ready to explain:

- Join keys: `obs_join_key`, `table_join_key`.
- Required age-analysis fields: age, donor ID, cell type.
- Why batch and sex can be covariates.
- Why metadata conflicts are important.

### 4. What scVI Contributes

Good answer:

> scVI learns a lower-dimensional latent representation of single-cell
> expression. In this workflow, that representation is stored as `X_scVI` and
> used for neighbor graph construction, UMAP/clustering, and age-prediction
> features.

Where this appears:

- `src/scvi_train.py`
- `src/cluster.py`
- `src/age_prediction.py`

Decision rationale:

- Observed: the project briefing names `scanpy` + `scvi-tools` as the core
  stack.
- Inferred: scVI is used to produce a compact representation that can absorb
  batch/covariate structure better than raw expression features for downstream
  embedding-based analyses.

Tradeoff:

- scVI is powerful and modern, but adds model-training complexity.
- It needs careful environment handling, especially CPU versus GPU.
- The latent space is useful for prediction and visualization, but it is not a
  direct gene-level differential expression result.

### 5. What CellTypist Contributes

Good answer:

> CellTypist provides automated immune-cell label transfer. The workflow uses it
> after clustering/scVI to assign `cell_type`, optionally in chunks for larger
> datasets, and stores confidence values when available.

Where this appears:

- `src/annotate_celltypist.py`
- `config/config.real.yml`
- `config/config.real.full.yml`

Tradeoff:

- Automated annotation is reproducible and fast.
- It depends on the selected CellTypist model and may need biological review.
- In large runs, chunking reduces memory pressure but disables majority voting
  to avoid chunk-wise voting artifacts.

### 6. What Composition Analysis Measures

Good answer:

> The composition analysis asks whether donor-level fractions of each annotated
> cell type change with age. It computes fractions per donor, applies minimum
> donor/cell thresholds, can adjust for covariates like sex and batch, reports
> Spearman correlations, slopes, bootstrapped confidence intervals, and FDR.

Where this appears:

- `src/composition_age.py`
- Outputs:
  - `age_celltype_fraction_by_donor.csv`
  - `age_celltype_trend_stats.csv`
  - `age_celltype_top_trends.png`

Tradeoff:

- This directly answers abundance-shift questions.
- It does not identify within-cell-type gene-expression changes.
- It depends on correct cell-type annotation and donor metadata.

### 7. What Signature Analysis Measures

Good answer:

> Signature analysis computes mean expression scores for predefined gene sets
> and tests whether those scores vary with age within donor-cell-type groups.
> It is a targeted gene-program analysis, not genome-wide differential
> expression.

Where this appears:

- `src/signature_age.py`
- Outputs:
  - `signature_scores_by_donor_celltype.csv`
  - `signature_age_associations.csv`
  - `signature_gene_coverage.csv`
  - `signature_age_heatmap.png`

Decision rationale:

- Observed: the project goal includes age-linked programs/signatures by cell
  type.
- Inferred: predefined signatures provide interpretable, lower-dimensional
  biology before implementing genome-wide pseudobulk DE.

Tradeoff:

- Easier to interpret and visualize than genome-wide DE.
- Limited to the selected signatures.
- Does not discover new genes outside the predefined sets.

### 8. What Age Prediction Evaluates

Good answer:

> The age-prediction stage asks whether donor age is recoverable from scVI
> latent embeddings aggregated by donor and cell type. It compares candidate
> regression models under biological-replicate-grouped cross-validation and
> reports prediction tables, metrics, model comparisons, and plots. The
> existing `X_scVI` representation was trained once on the full cohort, so this
> is a transductive internal evaluation rather than end-to-end validation on
> unseen donors.

Where this appears:

- `src/age_prediction.py`
- Outputs:
  - `age_pred_cv_predictions.csv`
  - `age_pred_metrics.csv`
  - `age_pred_model_comparison_summary.csv`
  - `age_pred_observed_vs_predicted.png`

What to emphasize:

- Biological-replicate-grouped CV reduces leakage across replicate-level
  labels.
- The result is a predictive signal assessment, not a causal model or a fully
  inductive unseen-donor evaluation.
- The corrected internal analysis reports 12.39 years donor-level MAE
  (95% interval 11.63-13.15) across 316 biological replicates, versus
  15.13 years for a fold-specific mean-age baseline. These metrics remain
  internal evidence rather than external validation.

Tradeoff:

- Predictive performance is useful for showing age signal in the latent space.
- It is less directly interpretable than composition or signature analyses.
- Prediction can be sensitive to donor distribution, model class, and
  cross-validation design.

### 9. What Sensitivity Analysis Does

Good answer:

> Sensitivity analysis reruns selected composition and signature analyses under
> alternative settings, such as covariate adjustment, sparse age-bin handling,
> and minimum cell thresholds. It helps check whether major conclusions depend
> on one arbitrary parameter choice.

Where this appears:

- `src/sensitivity_age.py`
- `src/supplementary_age_plots.py`
- `config/config.real.full.yml`

Tradeoff:

- Improves robustness and reviewer confidence.
- Adds runtime and more outputs to interpret.
- It is not a substitute for external validation.

### 10. What Remains Incomplete

Good answer:

> The main missing scientific layer is donor-aware, cell-type-specific
> pseudobulk differential expression. Current signature analyses test selected
> programs, but they do not provide genome-wide within-cell-type transcriptional
> change analysis. The real-data acquisition path is also not fully productized.

Where this appears:

- `IMMUNE-AGING_NOTES.md`
- `docs/phase2_next_steps_checklist.md`
- `README.md`
- `docs/repo_retrospective.md`

Tradeoff:

- Freezing Phase 1 keeps the repo coherent and reviewable.
- Adding pseudobulk DE would make the biology stronger but expands scope and
  requires careful design.

## Decision Rationale And Tradeoffs

| Decision | Source | Rationale | Tradeoff |
| --- | --- | --- | --- |
| Use Snakemake + YAML configs instead of a monolithic script | Observed in `IMMUNE-AGING_NOTES.md`, `workflows/Snakefile`, `config/` | Reproducible workflow graph and parameterized runs. | More files and structure to maintain. |
| Support demo, pilot, and scaled/full profiles | Observed in configs and briefing | Enables smoke testing, smaller real-data runs, and larger analysis profiles. | Multiple configs can drift if not audited. |
| Keep raw data and generated results out of version control | Observed in `.gitignore`, README, briefing | Avoids committing large/private/generated artifacts while preserving reproducibility through code/config. | Reviewers need instructions or curated assets to inspect results. |
| CPU-first environment with GPU opt-in | Observed in `environment.yml`, configs, README | Portable default across machines; GPU remains possible after local CUDA validation. | CPU runs can be slow for scVI and large data. |
| Merge metadata inside the workflow | Observed in briefing and `src/metadata_integrate.py` | Keeps donor age/batch/sex integration reproducible and auditable. | Requires careful join-key configuration. |
| Use scVI latent embeddings | Observed in stack and workflow | Compact representation for clustering and prediction. | Adds model complexity and less direct gene-level interpretability. |
| Use CellTypist for annotation | Observed in stack and `src/annotate_celltypist.py` | Reproducible automated immune-cell labeling. | Labels still require biological sanity checks. |
| Use donor-level composition analysis | Observed in `src/composition_age.py` | Appropriate for donor-level age effects. | Fewer independent units than per-cell analysis. |
| Use predefined signature scoring | Observed in `src/signature_age.py` and briefing | Interpretable targeted program analysis. | Not genome-wide discovery. |
| Use biological-replicate-grouped CV for age prediction | Observed in `src/age_prediction.py` and corrected validation outputs | Prevents reuse of replicate-level labels across training and test folds. | More stringent and may reduce apparent performance. |
| Screen corrected results before public promotion | Observed in `src/review_corrected_results.py` and decision card 003 | Keeps automated ranking separate from biological interpretation and requires traceable human sign-off. | Delays promotion and still requires manual review of confounding, support, and wording. |
| Use donor-aware pseudobulk DE for gene-level age effects | Observed in decision card 001 and the corrected workflow | Preserves the biological replicate and supports count-based cell-type-specific testing. | Technical acceptance does not substitute for biological review or external validation. |

## Interview Q&A Prompts

### Q: What is the biological question?

A: The project asks how the PBMC immune landscape changes with donor age. It
looks at cell-type abundance shifts, selected age-linked signature changes, and
whether donor age is predictable from scVI latent embeddings.

### Q: What makes this more than a standard Scanpy notebook?

A: The workflow is modular and reproducible: Snakemake defines dependencies,
YAML configs define run profiles, source modules write structured outputs, and
large raw/generated files are excluded from version control while curated
figures are stored under `docs/assets/`.

### Q: Why is donor-level analysis important?

A: Donor age is not a cell-level measurement. Treating cells as independent
replicates would inflate evidence. The pipeline therefore aggregates or
evaluates at donor or donor-cell-type level for age-related analyses.

### Q: How do you avoid data leakage in age prediction?

A: The age-prediction module aggregates latent embeddings by biological
replicate and cell type and uses replicate-grouped cross-validation for the
regressor. That prevents replicate-level labels from crossing regressor folds.
The current `X_scVI` embedding was learned on the full cohort, however, so a
strict end-to-end generalization claim would require fold-specific scVI
training or an independent cohort.

### Q: What did scVI add here?

A: scVI provides a learned latent representation used for downstream neighbor
graphs, UMAP/clustering, and age-prediction features. It is a representation
learning step, not the final biological interpretation by itself.

### Q: What did CellTypist add here?

A: CellTypist converts clusters/cells into immune-cell labels in a reproducible
way. Those labels make the composition, signature, and age-prediction outputs
cell-type aware.

### Q: What is the difference between signature analysis and differential expression?

A: Signature analysis tests predefined gene programs by summarizing expression
of selected genes. Differential expression tests genes individually. This
repository now uses biological-replicate-aware pseudobulk profiles and edgeR by
cell type; the full outputs passed technical validation, but their gene-level
biological interpretation remains pending.

### Q: What are the most important results to discuss cautiously?

A: After correcting the biological-replicate identifier, the public summary
focuses on lower Tcm/Naive cytotoxic T-cell and MAIT-cell fractions and a
higher CD16+ NK-cell fraction at older ages. All three associations were stable
across five sensitivity scenarios, but they remain cross-sectional and may
still be affected by residual confounding or annotation uncertainty. The
age-prediction model achieved 12.39 years donor-level MAE across 316 biological
replicates, compared with 15.13 years for a fold-specific mean-age baseline.
That result reflects modest internal predictive signal, not an externally
validated clock or clinical biomarker. Signature-score findings remain
exploratory.

### Q: What would you improve next?

A: Review the technically accepted pseudobulk gene-level results under a
pre-specified interpretation protocol, productize real-data acquisition and
standardization, and validate retained findings in an independent or
integrated PBMC aging atlas. A stricter predictive evaluation would retrain the
representation within each outer fold or use an independent cohort.

## Red Flags To Avoid In Interviews

- Do not claim that the workflow proves causal aging mechanisms.
- Do not describe signature scores as genome-wide differential expression.
- Do not imply that every cell is an independent age replicate.
- Do not overstate the GSE164378 example as external validation.
- Do not present technically accepted pseudobulk output as biologically
  interpreted or externally validated evidence.
- Do not imply raw data or generated results are meant to be committed.
- Do not claim the full profile always means every available cell; the current
  config uses `max_cells: 1000000`.

## Fast Review Checklist

Before an interview, review:

- `README.md`: project pitch, workflow, validation command, example results.
- `docs/repo_retrospective.md`: sections 2, 3, 6, 7, 8, and 10.
- `workflows/Snakefile`: rule order and output contract.
- `config/config.real.yml` and `config/config.real.full.yml`: real-data
  assumptions and thresholds.
- `src/metadata_integrate.py`: metadata join and conflict/source logic.
- `src/composition_age.py`: donor-level fraction analysis.
- `src/signature_age.py`: signature scoring versus pseudobulk DE distinction.
- `src/age_prediction.py`: donor-grouped CV and model comparison.
- `docs/phase2_next_steps_checklist.md`: next-step roadmap and current gaps.

## One-Minute Closing Summary

> This project demonstrates how I would turn a single-cell immune-aging idea
> into a reproducible computational biology workflow. The strongest parts are
> the Snakemake/config design, metadata integration, memory-aware `.h5ad`
> handling, scVI/CellTypist integration, donor-level age analyses, leakage-aware
> age prediction, and sensitivity reporting. The main scientific gap I would
> address next is donor-aware pseudobulk differential expression by cell type,
> followed by cleaner real-data acquisition and cross-dataset validation.
