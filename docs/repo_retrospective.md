# Repository Retrospective

This audit is based on the repository files, `IMMUNE-AGING_NOTES.md`,
`README.md`, `workflows/Snakefile`, configuration files, and source modules as
present in this working tree on 2026-06-28. It does not assume motivations not
documented in the repository. When a statement is inferred from code or local
outputs rather than explicitly documented, it is marked as inferred.

## 1. Project Purpose In Plain Language

This repository builds a reproducible analysis workflow for single-cell RNA-seq
data from peripheral blood mononuclear cells (PBMCs). In plain language, it
takes raw or pre-standardized single-cell data, cleans and annotates cells,
learns a compact scVI representation, and then asks whether immune-cell
composition, selected gene-program scores, and scVI latent embeddings are
associated with donor age.

The project is also structured as an independently reviewable computational
biology repository: it uses Snakemake, YAML configuration files, structured
Python modules, versioned documentation assets, and generated tables and
figures rather than a single manual notebook.

## 2. Main Biological/Longevity Question

The main biological question is: how does the PBMC immune landscape change with
donor age?

The repository breaks that into three practical questions documented in
`IMMUNE-AGING_NOTES.md`:

- Which immune cell types or states shift with age?
- Which age-linked programs or signatures change by cell type?
- Can donor age be predicted from scVI latent embeddings?

The longevity framing is therefore immune aging in human PBMCs. Observed: the
reviewed repository documents do not currently claim causal mechanisms,
lifespan effects, intervention targets, or validated biomarkers.

## 3. Workflow Map

Inferred from `workflows/Snakefile`, module command-line interfaces, and config
profiles: the workflow is defined in `workflows/Snakefile` and parameterized by
config profiles under `config/`.

| Stage | Input | Output | Purpose |
| --- | --- | --- | --- |
| `ingest` | Demo PBMC3K data or configured custom `.h5ad` | `results/.../01_raw.h5ad` | Load demo data or link/copy/subsample real data. |
| `metadata` | `01_raw.h5ad` plus optional metadata table | `01_meta.h5ad` | Parse cell/sample/donor fields and merge donor metadata such as age, sex, and batch. |
| `qc` | `01_meta.h5ad` | `02_qc.h5ad` | Apply gene/cell QC filters and preserve raw counts in a sparse counts layer when enabled. |
| `doublets` | `02_qc.h5ad` | `03_nodoublets.h5ad` | Optionally run scVI/SOLO doublet detection. |
| `scvi_train` | `03_nodoublets.h5ad` | `04_scvi.h5ad`, `results/.../models/scvi_model/` | Train scVI and store latent representation in `obsm["X_scVI"]`. |
| `cluster` | `04_scvi.h5ad` | `05_clustered.h5ad` | Build neighbors on scVI latent space, compute UMAP, and run Leiden clustering. |
| `annotate` | `05_clustered.h5ad` | `06_annotated.h5ad` | Assign immune-cell labels using CellTypist. |
| `report` | `06_annotated.h5ad` | UMAP, QC, composition figures and summary tables | Produce baseline descriptive outputs. |
| `composition_age` | `06_annotated.h5ad` | Age-composition figures and tables | Test donor-level cell-type fraction associations with age. |
| `signature_age` | `06_annotated.h5ad` | Signature-score figures and tables | Score configured gene signatures and test donor-cell-type associations with age. |
| `age_prediction` | `06_annotated.h5ad` | Prediction figures, CV predictions, metrics, model summary | Aggregate full-cohort scVI latent embeddings by biological replicate and cell type, then evaluate regressors with replicate-grouped CV; this is transductive internal evaluation. |
| `sensitivity_age` | `06_annotated.h5ad` | `results/.../sensitivity_age/*` | Re-run selected age-analysis settings to assess robustness when enabled. |
| `supplementary_age` | Age-analysis tables and optional sensitivity manifest | Supplementary forest/stability figures | Summarize effect sizes and sensitivity stability. |

Observed configs:

- `config/config.demo.yaml`: demo or smoke-test profile.
- `config/config.real.yml`: pilot real-data profile with `max_cells: 50000`.
- `config/config.real.full.yml`: scaled/full profile with output directory
  `results/gse164378_full` and `max_cells: 1000000`.
- CPU is the tracked default; GPU is documented as an opt-in after CUDA-enabled
  PyTorch validation.

## 4. Main Files And What They Do

File roles below are observed from repository structure and inferred from source
code where the role is not explicitly described in `README.md` or
`IMMUNE-AGING_NOTES.md`.

| File | Role |
| --- | --- |
| `README.md` | Public-facing project summary, quickstart commands, validation commands, example GSE164378 result narrative, and output inventory. |
| `IMMUNE-AGING_NOTES.md` | Project briefing: goals, constraints, important commands, files, results, and current unresolved gaps. |
| `AGENTS.md` | Repo-specific working rules for computational biology, reproducibility, raw-data policy, workflow control, and scientific caution. |
| `environment.yml` | CPU-first conda environment with Python, Snakemake, Scanpy, scVI tools, CellTypist, and scientific Python dependencies. |
| `.gitignore` | Excludes raw data, generated results, `.h5ad` files, Python caches, `.snakemake/`, and local working documents; allows selected images under `docs/assets/`. |
| `workflows/Snakefile` | Snakemake workflow contract and output graph. |
| `config/config.demo.yaml` | Small demo profile using PBMC3K. |
| `config/config.real.yml` | Pilot real-data profile. |
| `config/config.real.full.yml` | Scaled/full real-data profile. |
| `docs/real_data.md` | Real-data schema, memory-safety notes, metadata requirements, and optional Synapse metadata discovery notes. |
| `docs/troubleshooting.md` | GPU and memory troubleshooting notes. |
| `docs/phase2_next_steps_checklist.md` | Local roadmap for Phase 2 work; this file is ignored by `.gitignore`. |
| `src/utils.py` | Shared config loading and directory creation helpers. |
| `src/ingest.py` | Loads demo data or custom `.h5ad`; supports backed subsampling for large data. |
| `src/metadata_integrate.py` | Parses obs names, normalizes metadata columns, merges external metadata, and records source/conflict fields. |
| `src/qc.py` | Computes QC metrics, filters cells/genes, ensures sparse matrix representation, and stores counts layer when enabled. |
| `src/doublets_solo.py` | Optional SOLO doublet detection using scVI. |
| `src/scvi_train.py` | Trains scVI, stores latent representation, and saves model weights. |
| `src/cluster.py` | Builds neighbors, UMAP, and Leiden clusters from the configured representation. |
| `src/annotate_celltypist.py` | Runs CellTypist annotation, with chunking support for larger datasets. |
| `src/report.py` | Generates UMAPs, QC plots, cell-type composition plots, and baseline count/fraction tables. |
| `src/composition_age.py` | Donor-level cell-type fraction versus age analysis with optional covariate adjustment, bootstrapped confidence intervals, and FDR correction. |
| `src/signature_age.py` | Gene-signature scoring by donor/cell type with sampling, covariate adjustment, bootstrapped intervals, and FDR correction. |
| `src/age_prediction.py` | Biological-replicate/cell-type aggregation of full-cohort scVI latent embeddings, replicate-grouped regressor cross-validation, model comparison, bootstrap confidence intervals, and prediction plots. |
| `src/sensitivity_age.py` | Runs predefined sensitivity scenarios for composition and signature analyses. |
| `src/supplementary_age_plots.py` | Builds supplementary effect-size and sensitivity-stability figures. |
| `src/fetch_synapse_metadata.py` | Optional helper for credentialed Synapse metadata discovery/download. |
| `src/build_donor_metadata.py` | Helper to build a replicate-level audit table from an explicit authoritative replicate field and validated cell-level join. |
| `src/plot_style.py` | Shared plotting style and placeholder-figure helpers. |

## 5. Main Outputs Generated

Inferred from `workflows/Snakefile`, README output lists, and local result
directories: the workflow writes outputs under the configured `project.out_dir`,
usually `results/` for demo runs, `results/gse164378_pilot` for pilot real-data
runs, or `results/gse164378_full` for the full-profile example.

Main intermediate objects:

- `01_raw.h5ad`
- `01_meta.h5ad`
- `02_qc.h5ad`
- `03_nodoublets.h5ad`
- `04_scvi.h5ad`
- `05_clustered.h5ad`
- `06_annotated.h5ad`
- `models/scvi_model/`

Baseline report outputs:

- `figures/umap_leiden.png`
- `figures/umap_cell_type.png`
- `figures/qc_distributions.png`
- `figures/cell_type_composition.png`
- `tables/cell_type_counts.csv`
- `tables/cell_type_fractions.csv`
- `tables/cluster_celltype_crosstab.csv`

Age-composition outputs:

- `figures/age_celltype_composition_by_bin.png`
- `figures/age_celltype_top_trends.png`
- `tables/age_celltype_fraction_by_donor.csv`
- `tables/age_celltype_trend_stats.csv`

Signature-age outputs:

- `figures/signature_age_heatmap.png`
- `figures/signature_age_top_associations.png`
- `tables/signature_scores_by_donor_celltype.csv`
- `tables/signature_age_associations.csv`
- `tables/signature_gene_coverage.csv`

Age-prediction outputs:

- `figures/age_pred_observed_vs_predicted.png`
- `figures/age_pred_mae_by_celltype.png`
- `tables/age_pred_cv_predictions.csv`
- `tables/age_pred_metrics.csv`
- `tables/age_pred_model_comparison_summary.csv`

Sensitivity and supplementary outputs:

- `sensitivity_age/sensitivity_manifest.csv`
- `sensitivity_age/sensitivity_summary.csv`
- `sensitivity_age/done.txt`
- `figures/supp_age_effect_ci_forest.png`
- `figures/supp_age_sensitivity_stability.png`

Curated README assets:

- `docs/assets/gse164378_corrected_composition_core_trends.png`
- `docs/assets/gse164378_corrected_composition_effect_forest.png`
- `docs/assets/gse164378_corrected_age_prediction_internal_cv.png`

These assets come from the biological-replicate-corrected review. The signature
heatmap remains an exploratory generated result and is not a curated README
asset.

## 6. Boilerplate Versus Scientifically Important Stages

The classification below is inferred from the workflow structure and source
modules. It is an audit judgment, not a stated repository requirement.

Boilerplate here does not mean disposable. It means the stage is mostly
infrastructure needed to make the scientific analysis reproducible.

| Category | Stages/files | Why |
| --- | --- | --- |
| Boilerplate / infrastructure | `environment.yml`, `.gitignore`, `src/utils.py`, CLI argument parsing, directory creation, Snakemake shell wrappers | These make the project runnable, reviewable, and reproducible but are not themselves the biological analysis. |
| Reproducibility-critical engineering | `workflows/Snakefile`, `config/*.yaml`, `config/*.yml`, `src/ingest.py`, `docs/real_data.md` | These define the analysis contract, input/output paths, run profiles, and raw/generated data policy. |
| Data-integrity critical | `src/metadata_integrate.py`, `src/qc.py`, `src/build_donor_metadata.py` | These determine whether donor IDs, age, batch, and count data are valid enough for downstream biological analysis. |
| Standard single-cell analysis | `src/doublets_solo.py`, `src/scvi_train.py`, `src/cluster.py`, `src/annotate_celltypist.py`, `src/report.py` | These are common single-cell workflow components, but choices here affect all downstream results. |
| Scientifically central | `src/composition_age.py`, `src/signature_age.py`, `src/age_prediction.py`, `src/sensitivity_age.py`, `src/supplementary_age_plots.py` | These directly address age-associated immune composition, age-linked signatures, predictive age signal, robustness, and interpretation. |
| Scientifically central extension | `src/pseudobulk_aggregate.py`, `src/pseudobulk_edger.R`, `src/validate_pseudobulk_de.py` | Donor-aware cell-type-specific differential expression is implemented and technically accepted; gene-level biological review remains pending. |

## 7. Stages Most Relevant For Employability

The assessment below is inferred from common bioinformatics and computational
biology engineering expectations and the repository's documented
reproducibility and review goals. The most employability-relevant parts are:

- Snakemake orchestration in `workflows/Snakefile`: shows workflow thinking,
  dependency structure, and reproducible execution.
- Config-driven design in `config/`: shows separation of code from run
  parameters and support for demo, pilot, and scaled profiles.
- AnnData and large `.h5ad` handling in `src/ingest.py`, `src/qc.py`, and the
  age-analysis modules: shows practical single-cell data engineering and memory
  caution.
- Metadata harmonization in `src/metadata_integrate.py`: shows awareness that
  donor age, batch, sex, sample, and cell IDs are core scientific assumptions.
- scVI and CellTypist integration: shows use of modern single-cell modeling and
  annotation tools.
- Donor-level composition/signature analyses: shows the right instinct to avoid
  naive per-cell inference for donor-level covariates.
- Age prediction with donor-grouped CV in `src/age_prediction.py`: shows
  applied machine learning with leakage-aware evaluation.
- Sensitivity analysis and supplementary plots: shows robustness and
  communication beyond a single best-looking result.
- README and `docs/assets/`: show the ability to document analysis code and
  results clearly and reproducibly.

## 8. Parts To Understand For Interviews

Inferred from the workflow, code, and documented project goals: you should be
able to explain:

- What problem the workflow solves: PBMC single-cell immune-aging analysis from
  input `.h5ad` to age-associated tables and figures.
- Why donor-level inference matters: age is a donor-level covariate, so treating
  individual cells as independent biological replicates would overstate
  evidence.
- How metadata enters the pipeline: obs-name parsing, external metadata merge,
  canonical donor/batch fields, and conflict/source tracking.
- What scVI contributes: a latent representation used for clustering,
  visualization, and age-prediction features.
- What CellTypist contributes: automated immune-cell label transfer, with
  confidence scores when available.
- What the composition analysis measures: donor-level cell-type fractions
  associated with age, with optional sex/batch adjustment, FDR correction, and
  bootstrap intervals.
- What the signature analysis measures: mean expression of predefined gene
  signatures by donor and cell type, not genome-wide differential expression.
- What the age-prediction analysis evaluates: whether age signal is recoverable
  from replicate-cell-type scVI latent embeddings under grouped regressor
  cross-validation. Because scVI was fitted once on the full cohort, this does
  not establish end-to-end generalization to unseen donors.
- What the sensitivity analysis does: tests whether selected results are stable
  across reasonable analysis-setting changes.
- What remains incomplete: biological interpretation and external validation
  of pseudobulk results, plus more productized real-data acquisition and
  standardization.
- How the project handles large files: raw data and generated results are
  ignored, configs point to local data, and selected figures are copied to
  `docs/assets/`.

## 9. Parts That Can Remain Delegated To AI

This delegation split is inferred from the project structure and scientific
guardrails in `AGENTS.md`.

Appropriate AI-delegated work:

- Drafting or revising README/doc sections from already verified results.
- Generating checklists, PR summaries, and plain-language explanations.
- Creating boilerplate Snakemake rules after the intended inputs/outputs are
  specified by a human.
- Refactoring repetitive plotting or table-writing code without changing
  scientific meaning.
- Proposing config variants for smoke, pilot, and scaled runs.
- Running mechanical consistency checks, such as verifying that Snakefile
  outputs match README output lists.
- Summarizing result tables, as long as the source table and run profile are
  named and overclaiming is avoided.

Work that should not be fully delegated without human review:

- Choosing biological interpretation or claiming immune-aging mechanisms.
- Deciding whether a metadata source is valid and appropriate.
- Choosing statistical thresholds and covariates for final scientific claims.
- Deciding whether to treat a cell type as reliable when donor/cell counts are
  sparse.
- Approving gene-level pseudobulk interpretation and public promotion.
- Approving final public results narrative.
- Making credentialed downloads or externally impactful data operations.

## 10. Open Risks, Missing Documentation, Or Unclear Assumptions

- Donor-aware, cell-type-specific pseudobulk differential expression is
  implemented and technically accepted. Its gene-level biological review and
  any public promotion remain incomplete.
- Real-data acquisition and standardization are not fully productized. The repo
  documents optional Synapse metadata discovery, but the end-to-end raw-data
  acquisition path remains partly manual.
- Metadata provenance needs more formal capture. The configs name local
  metadata files, but a durable source/version/retrieval-date record would make
  the real-data path stronger.
- The full-profile config is named "full" but currently uses `max_cells:
  1000000`; this is a controlled large subset rather than necessarily every
  available cell. The README states this, but the naming could still confuse a
  reviewer.
- Several large-data stages still materialize AnnData in memory. The repo has
  sparse/backed safeguards and pilot/full profiles, but scVI training,
  clustering, metadata integration, and CellTypist annotation can remain
  memory-intensive.
- The active `immune-aging-scvi` environment is stale and currently lacks
  scVI, CellTypist, and XGBoost even though the tracked environment declares
  them. Rebuild or update it before running those stages.
- A focused unittest suite now covers biological-replicate handling,
  pseudobulk contracts, result review, metadata joins, dependency failures,
  covariate policy, and provenance boundaries. Heavy scVI/CellTypist execution
  and the optional Docker edgeR integration test remain separate checks.
- The working tree has an empty `.git` directory. `.gitignore` expresses the
  intended raw/results policy, but history, branch, remote, and tracked state
  cannot be verified until Git metadata is recovered.
- `docs/ai_skills/` contains generated local context dumps and is ignored by
  `.gitignore`. It should not be treated as authoritative project
  documentation.
- The README example results are from one configured GSE164378 run. They are
  descriptive associations and model-performance summaries, not external
  validation.
- Phase 2 validation on the integrated `npj Aging` PBMC atlas is documented as
  future work, not implemented here.
