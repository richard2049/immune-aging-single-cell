# IMMUNE-AGING Project Brief

## Goal
- Build a modular, reproducible single-cell RNA-seq pipeline for PBMC immune-aging analysis.
- Demonstrate production-oriented engineering with `Snakemake`, config-driven execution, and memory-safe handling of large `.h5ad` inputs.
- Answer three practical analysis questions:
  - which immune cell types/states shift with age,
  - which age-linked programs/signatures change by cell type,
  - whether donor age can be predicted from scVI latent embeddings.
- Keep the repository independently reviewable: use a clear structure,
  reproducible commands, interpretable figures and tables, and no committed raw
  data.

## Constraints
- Large real datasets must be handled without dense conversions or unnecessary full-memory loads.
- Raw data, generated `results/`, and `.snakemake/` artifacts stay out of version control.
- The workflow must support both a small demo path and a real-data path.
- Windows/PowerShell use matters; commands and docs should not assume Bash-only behavior.
- Default environment should remain portable; tracked dependencies should reflect the scientific workflow, not local maintainer tooling.

## Key Decisions Already Made
- Core stack: `scanpy` + `scvi-tools` + `CellTypist`.
- Orchestration: `workflows/Snakefile` with YAML configs rather than a monolithic script.
- Real-data metadata integration happens inside the workflow (`metadata` step), not as an external manual side process.
- Real-data execution is split into profiles:
  - `config/config.real.yml` for pilot/subset runs,
  - `config/config.real.full.yml` for scaled/full runs,
  - `config/config.demo.yaml` for smoke tests.
- Phase 1 already includes:
  - preprocessing + annotation,
  - composition-vs-age analysis,
  - signature-vs-age analysis,
  - age prediction,
  - sensitivity analysis,
  - supplementary summary plots.
- `environment.yml` is limited to runtime and analysis dependencies;
  developer-only helper tools are outside the tracked environment contract.

## Unresolved Questions
- The main remaining scientific gap is biological interpretation and external
  validation of the completed donor-aware, cell-type-specific transcriptional
  change analysis.
- Pseudobulk differential expression by biological replicate x cell type is
  implemented as a targeted corrected-profile stage; naive per-cell DE remains
  excluded from the analysis contract.
- A high-priority feasibility audit found that `donor_id` contains 166 labels
  reused across pools, whereas source `Tube_id` identifies 317 conflict-free
  biological replicates. Remediation and targeted downstream regeneration are
  complete under `results/gse164378_full_replicate_corrected/`. Human review
  retained a conservative subset for public presentation.
- A pre-specified corrected-results review is implemented. Its evidence table
  contains 109 rows and prioritizes 12 associations plus one internal
  prediction result. Human dispositions are recorded for all 13 prioritized
  rows, and the approved conservative subset is reflected in the README and
  curated assets.
- The raw CSR count matrix passed a complete integrity scan, and 13 cell types
  pass candidate pseudobulk support/design criteria. The analysis contract is
  approved: 12 populations are primary and NK cells are exploratory. The
  pinned R 4.6.1/Bioconductor 3.23/edgeR 4.10.1 Docker runtime passed a bounded
  integration test. Full-data execution and automated technical acceptance
  completed on 2026-07-31. The pre-specified robustness audit completed on
  2026-08-02 with 212 completed and three non-estimable models, no model
  failures, and a bounded 249-row human-review queue. Biological disposition
  of the pseudobulk gene-level results remains pending.
- The public-facing Phase 1 story now has a concise GSE164378 example section
  in `README.md` and curated figures in `docs/assets/`.
- The real-data path works, but acquisition/standardization is not yet fully productized end-to-end.
- Phase 2 validation on the integrated `npj Aging` PBMC atlas remains future work after Phase 1 is frozen.

## Important Commands
- Demo run:
  ```powershell
  python -m snakemake -s workflows/Snakefile -c 1 --configfile config/config.demo.yaml
  ```
- Pilot real-data run:
  ```powershell
  python -m snakemake -s workflows/Snakefile -c 1 --configfile config/config.real.yml
  ```
- Full/scaled real-data run:
  ```powershell
  python -m snakemake -s workflows/Snakefile -c 1 --configfile config/config.real.full.yml
  ```
- Targeted age-prediction run:
  ```powershell
  python -m snakemake -s workflows/Snakefile -c 1 age_prediction --configfile config/config.real.yml
  ```
- Syntax sanity check:
  ```powershell
  python -m compileall src workflows
  ```
- Pseudobulk DE technical validation without recomputation:
  ```powershell
  python -u -m src.validate_pseudobulk_de `
    --config config/config.real.full.replicate_corrected.yml `
    --outdir results/gse164378_full_replicate_corrected/pseudobulk_de `
    --report results/gse164378_full_replicate_corrected/pseudobulk_de/technical_validation.json
  ```
- Pseudobulk robustness evidence preparation from completed sensitivity models:
  ```powershell
  python -m src.prepare_pseudobulk_evidence `
    --config config/config.real.full.replicate_corrected.yml `
    --pseudobulk-dir results/gse164378_full_replicate_corrected/pseudobulk_de `
    --report-out docs/pseudobulk_robustness_audit.md
  ```
- Targeted corrected-results review:
  ```powershell
  python -m src.review_corrected_results --config config/config.real.full.replicate_corrected.yml --corrected-out results/gse164378_full_replicate_corrected --evidence-out results/gse164378_full_replicate_corrected/review/evidence_table.csv --summary-out results/gse164378_full_replicate_corrected/review/review_summary.json --promotion-out results/gse164378_full_replicate_corrected/review/promotion_manifest.json --report-out docs/corrected_results_review.md
  ```

## Important Files
- `workflows/Snakefile`: current workflow contract and target graph.
- `README.md`: public-facing project summary and run instructions.
- `docs/real_data.md`: real-data schema, ingestion notes, and operational guidance.
- `config/config.demo.yaml`: smoke-test profile.
- `config/config.real.yml`: pilot real-data profile.
- `config/config.real.full.yml`: scaled/full real-data profile.
- `src/metadata_integrate.py`: metadata merge and obs parsing.
- `src/composition_age.py`: donor-level composition analysis.
- `src/signature_age.py`: donor-cell type signature analysis.
- `src/age_prediction.py`: scVI latent-space age prediction and model comparison.
- `src/sensitivity_age.py`: sensitivity runs over analysis settings.
- `src/supplementary_age_plots.py`: supplementary effect/stability figures.
- `src/validate_pseudobulk_de.py`: technical acceptance checks for the full
  donor-aware pseudobulk output contract; it does not interpret biology.
- `src/pseudobulk_robustness.R`: approved edgeR sensitivity scenarios over
  frozen pseudobulk counts.
- `src/run_pseudobulk_robustness.py`: pinned-runtime wrapper for the robustness
  models.
- `src/prepare_pseudobulk_evidence.py`: validates scenario outputs and creates
  bounded evidence tables without biological dispositions.
- `docs/pseudobulk_robustness_audit.md`: tracked technical audit report and
  interpretation boundary.
- `docs/validation/pseudobulk_de_full_run_acceptance.md`: execution inventory,
  validation outcome, and remaining interpretation gates for the full run.
- `src/review_corrected_results.py`: traceable evidence screening and blocked
  promotion manifest for corrected donor-level results.
- `docs/corrected_results_review.md`: human-review report; it is not an
  approved public-results narrative.
- `docs/reviews/corrected_results_dispositions.yml`: tracked human
  dispositions and presentation tiers for prioritized corrected results.
- `src/draft_corrected_public_assets.py`: generates non-public,
  disposition-aware figure and wording drafts for final review.
- `docs/reviews/corrected_public_promotion_record.md`: approved README wording,
  captions, assets, and final sign-off record.

## Important Results
- Base reporting:
  - `results/.../figures/umap_leiden.png`
  - `results/.../figures/umap_cell_type.png`
  - `results/.../tables/cell_type_counts.csv`
- Composition analysis:
  - `results/.../figures/age_celltype_top_trends.png`
  - `results/.../tables/age_celltype_trend_stats.csv`
- Signature analysis:
  - `results/.../figures/signature_age_heatmap.png`
  - `results/.../tables/signature_age_associations.csv`
- Age prediction:
  - `results/.../figures/age_pred_observed_vs_predicted.png`
  - `results/.../tables/age_pred_model_comparison_summary.csv`
- Sensitivity / supplementary:
  - `results/.../sensitivity_age/sensitivity_summary.csv`
  - `results/.../figures/supp_age_effect_ci_forest.png`
  - `results/.../figures/supp_age_sensitivity_stability.png`
- Curated README assets:
  - `docs/assets/gse164378_corrected_composition_core_trends.png`
  - `docs/assets/gse164378_corrected_composition_effect_forest.png`
  - `docs/assets/gse164378_corrected_age_prediction_internal_cv.png`
- Corrected-results review:
  - `results/gse164378_full_replicate_corrected/review/evidence_table.csv`
  - `results/gse164378_full_replicate_corrected/review/review_summary.json`
  - `results/gse164378_full_replicate_corrected/review/promotion_manifest.json`
  - `docs/corrected_results_review.md`
  - `results/gse164378_full_replicate_corrected/review/figures/`
  - `results/gse164378_full_replicate_corrected/review/draft_asset_manifest.json`
- Donor-aware pseudobulk differential expression:
  - `results/gse164378_full_replicate_corrected/pseudobulk_de/gene_level_results.csv.gz`
  - `results/gse164378_full_replicate_corrected/pseudobulk_de/celltype_result_manifest.csv`
  - `results/gse164378_full_replicate_corrected/pseudobulk_de/design_diagnostics.csv`
  - `results/gse164378_full_replicate_corrected/pseudobulk_de/technical_validation.json`
- Pseudobulk robustness and evidence preparation:
  - `results/gse164378_full_replicate_corrected/pseudobulk_de/robustness/candidate_evidence.csv.gz`
  - `results/gse164378_full_replicate_corrected/pseudobulk_de/robustness/review_queue.csv`
  - `results/gse164378_full_replicate_corrected/pseudobulk_de/robustness/celltype_robustness_summary.csv`
  - `results/gse164378_full_replicate_corrected/pseudobulk_de/robustness/audit_summary.json`
  - `docs/pseudobulk_robustness_audit.md`

## Current State
- Phase 1 workflow is implemented and consolidated.
- The branch state was previously cleaned and synchronized around the Phase 1 consolidation point.
- The pseudobulk feasibility stage is implemented in
  `src/pseudobulk_de_feasibility.py`, with generated diagnostics under
  `results/gse164378_full/pseudobulk_de_feasibility/`.
- The approved pseudobulk implementation is in
  `src/pseudobulk_aggregate.py`, `src/pseudobulk_edger.R`, and
  `src/run_edger.py`, with targeted rules in
  `workflows/Snakefile.replicate_corrected`. Its pinned Docker runtime and
  bounded integration test are validated. The full run produced 2,944 profiles
  from 317 biological replicates and 107,948 gene-by-cell-type tests across 13
  populations. `src/validate_pseudobulk_de.py` passed the technical output
  contract; no gene-level interpretation has been approved.
- The decision-card-005 robustness stage is implemented in
  `src/pseudobulk_robustness.R`, `src/run_pseudobulk_robustness.py`, and
  `src/prepare_pseudobulk_evidence.py`. The full audit preserved all 7,970
  global-FDR candidate rows, separated nine exploratory NK-cell rows, and
  created a 249-row queue for later D-class human review without assigning
  automated biological dispositions.
- Phase 1 donor-level results were regenerated with source `Tube_id` exposed as
  `biological_replicate_id`, using the existing annotated H5AD read-only.
- Technical validation passed for 1,000,000 mapped cells, 317 biological
  replicates, zero metadata conflicts, zero grouped-CV leakage, five
  sensitivity scenarios, and eight non-empty corrected figures.
- A descriptive before/after comparison exists under
  `results/gse164378_full_replicate_corrected/comparison/`; it must not be
  treated as biological interpretation.
- Decision card
  `docs/decisions/002_biological_replicate_remediation.md` is implemented.
- The D-class human review documented in
  `docs/decisions/003_corrected_results_interpretation_and_promotion.md`.
  Its evidence assessment, candidate dispositions, professional wording, and
  three disposition-aware figures are implemented, validated, approved, and
  reflected in the README. The four provisional curated assets were removed.
