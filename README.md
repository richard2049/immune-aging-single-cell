# Immune Aging Single-Cell

[![CI](https://github.com/richard2049/immune-aging-single-cell/actions/workflows/ci.yml/badge.svg)](https://github.com/richard2049/immune-aging-single-cell/actions/workflows/ci.yml)

Modular Snakemake pipeline for PBMC single-cell RNA-seq analysis focused on immune aging.

Stack: `scanpy` + `scvi-tools` + `CellTypist`.

## Key Properties
- Reproducible, config-driven workflow design (`Snakemake` + YAML configs).
- Memory-safe single-cell processing (sparse counts, no dense graph conversion).
- Standard single-cell analysis stages: QC, doublet removal, scVI latent
  modeling, clustering, and automated annotation.
- Structured, inspectable outputs: `.h5ad` intermediates, figures, and summary
  tables.

## Documentation
- [Documentation index](docs/README.md)
- [Analysis decisions](docs/analysis_decisions.md)
- [Real-data guide](docs/real_data.md)
- [Scientific verification contract](docs/verification_contract.md)

## Pipeline
1. `ingest` -> `results/01_raw.h5ad`
2. `metadata` -> `results/01_meta.h5ad` (obs ID parsing + optional external metadata merge)
3. `qc` -> `results/02_qc.h5ad`
4. `doublets` (SOLO) -> `results/03_nodoublets.h5ad`
5. `scvi_train` -> `results/04_scvi.h5ad` + `results/models/scvi_model/`
6. `cluster` -> `results/05_clustered.h5ad`
7. `annotate` (CellTypist) -> `results/06_annotated.h5ad`
8. `report` -> figures + tables + `results/reports.done`
9. `composition_age` -> age-stratified composition figures + trend tables
10. `signature_age` -> sample-level signature trends with subject-clustered inference
11. `age_prediction` -> chronological age prediction with subject-grouped cross-validation
12. `sensitivity_age` -> parameter-sensitivity runs for age analyses
13. `supplementary_age` -> supplementary summary figures for effect sizes and stability

## Quickstart (Demo)
```bash
conda env create -f environment.yml
conda activate immune-aging-scvi

snakemake -s workflows/Snakefile -c 1 --configfile config/demo.yaml
```

For an existing environment, apply the tracked dependency contract explicitly:

```powershell
$mamba = Join-Path ((conda info --base).Trim()) "Library\bin\mamba.exe"
if (-not (Test-Path $mamba)) { throw "mamba.exe not found at $mamba" }
& $mamba env update -n immune-aging-scvi -f environment.yml --prune
if ($LASTEXITCODE -ne 0) { throw "Environment update failed; do not run validation in the stale environment." }
conda activate immune-aging-scvi
python -c "import scanpy, scvi, celltypist, xgboost; print(scanpy.__version__, scvi.__version__, celltypist.__version__, xgboost.__version__)"
if ($LASTEXITCODE -ne 0) { throw "Required scientific packages are unavailable." }
```

The update can require a full dependency solve. Do not run it concurrently with
Snakemake or another conda operation. If channel downloads fail, follow the
network checks in `docs/troubleshooting.md`; do not continue with workflow
validation in the stale environment.

### Optional NVIDIA GPU Environment

`environment.yml` is the portable CPU contract used by CI. On a supported
NVIDIA system, create the separately qualified CUDA 12.1 environment:

```powershell
conda env create -n immune-aging-scvi-gpu -f environment-gpu.yml
conda activate immune-aging-scvi-gpu
python -c "import torch, scvi; print(torch.__version__, torch.version.cuda, torch.cuda.is_available())"
```

The CUDA-enabled PyTorch build can execute either backend. Select `cpu` or
`gpu` through `scvi.accelerator` in the local config used for the run; changing
the accelerator does not replace the accepted seed, covariates, or input
contract. Qualify both modes with a bounded smoke test before a large run and
use a CPU profile when CUDA is unavailable. See
[`docs/troubleshooting.md`](docs/troubleshooting.md#gpu-training) for the
driver and runtime checks.

Windows/PowerShell equivalent:
```powershell
python -m snakemake -s workflows/Snakefile -c 1 --configfile config/demo.yaml
```

Run a specific target:
```powershell
python -m snakemake -s workflows/Snakefile -c 1 age_prediction --configfile config/blood_age_atlas_pilot.yaml
```

## Maintenance Rebuilds

The demo command above is a smoke-test path, not a rebuild of the Blood Age
Atlas analysis. The maintained real-data path has two sequential stages:

1. `workflows/Snakefile` reconstructs and validates the one-million-cell scVI
   checkpoint, creates and validates the clustered checkpoint, and only then
   runs annotation under `results/blood_age_atlas_1m/`.
2. `workflows/Snakefile.longitudinal` applies the explicit subject, sample-unit,
   and technical-library contract and writes downstream outputs under
   `results/blood_age_atlas_longitudinal/`.

Directories beginning with `results/gse164378` are historical checkpoints from
before the study-identity correction. They are not output targets and must not
be overwritten or treated as current evidence.

Run the following from the repository root in an activated
`immune-aging-scvi` environment. Keep `-c 1` for the conservative memory-safe
route, do not run the two Snakemake commands concurrently, and ensure Docker
Desktop is available for the pinned `dream` stage.

```powershell
$numbaCache = Join-Path (Resolve-Path ".tmp").Path "numba-cache"
New-Item -ItemType Directory -Force -Path $numbaCache | Out-Null
$env:NUMBA_CACHE_DIR = $numbaCache

$requiredInputs = @(
    "data/raw/raw_counts_h5ad/pbmc_gex_raw_with_var_obs.h5ad",
    "data/raw/raw_counts_h5ad/all_pbmcs/all_pbmcs_metadata.csv"
)

foreach ($path in $requiredInputs) {
    if (-not (Test-Path -LiteralPath $path)) {
        throw "Required input not found: $path"
    }
}

docker --context desktop-linux image inspect `
    immune-aging-dream:bioc-3.23-dream-1.42.0 | Out-Null
if ($LASTEXITCODE -ne 0) {
    throw "The pinned dream Docker image is unavailable."
}

$logDir = "results/logs"
New-Item -ItemType Directory -Force -Path $logDir | Out-Null
$logPath = Join-Path $logDir (
    "maintenance_rebuild_{0}.log" -f (Get-Date -Format "yyyyMMdd_HHmmss")
)

Start-Transcript -Path $logPath
try {
    python -m snakemake `
        -s workflows/Snakefile `
        -c 1 `
        scvi_train `
        --configfile config/blood_age_atlas_1m.yaml `
        --rerun-incomplete `
        --printshellcmds

    if ($LASTEXITCODE -ne 0) {
        throw "scVI reconstruction failed; do not continue."
    }

    python -m src.validate_scvi_checkpoint `
        --config config/blood_age_atlas_1m.yaml `
        --inp results/blood_age_atlas_1m/03_nodoublets.h5ad `
        --checkpoint results/blood_age_atlas_1m/04_scvi.h5ad `
        --report results/blood_age_atlas_1m/validation/04_scvi_checkpoint.json

    if ($LASTEXITCODE -ne 0) {
        throw "The scVI checkpoint gate failed; do not run clustering."
    }

    python -m snakemake `
        -s workflows/Snakefile `
        -c 1 `
        cluster `
        --configfile config/blood_age_atlas_1m.yaml `
        --rerun-incomplete `
        --printshellcmds

    if ($LASTEXITCODE -ne 0) {
        throw "Clustering failed; do not run annotation."
    }

    python -m src.validate_cluster_checkpoint `
        --config config/blood_age_atlas_1m.yaml `
        --inp results/blood_age_atlas_1m/04_scvi.h5ad `
        --checkpoint results/blood_age_atlas_1m/05_clustered.h5ad `
        --report results/blood_age_atlas_1m/validation/05_cluster_checkpoint.json

    if ($LASTEXITCODE -ne 0) {
        throw "The clustered checkpoint gate failed; do not run annotation."
    }

    python -m snakemake `
        -s workflows/Snakefile `
        -c 1 `
        annotate `
        --configfile config/blood_age_atlas_1m.yaml `
        --rerun-incomplete `
        --printshellcmds

    if ($LASTEXITCODE -ne 0) {
        throw "Annotation failed; do not run the longitudinal stage."
    }

    python -m snakemake `
        -s workflows/Snakefile.longitudinal `
        -c 1 `
        --configfile config/blood_age_atlas_longitudinal.yaml `
        --rerun-incomplete `
        --printshellcmds

    if ($LASTEXITCODE -ne 0) {
        throw "The longitudinal workflow failed."
    }
}
finally {
    Stop-Transcript
}
```

The workspace-scoped Numba cache avoids writes to protected environment or
user-cache directories on Windows. It contains disposable compiled cache
files only and is excluded from version control. See
[`docs/troubleshooting.md`](docs/troubleshooting.md#windows-cache-permissions)
if Scanpy or scVI stalls during import.

The longitudinal `all` target includes unit mapping, scientific checkpoint
audit, outcome-blind annotation and design qualification, subject-aware
composition and signature models, subject-grouped age prediction, pseudobulk
aggregation, `dream` repeated-measures inference, the pre-specified
one-sample-per-subject sensitivity, and technical validation.
Snakemake normally reruns only incomplete or out-of-date jobs. Do not use
`--forceall` for checkpoint requalification; it would exceed the localized
reconstruction boundary.

Before the first inferential run, generate only the review evidence:

```powershell
python -m snakemake `
    -s workflows/Snakefile.longitudinal `
    -c 1 `
    annotation_design_qualification `
    --configfile config/blood_age_atlas_longitudinal.yaml `
    --printshellcmds
```

Follow the human-review procedure in
[`docs/annotation_mapping_review_guide.md`](docs/annotation_mapping_review_guide.md).
Use `docs/reviews/annotation_mapping_review.yml` to record conclusions and the
generated `validation/annotation_mapping_review_table.csv` as the compact
evidence index. The full workflow fails safely until every mapping row has
received human approval and the fingerprint-matched approval sentinel has been
generated.

Run the longitudinal workflow only after the reconstructed scVI and clustered
checkpoint gates pass and `results/blood_age_atlas_1m/06_annotated.h5ad` is
complete. Its first rules run the configured scientific checkpoint audit and
annotation qualification. Inspect all validation reports and do not publish
regenerated results or claims without independent scientific review.

## Validation
Use the project conda environment before running workflow checks. The tracked
environment is CPU-first for portability.

```powershell
python -m compileall src workflows
python -m unittest discover -s tests -v
python -m snakemake -s workflows/Snakefile -c 1 -n --configfile config/demo.yaml
```

GitHub Actions applies the same Ruff and compilation checks, runs the unit-test
suite in the tracked scientific environment, and constructs the demo workflow
DAG. Docker-qualified `dream` tests, real-data execution, and full analysis runs
remain explicit validation steps outside routine CI. A passing CI run confirms
the automated technical checks only; it does not establish biological validity.

Development-only linting is kept outside the scientific runtime environment:

```powershell
$base_python = Join-Path ((conda info --base).Trim()) "python.exe"
& $base_python -m venv .venv-dev
$dev_python = Join-Path (Resolve-Path ".venv-dev") "Scripts\python.exe"
& $dev_python -m pip install -r requirements-dev.txt
& $dev_python -m ruff check src tests workflows
& $dev_python -m ruff format --check src tests workflows
```

The optional GPU profile is not exercised by GitHub Actions because hosted and
local GPU drivers are machine-specific. Its scientific Python versions match
the CPU contract, while the PyTorch backend is qualified separately.

## Blood Age Atlas Requalification Status

The analysed study is the Terekhova et al. Blood Age Atlas distributed through
Synapse as `syn49637038`, not GEO accession `GSE164378`. Source metadata describe
317 sample units from 166 people. The workflow assigns distinct roles to
`subject_id` (`Donor_id`), `sample_unit_id` (`Tube_id`), and
`technical_library_id` (`File_name`).

Historical downstream outputs treated sample units as independent people and
used the wrong study name. They remain available only for provenance and must
not support current associations, prediction metrics, genes, pathways, or
figures. New results require an isolated rebuild, L3 technical validation, and
human scientific review.

Composition and signature models use subject-clustered uncertainty. Prediction
folds are grouped by subject and remain transductive because scVI is trained on
the full cohort. Pseudobulk inference uses a subject random intercept and a
separate deterministic one-sample-per-subject sensitivity. No gene-, pathway-,
or cell-population claim is currently approved.

## Report Outputs
Generated by `src/report.py`:
- `results/figures/umap_leiden.png`
- `results/figures/umap_cell_type.png`
- `results/figures/qc_distributions.png`
- `results/figures/cell_type_composition.png`
- `results/tables/cell_type_counts.csv`
- `results/tables/cell_type_fractions.csv`
- `results/tables/cluster_celltype_crosstab.csv`

Generated by `src/composition_age.py`:
- `results/figures/age_celltype_composition_by_bin.png`
- `results/figures/age_celltype_top_trends.png`
- `results/tables/age_celltype_fraction_by_sample_unit.csv`
- `results/tables/age_celltype_trend_stats.csv`

Generated by `src/signature_age.py`:
- `results/figures/signature_age_heatmap.png`
- `results/figures/signature_age_top_associations.png`
- `results/tables/signature_scores_by_sample_unit_celltype.csv`
- `results/tables/signature_age_associations.csv`
- `results/tables/signature_gene_coverage.csv`

Generated by `src/age_prediction.py`:
- `results/figures/age_pred_observed_vs_predicted.png`
- `results/figures/age_pred_mae_by_celltype.png`
- `results/tables/age_pred_cv_predictions.csv`
- `results/tables/age_pred_metrics.csv`
- `results/tables/age_pred_model_comparison_summary.csv`
  - `age_pred_metrics.csv` compares candidate models under identical
    subject-grouped cross-validation folds and subject-clustered bootstrap
    intervals.
  - `age_pred_model_comparison_summary.csv` adds fold-level winner counts and bootstrap CI for paired MAE deltas vs best Ridge.

Generated by `src/sensitivity_age.py`:
- `results/sensitivity_age/sensitivity_manifest.csv`
- `results/sensitivity_age/sensitivity_summary.csv`

Generated by `src/supplementary_age_plots.py`:
- `results/figures/supp_age_effect_ci_forest.png`
- `results/figures/supp_age_sensitivity_stability.png`

## Real-Data Mode
Copy `config/custom.example.yaml` to an ignored local profile before adapting
it to a new dataset:

```powershell
Copy-Item config/custom.example.yaml config/custom.local.yaml
```

Set `cfg_path: config/custom.local.yaml`, then review the input path, metadata
schema, subject/sample/library fields, covariates, thresholds, and output
directory. The tracked Blood Age Atlas profiles record fixed study runs and are not
generic templates.

Detailed schema and ingestion notes are in `docs/real_data.md`; profile roles
and inheritance are documented in `docs/configuration.md`.

## Notes
- Default configs use CPU for portability on Windows/macOS/Linux and match the
  CPU-only PyTorch environment in `environment.yml`. The optional
  `environment-gpu.yml` can run those same profiles on CPU or a locally
  qualified NVIDIA GPU.
- Large data files and `results/` are gitignored.
- Use `config/blood_age_atlas_pilot.yaml` for the bounded 50,000-cell Blood Age
  Atlas smoke run and `config/blood_age_atlas_1m.yaml` for the one-million-cell
  checkpoint. An uncapped analysis is not currently a tracked run profile.
- Store figures referenced by the README under `docs/assets/`
  (`.png/.jpg/.jpeg/.webp`).
- Keep `.snakemake/`, raw downloads, and generated artifacts out of version control.

## Citation
Citation metadata for release `v0.1.0` is provided in
[`CITATION.cff`](CITATION.cff).

## License
Unless otherwise noted, original code and documentation in this repository are
available under the [BSD 3-Clause License](LICENSE). External datasets,
pretrained models, and third-party software remain subject to their original
licenses and terms.
