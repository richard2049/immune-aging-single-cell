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
10. `signature_age` -> donor-level signature trends vs age (figures + stats)
11. `age_prediction` -> chronological age prediction from donor-celltype scVI embeddings
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
python -m snakemake -s workflows/Snakefile -c 1 age_prediction --configfile config/gse164378_pilot.yaml
```

## Maintenance Rebuilds

The demo command above is a smoke-test path, not a rebuild of the current
GSE164378 analysis. The maintained one-million-cell result path has two
sequential stages:

1. `workflows/Snakefile` creates the one-million-cell annotated checkpoint and provisional
   outputs under `results/gse164378_full/`.
2. `workflows/Snakefile.replicate_corrected` reuses that checkpoint and creates
   the authoritative replicate-corrected outputs under
   `results/gse164378_full_replicate_corrected/`.

Run the following from the repository root in an activated
`immune-aging-scvi` environment. Keep `-c 1` for the conservative memory-safe
route, do not run the two Snakemake commands concurrently, and ensure Docker
Desktop is available for the pinned edgeR stages.

```powershell
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
    immune-aging-edger:bioc-3.23-edger-4.10.1 | Out-Null
if ($LASTEXITCODE -ne 0) {
    throw "The pinned edgeR Docker image is unavailable."
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
        --configfile config/gse164378_1m.yaml `
        --rerun-incomplete `
        --printshellcmds

    if ($LASTEXITCODE -ne 0) {
        throw "The one-million-cell workflow failed; do not run the corrected stage."
    }

    python -m snakemake `
        -s workflows/Snakefile.replicate_corrected `
        -c 1 `
        all pseudobulk_robustness_evidence `
        --configfile config/gse164378_corrected.yaml `
        --rerun-incomplete `
        --printshellcmds

    if ($LASTEXITCODE -ne 0) {
        throw "The replicate-corrected workflow failed."
    }
}
finally {
    Stop-Transcript
}
```

The explicit `pseudobulk_robustness_evidence` target extends the corrected
workflow's default `all` target through pseudobulk aggregation, edgeR modeling,
technical validation, sensitivity models, and bounded evidence preparation.
Snakemake normally reruns only incomplete or out-of-date jobs. Add `--forceall`
to both commands only for a deliberate complete recomputation, including scVI
training and annotation.

Run only the second Snakemake command when
`results/gse164378_full/06_annotated.h5ad` has passed the configured scientific
checkpoint audit. After either route, inspect validation reports and
`git diff -- docs`; do not publish regenerated reports or claims without the
required scientific review.

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
DAG. Docker-qualified edgeR tests, real-data execution, and full analysis runs
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

## GSE164378 Requalification Status

The one-million-cell workflow and donor-aware pseudobulk implementation remain
available, but numerical examples and curated result figures are temporarily
withheld from the main project presentation. A safeguards audit found that the
historical annotated checkpoint predates the current scVI, clustering, and
CellTypist provenance records. The composition workflow also now restores
zero-abundance replicate-cell-type combinations explicitly.

The checkpoint and all affected downstream outputs must therefore be rebuilt,
technically validated, and reviewed before any association, prediction metric,
gene, pathway, or figure is presented as a current result. This status does not
establish that earlier result directions were incorrect; it means that they do
not yet satisfy the current reproducibility and acceptance contract.

The accepted biological replicate remains source `Tube_id`, exposed as
`biological_replicate_id`. Gene-level pseudobulk interpretation remains pending
and no gene- or pathway-level claim is approved.

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
- `results/tables/age_celltype_fraction_by_donor.csv`
- `results/tables/age_celltype_trend_stats.csv`

Generated by `src/signature_age.py`:
- `results/figures/signature_age_heatmap.png`
- `results/figures/signature_age_top_associations.png`
- `results/tables/signature_scores_by_donor_celltype.csv`
- `results/tables/signature_age_associations.csv`
- `results/tables/signature_gene_coverage.csv`

Generated by `src/age_prediction.py`:
- `results/figures/age_pred_observed_vs_predicted.png`
- `results/figures/age_pred_mae_by_celltype.png`
- `results/tables/age_pred_cv_predictions.csv`
- `results/tables/age_pred_metrics.csv`
- `results/tables/age_pred_model_comparison_summary.csv`
  - `age_pred_metrics.csv` compares candidate models under identical
    donor-grouped cross-validation folds.
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
schema, biological-replicate field, covariates, thresholds, and output
directory. The tracked GSE164378 profiles record fixed study runs and are not
generic templates.

Detailed schema and ingestion notes are in `docs/real_data.md`; profile roles
and inheritance are documented in `docs/configuration.md`.

## Notes
- Default configs use CPU for portability on Windows/macOS/Linux and match the
  CPU-only PyTorch environment in `environment.yml`. The optional
  `environment-gpu.yml` can run those same profiles on CPU or a locally
  qualified NVIDIA GPU.
- Large data files and `results/` are gitignored.
- Use `config/gse164378_pilot.yaml` for the bounded 50,000-cell GSE164378
  smoke run and `config/gse164378_1m.yaml` for the accepted one-million-cell
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
