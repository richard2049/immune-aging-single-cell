# Demo Workflow End-to-End Qualification

## Scope

This record covers a clean execution of the public PBMC3K demo from ingestion
through the default `all` target. It qualifies the current software path and
its output contract; it does not validate immune-aging associations.

PBMC3K represents one donor and does not provide donor-level age metadata.
Consequently, the age-association modules are expected to emit empty tables and
clearly labelled placeholder figures. Their successful execution tests graceful
handling of unavailable metadata, not biological inference.

## Execution

- Date: 2026-08-04
- Git commit: `f045a8a424186017afb4e9e0525c839a7668e3c0`
- Platform: Windows, CPU execution
- Python: 3.10.20
- Snakemake: 7.32.4
- scanpy: 1.11.5
- scvi-tools: 1.3.3
- CellTypist: 1.7.1
- Output root: `results/demo_e2e_validation/`
- Duration: approximately 39 minutes

The run used an untracked copy of `config/demo.yaml`. Only `cfg_path` and
`project.out_dir` were changed so that the run used its own configuration path
and a fresh output directory. The effective command was:

```powershell
python -u -m snakemake `
  -s workflows/Snakefile `
  -c 1 `
  all `
  --configfile <local-demo-config> `
  --printshellcmds `
  --rerun-incomplete `
  --latency-wait 30
```

## Input Provenance

| Resource | SHA-256 |
| --- | --- |
| scanpy PBMC3K H5AD | `89a96f1beaa2dd83a687666d3f19a4513ac27a2a2d12581fcd77afed7ea653a1` |
| CellTypist `models.json` | `f617fcd7863b08e6ba2dca1db959d60d6f82f06a0833897e8968c9164a18b132` |
| CellTypist `Immune_All_Low.pkl` | `290874d35dac039d4c9218c343fde4aac1077709b72a331ce7266f6828c36502` |

The CellTypist provenance embedded in the final AnnData identifies model
version `v2`, dated 2022-07-16, and CellTypist 1.7.1.

## Observed Run Inventory

| Checkpoint | Cells | Genes |
| --- | ---: | ---: |
| Raw and metadata-integrated | 2,700 | 32,738 |
| After QC | 2,556 | 13,714 |
| After doublet removal | 2,421 | 13,714 |
| Final annotated AnnData | 2,421 | 13,714 |

The final object contains `X_scVI` and `X_umap`, twelve CellTypist populations,
and model provenance. The report completion marker contains
`report_complete`. All twelve expected PNG files were non-empty; the annotated
UMAP and QC distributions were also inspected visually.

## Technical Acceptance

- All 13 Snakemake jobs completed.
- The workflow exited with status 0.
- No traceback or failed rule was recorded.
- A post-run dry-run reported: `Nothing to be done (all requested files are
  present and up to date).`
- The Git worktree remained clean apart from ignored `.snakemake/`, `data/`,
  `results/`, and Python cache directories.

Four non-fatal dependency warnings were recorded: two PyTorch distribution
support warnings, one scanpy Leiden backend future warning, and one CellTypist
scanpy-version deprecation warning. None interrupted a rule or changed the
expected output contract.

## Interpretation Boundary

This qualification demonstrates that the documented demo can build the
environment-dependent analysis chain and produce its declared outputs on a
clean checkout. It does not test donor-aware age inference because the demo
dataset lacks the required experimental design. Generated inputs, models,
intermediates, figures, tables, and logs remain excluded from version control.
