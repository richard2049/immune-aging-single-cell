# Scientific Reproducibility Remediation

## Scope

- Date: 2026-08-01
- Decision: `docs/decisions/004_scientific_reproducibility_remediation.md`
- Stage classification: C - scientifically important
- Verification level: L3 - claim-critical
- Interpretation boundary: no result was regenerated or promoted as a new
  biological claim.

## Acceptance Contract

The implementation was checked against `docs/verification_contract.md`, the
validated biological-replicate contract in decision card 002, the global and
repository `AGENTS.md` files, and the global
`scientific-code-verification` Skill. Existing replicate and pseudobulk tests
were preserved without changing assertions or scientific thresholds.

## Finding Resolution

| Audit finding | Resolution | Verification | Status |
| --- | --- | --- | --- |
| Empty `.git` directory prevents version-control audit | No replacement history was invented. The stale roadmap claims were marked blocked. | `.git` exists with zero entries; `git status` still fails. | Blocked: trusted remote or backup required |
| Metadata join silently dropped duplicate keys | Exact duplicate rows may collapse; conflicting duplicate payloads fail. The merge validates `many_to_one`, row count, order, key presence, and unmatched cells. | Positive, duplicate-conflict, and unmatched-cell unit tests pass. | Resolved |
| Metadata helper grouped reused `donor_id` labels | The helper now requires an explicit authoritative replicate field and groups by `biological_replicate_id`; it rejects replicate-level conflicts and writes outside `data/raw`. | Reused-label, conflict, and raw-output unit tests pass. | Resolved |
| Configured prediction model could disappear when a dependency was missing | Requested XGBoost or LightGBM models now raise an actionable import error. XGBoost is declared in `environment.yml`. | Bounded missing-dependency test passes; environment YAML parses; a prior network-enabled mamba dry-run resolved the specification. | Repository contract resolved |
| Active conda environment lacks scVI, CellTypist, and XGBoost | Direct dependencies are version constrained and the README documents update and import checks. | Module discovery confirms the active environment remains stale. | Manual environment update required |
| Missing scVI covariates were silently skipped | Required by default; an alternative model must set `allow_missing_covariates: true` explicitly. | Strict and explicit-alternative unit tests pass. | Resolved |
| Stochastic stages did not consistently use the configured seed | scVI, SOLO, neighbors, UMAP, and Leiden now consume `run.seed`; provenance is saved in AnnData. Every primary config defines the seed. | Config tests, source inspection, and compilation pass. | Resolved in code; heavy retraining not run |
| CellTypist model provenance was incomplete | Annotation loads one resolved model object and records package version, model metadata, path, SHA-256, and parameters; optional configured checksums are enforced. | Bounded hash/provenance and mismatch tests pass. | Resolved in code; real model run pending environment repair |
| Prediction wording implied a stronger validation scope | README, review, interview guide, retrospective, and verification contract state that grouped regressor CV uses a full-cohort `X_scVI` representation and is transductive. | Targeted text inspection. | Resolved for current claims |
| Global verification Skill was not discoverable | Replaced the malformed frontmatter delimiter with `---`. The Unicode dashes were valid UTF-8, not corrupted file content. | Frontmatter parses as YAML with the expected name and description. | Resolved |
| No repository lint contract | Added a minimal Ruff configuration and a separate development-only requirements file. | Configuration inspection; Ruff was not installed in the active environment. | Configured; manual lint installation required |
| Large legacy functions reduce maintainability | No broad refactor was mixed into this L3 repair because it would alter validated scientific modules without a behavioral need. New validation logic was kept in small functions. | Function-scope review. | Residual low-risk technical debt |

## Commands And Results

- `python -m unittest discover -s tests -v`: 33 tests ran; 32 passed and the
  optional Docker edgeR integration test was skipped because its opt-in
  variable was not set.
- `python -m py_compile ...`: passed for every changed Python module and test.
- YAML parsing: `environment.yml` and all five tracked config profiles passed.
- A prior `mamba env update ... --dry-run` resolved successfully through a
  network-enabled execution route and did not modify the environment. A later
  user-terminal update and diagnostic retries could not reach any Anaconda
  channel over HTTPS. The specification is resolvable, but the installed
  environment remains stale until local network access is restored.
- Existing corrected-result acceptance JSON files still report `passed: true`.
- Snakemake DAG construction was not repeated: it has previously stalled on
  this Windows checkout, and the active environment lacks three required
  runtime packages. No workflow lock or incomplete marker was present.
- Ruff execution was attempted from an isolated temporary installation, but
  package download access to PyPI was blocked. The temporary directory was
  removed; the manual command below remains required.

## Manual Recovery Actions

### Rebuild The Local Environment

Run from the repository root when no workflow or conda process is active:

```powershell
$mamba = Join-Path ((conda info --base).Trim()) "Library\bin\mamba.exe"
if (-not (Test-Path $mamba)) { throw "mamba.exe not found at $mamba" }
& $mamba env update -n immune-aging-scvi -f environment.yml --prune
if ($LASTEXITCODE -ne 0) { throw "Environment update failed; stop here." }
conda activate immune-aging-scvi
python -c "import scanpy, scvi, celltypist, xgboost; print(scanpy.__version__, scvi.__version__, celltypist.__version__, xgboost.__version__)"
if ($LASTEXITCODE -ne 0) { throw "Required scientific packages are unavailable." }
python -m unittest discover -s tests -v
```

Do not run the broad workflow until the import check completes. A future
CellTypist annotation run will store the real model checksum; existing H5AD
checkpoints cannot be assigned a checksum retroactively without the exact model
file used for their creation.

### Recover Git Metadata

Obtain the authoritative remote URL and default branch or a trusted backup.
Preserve the empty directory until recovery is verified:

```powershell
Rename-Item .git .git.empty-backup
git clone <authoritative-remote-url> ..\immune-aging-scvi-git-recovery
git -C ..\immune-aging-scvi-git-recovery status
```

Compare the clone with this working tree before restoring `.git`. If it is the
same repository and intended branch, copy the clone's `.git` directory into
this working tree, run `git status --short`, and inspect every reported change.
Do not use `git init` as if it recovered the missing history.

### Run Lint

```powershell
python -m pip install -r requirements-dev.txt
python -m ruff check src tests workflows
python -m ruff format --check src tests workflows
```

Lint installation is development tooling and should not be added to the
scientific runtime environment.
