# Configuration Profiles

The files under `config/` have distinct execution roles. Tracked study
profiles record reproducible contracts; they are not a menu of settings to
edit in place.

| Profile | Role | Workflow | Output | Edit policy |
| --- | --- | --- | --- | --- |
| `config/demo.yaml` | Small fixture run and CI DAG check | `workflows/Snakefile` | `results/` | Fixed demo contract |
| `config/gse164378_pilot.yaml` | Reproducible 50,000-cell GSE164378 smoke run | `workflows/Snakefile` | `results/gse164378_pilot/` | Fixed study profile |
| `config/gse164378_1m.yaml` | Accepted 1,000,000-cell GSE164378 run | `workflows/Snakefile` | `results/gse164378_full/` (legacy path) | Fixed study profile |
| `config/gse164378_corrected.yaml` | Replicate-aware downstream analyses and pseudobulk evidence | `workflows/Snakefile.replicate_corrected` | `results/gse164378_full_replicate_corrected/` | Fixed study profile |
| `config/custom.example.yaml` | Starting point for another standardized H5AD | `workflows/Snakefile` | `results/custom/` | Copy before editing |

## Custom Data

Create an ignored local profile rather than changing a tracked study profile:

```powershell
Copy-Item config/custom.example.yaml config/custom.local.yaml
```

In the copy, set `cfg_path: config/custom.local.yaml` and choose a distinct
`project.out_dir`. Review every metadata field and analysis threshold against
the new study design. In particular, `biological_replicate_id` must identify
the independent experimental unit from authoritative metadata; a donor-like
label or parsed filename is not sufficient evidence.

The example assumes required fields are already present in `adata.obs` because
`metadata.enabled` is `false`. Enable the metadata integration step only after
defining and validating the join key, source columns, cardinality, and unmatched
row policy described in [the real-data guide](real_data.md).

## Inheritance

`config/gse164378_corrected.yaml` extends `gse164378_1m.yaml`. The loader in
`src/utils.py` resolves the base path relative to the child profile and deep
merges mapping values. The corrected profile changes only the checkpoint,
output, biological-replicate contract, and downstream settings required for
the accepted replicate-aware analyses.

Each profile repeats its own path in `cfg_path` because the current workflow
passes that file to Python modules. Automated tests require this value to match
the tracked filename and inspect the fully resolved configuration, including
inherited values.

## Choosing A Profile

- Use the demo profile for installation checks and CI-like validation.
- Use the pilot profile to test the GSE164378 path with bounded compute; do not
  treat its outputs as the maintained one-million-cell results.
- Use the one-million-cell profile to rebuild the accepted provisional
  GSE164378 checkpoint.
- Use the corrected profile only with the validated full annotated checkpoint.
- Use a local copy of the custom example for any other dataset.

Changing input inclusion, replicate identity, covariates, model formulas, or
inferential thresholds is an analysis-design change and requires scientific
review. A new output directory should be used whenever such a change could
make results non-comparable.

The current source H5AD contains 1,916,367 cells, whereas the accepted
checkpoint contains 1,000,000. The historical
`results/gse164378_full/` directory name is retained to avoid moving or
invalidating that large checkpoint; it does not mean that all source cells
were included. A future uncapped run requires a separate profile, output
directory, compute qualification, and scientific validation.
