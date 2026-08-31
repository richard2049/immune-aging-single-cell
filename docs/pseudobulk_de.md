# Longitudinal Pseudobulk Workflow

The maintained pseudobulk path is part of `workflows/Snakefile.longitudinal`
and uses `config/blood_age_atlas_longitudinal.yaml`.

## Stages

1. `src.pseudobulk_aggregate` validates integer sparse counts and writes a
   profiles-by-genes Matrix Market matrix plus profile and gene metadata.
2. `src.run_pseudobulk_dream` invokes the pinned R/Bioconductor runtime.
3. `src.pseudobulk_dream.R` fits the repeated-measures primary model and the
   one-sample-per-subject sensitivity.
4. `src.validate_pseudobulk_de` checks count conservation, metadata, designs,
   model scopes, FDR calculations, versions, logs, and output containment.

## Outputs

Under `results/blood_age_atlas_longitudinal/pseudobulk_de/`:

- `pseudobulk_counts.mtx.gz`
- `pseudobulk_profiles.csv`
- `pseudobulk_genes.csv`
- `celltype_eligibility.csv`
- `aggregation_audit.json`
- `gene_level_results.csv.gz`
- `one_sample_per_subject_results.csv.gz`
- `celltype_result_manifest.csv`
- `design_diagnostics.csv`
- `dream_diagnostic_plots.pdf`
- `runtime_versions.csv`
- `session_info.txt`
- `dream.log`
- `technical_validation.json`

## Targeted Commands

Build the pinned container once and inspect the build log:

```powershell
docker --context desktop-linux build `
  -t immune-aging-dream:bioc-3.23-dream-1.42.0 `
  containers/dream
```

Run only pseudobulk technical validation and its dependencies:

```powershell
python -m snakemake -s workflows/Snakefile.longitudinal -c 1 `
  --configfile config/blood_age_atlas_longitudinal.yaml `
  results/blood_age_atlas_longitudinal/pseudobulk_de/technical_validation.json `
  --rerun-incomplete --printshellcmds
```

Revalidate completed outputs without recomputation:

```powershell
python -m src.validate_pseudobulk_de `
  --config config/blood_age_atlas_longitudinal.yaml `
  --outdir results/blood_age_atlas_longitudinal/pseudobulk_de `
  --report results/blood_age_atlas_longitudinal/pseudobulk_de/technical_validation.json
```

These commands establish technical reviewability only. Gene- and pathway-level
interpretation requires a later D-stage decision and human review.
