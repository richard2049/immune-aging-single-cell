# Real Data Guide

This project expects a pre-standardized `.h5ad` for real runs.

## Minimum AnnData Requirements
- Matrix in `adata.X` with raw counts preferred.
- Unique identifiers:
  - `adata.obs_names` unique cell barcodes
  - `adata.var_names` unique gene symbols/IDs
- Recommended columns in `adata.obs`:
  - `batch` (strongly recommended for scVI integration)
  - `sample_id`
  - `donor_id`
  - `age` (required for age-stratified analyses)
  - `sex` (optional)
  - `condition` (optional)

## Configure Real Mode
Edit `config/config.real.yml`:
- `run.dataset: custom_h5ad`
- `paths.input_h5ad: data/raw/input.h5ad`
- `scvi.categorical_covariates`: include available batch-like columns (for example: `["batch", "donor_id"]`)

Run:
```bash
snakemake -s workflows/Snakefile -c 1 --configfile config/config.real.yml
```

## Standardization Checklist
Before running, verify:
- Counts are non-negative integers (or close to integer UMI counts).
- Mitochondrial genes use `MT-` prefix (or adjust QC code if your naming differs).
- Data remains sparse when possible.

## Memory-Safe Defaults
- Keep `adata.X` sparse; avoid `.toarray()` on large matrices.
- Keep raw counts in `adata.layers["counts"]` (sparse).
- Do not store dense full-matrix normalized layers unless required.
- Use CPU defaults first; switch to GPU only after environment validation.

## Optional Ingestion Extension (Next Step)
A future `src/ingest_real.py` can automate:
1. download from GEO/HCA/source bucket,
2. convert to standardized `.h5ad`,
3. validate required `obs`/`var` fields,
4. write provenance metadata (source URL, accession, date).
