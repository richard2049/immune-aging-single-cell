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
- `metadata.enabled: true` to run `src/metadata_integrate.py`
- `metadata.parse_obs_names: true` to derive `sample_id`, `donor_id`, `batch` from `obs_names`
- `metadata.table_path: ...` (optional) to merge donor-level metadata such as `age`, `sex`, `condition`
- `scvi.categorical_covariates`: include available batch-like columns (for example: `["batch", "donor_id"]`)

Run:
```bash
snakemake -s workflows/Snakefile -c 1 --configfile config/config.real.yml
```

PowerShell:
```powershell
python -m snakemake -s workflows/Snakefile -c 1 --configfile config/config.real.yml
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

## Build Donor Metadata
Use `src/build_donor_metadata.py` to create a donor-level table used by `metadata.table_path`.

Base run (derive donor/sample from `obs_names` only):
```bash
python -m src.build_donor_metadata \
  --h5ad data/raw/raw_counts_h5ad/pbmc_gex_raw_with_var_obs.h5ad \
  --out data/raw/raw_counts_h5ad/donor_metadata.csv
```

Run with supplementary clinical table (age/sex/cohort):
```bash
python -m src.build_donor_metadata \
  --h5ad data/raw/raw_counts_h5ad/pbmc_gex_raw_with_var_obs.h5ad \
  --supp data/raw/raw_counts_h5ad/clinical_metadata.csv \
  --out data/raw/raw_counts_h5ad/donor_metadata.csv
```

Then in `config/config.real.yml`:
- `metadata.table_path: data/raw/raw_counts_h5ad/donor_metadata.csv`
- `metadata.table_columns: ["age", "sex", "cohort"]`

## Data Source Notes (Important)
- Your current `raw_counts_h5ad.tar.gz` content is donor-rich but does not include age directly in `obs`.
- For age analyses, you need an additional donor-level clinical file (age/sex/etc.).
- The Immunity 2023 healthy blood atlas is distributed via Synapse/AWS resources; download the companion donor metadata/clinical file and pass it as `--supp`.

### Getting Supplementary Donor Metadata
1. Login and inspect Synapse folder content (replace token and IDs as needed):
```python
import synapseclient
syn = synapseclient.Synapse()
syn.login(authToken="YOUR_PAT")
for ch in syn.getChildren("syn56693935"):
    print(ch["id"], ch["name"])
```

2. Download only metadata-like files by ID (avoid downloading everything):
```python
meta_ids = ["synXXXXX", "synYYYYY"]  # files that contain donor age/sex/cohort
for sid in meta_ids:
    syn.get(sid, downloadLocation="data/raw/raw_counts_h5ad")
```

3. Build donor metadata table:
```bash
python -m src.build_donor_metadata \
  --h5ad data/raw/raw_counts_h5ad/pbmc_gex_raw_with_var_obs.h5ad \
  --supp data/raw/raw_counts_h5ad/<downloaded_metadata_file>.csv \
  --out data/raw/raw_counts_h5ad/donor_metadata.csv
```
