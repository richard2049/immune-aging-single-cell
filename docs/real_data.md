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
  - `biological_replicate_id` from authoritative source metadata
  - `donor_id` as a source label when distinct from the replicate identifier
  - `age` (required for age-stratified analyses)
  - `sex` (optional)
  - `condition` (optional)

## Configure Real Mode
Do not edit a tracked GSE164378 study profile. Copy the generic example to an
ignored local profile:

```powershell
Copy-Item config/custom.example.yaml config/custom.local.yaml
```

In `config/custom.local.yaml`, first set:
- `cfg_path: config/custom.local.yaml`
- `run.dataset: custom_h5ad`
- `paths.input_h5ad: data/raw/input.h5ad`
- a distinct `project.out_dir`

Then review:
- `metadata.enabled: true` to run `src/metadata_integrate.py`
- `metadata.parse_obs_names: true` only when its tokenization has been validated;
  parsed labels do not establish the biological replicate
- `metadata.table_path: ...` (optional) to merge donor-level metadata such as `age`, `sex`, `condition`
- `scvi.categorical_covariates`: include available batch-like columns (for example: `["batch", "donor_id"]`)

The profile roles and inheritance contract are described in
[`docs/configuration.md`](configuration.md).

Run:
```bash
snakemake -s workflows/Snakefile -c 1 --configfile config/custom.local.yaml
```

PowerShell:
```powershell
python -m snakemake -s workflows/Snakefile -c 1 --configfile config/custom.local.yaml
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
- Use CPU defaults first; switch to GPU only after installing and validating a
  CUDA-enabled PyTorch build.

## Optional Ingestion Extension (Next Step)
A future `src/ingest_real.py` can automate:
1. download from GEO/HCA/source bucket,
2. convert to standardized `.h5ad`,
3. validate required `obs`/`var` fields,
4. write provenance metadata (source URL, accession, date).

## Build Replicate Metadata

`src/build_donor_metadata.py` creates a derived replicate-level audit table. It
does not infer the independent unit from a donor-like token. Supply either an
authoritative AnnData observation field with `--replicate-obs-col`, or an
explicit cell-level metadata join with `--cell-col` and `--replicate-col`.

For GSE164378, `Tube_id` is the validated biological-replicate field:

```powershell
python -m src.build_donor_metadata `
  --h5ad data/raw/raw_counts_h5ad/pbmc_gex_raw_with_var_obs.h5ad `
  --supp data/raw/raw_counts_h5ad/all_pbmcs/all_pbmcs_metadata.csv `
  --cell-col "Unnamed: 0" `
  --replicate-col Tube_id `
  --sample-col File_name `
  --batch-col Batch `
  --donor-col Donor_id `
  --age-col Age `
  --sex-col Sex `
  --out data/derived/gse164378/biological_replicate_metadata.csv
```

The command rejects conflicting cell keys, within-replicate metadata
conflicts, unmatched cells, and output paths under `data/raw`. The main
GSE164378 workflow joins the authoritative cell-level source table directly;
this replicate-level file is an audit/export and is not a replacement for the
cell-level join table.

## Data Source Notes (Important)
- Your current `raw_counts_h5ad.tar.gz` content is donor-rich but does not include age directly in `obs`.
- For age analyses, you need an additional donor-level clinical file (age/sex/etc.).
- The Immunity 2023 healthy blood atlas is distributed via Synapse/AWS resources; download the companion donor metadata/clinical file and pass it as `--supp`.

### Getting Supplementary Donor Metadata
Synapse access is optional and credentialed. `synapseclient` is intentionally
not required for the core analysis environment unless you are acquiring data
from Synapse. Use a Synapse profile or a local `SYNAPSE_AUTH_TOKEN`
environment variable; do not store tokens in the repository.

Dry-run metadata discovery:
```powershell
python -m src.fetch_synapse_metadata --entity-id syn56693935 --recursive --limit 20
```

Download only metadata-like matches after reviewing the dry-run output:
```powershell
python -m src.fetch_synapse_metadata --entity-id syn56693935 --recursive --download --limit 5
```

Before using a newly downloaded table, identify and document its cell join key
and authoritative biological-replicate field. Do not substitute a donor-like
column based only on its name or apparent uniqueness.
