# Real Data Guide

This project expects a pre-standardized `.h5ad` for real runs.

## Minimum AnnData Requirements
- Matrix in `adata.X` with raw counts preferred.
- Unique identifiers:
  - `adata.obs_names` unique cell barcodes
  - `adata.var_names` unique gene symbols/IDs
- Recommended columns in `adata.obs`:
  - `batch` (strongly recommended for scVI integration)
  - `subject_id` from authoritative participant metadata
  - `sample_unit_id` for the specimen, visit, or extraction being analysed
  - `technical_library_id` for library/file provenance
  - `age` (required for age-stratified analyses)
  - `sex` (optional)
  - `condition` (optional)

## Configure Real Mode
Do not edit a tracked Blood Age Atlas study profile. Copy the generic example to an
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
  parsed labels do not establish the subject or independent unit
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

## Build Unit Metadata

`src.biological_replicates` builds the maintained cell-to-unit mapping from an
authoritative cell-level table. Configure source and canonical columns under
`biological_units`; do not infer independence from uniqueness or naming. The
mapping rejects duplicate cell keys, missing required values, sample-to-subject
conflicts, inconsistent age/sex/batch within a sample, and unexpected cohort
counts. Derived mappings must be written outside `data/raw`.

## Data Source Notes (Important)
- The current Blood Age Atlas count archive does not include every required
  longitudinal field directly in `obs`.
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

Before using a newly downloaded table, identify and document its cell join key,
subject, sample-unit, and technical-library fields. Do not substitute one role
for another based only on a name or apparent uniqueness.
