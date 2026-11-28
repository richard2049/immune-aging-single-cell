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

## Build Unit Metadata

`src.biological_replicates` builds the maintained cell-to-unit mapping from an
authoritative cell-level table. Configure source and canonical columns under
`biological_units`; do not infer independence from uniqueness or naming. The
mapping rejects duplicate cell keys, missing required values, sample-to-subject
conflicts, inconsistent age/sex/batch within a sample, and unexpected cohort
counts. Derived mappings must be written outside `data/raw`.

## Blood Age Atlas inputs

The [source study](https://pubmed.ncbi.nlm.nih.gov/37963457/) describes
317 blood-sample units from 166 people aged 25?85. The
[Synapse study record](https://www.synapse.org/Synapse:syn49637038) and
[AWS registry record](https://registry.opendata.aws/singlecellhumanbloodatlas/)
describe the distribution and access conditions.

| Archive | Required local file in the tracked profile |
| --- | --- |
| [Raw counts, syn56693935](https://www.synapse.org/Synapse:syn56693935) | `data/raw/raw_counts_h5ad/pbmc_gex_raw_with_var_obs.h5ad` |
| [Companion PBMC archive, syn51197006](https://www.synapse.org/Synapse:syn51197006) | `data/raw/raw_counts_h5ad/all_pbmcs/all_pbmcs_metadata.csv` |

Follow Synapse's access conditions and compare downloads to the checksums
provided for the selected file versions where available. Keep original
archives and credentials outside Git. The core workflow reads local files.
The optional download helper uses a Synapse profile or a local
`SYNAPSE_AUTH_TOKEN` and can list metadata candidates before any download.

```powershell
python -m src.fetch_synapse_metadata --entity-id syn51197006 --recursive --limit 20
```

The accepted profile joins cell identifiers to companion metadata to recover
`Donor_id`, `Tube_id`, and `File_name`. Library identifiers record technical
provenance. A multiplexed library may contain cells from multiple sample
units and people. Check the archive versions, join key, and unit columns
before using replacement inputs. [Profile and schema](configuration.md#longitudinal-unit-contract)

### Cell retention in the accepted profile

The archived count matrix contained 1,916,367 cells. The profile selected one
million cells and retained them through metadata, scVI, clustering, and
annotation. The source-study code uses genotype-based demultiplexing with
souporcell and includes downstream preprocessing. The exact upstream QC and
doublet filters applied to this specific deposited count matrix have
incomplete file-level provenance.
[Source demultiplexing workflow](https://github.com/teresho4/scRNA-seq_atlas_Hs_PBMC_aging/blob/main/Demultiplexing_pipeline/Snakefile)

The accepted local profile disabled QC filtering and SOLO doublet detection.
The intermediate checkpoint names describe workflow positions. The
composition denominator contains 999,766 retained cells across the 316
eligible visits, including cells labeled `Other/unresolved`. Signature,
pseudobulk, and prediction models apply their own support rules.
[Reviewed results and support](../RESULTS.md)
