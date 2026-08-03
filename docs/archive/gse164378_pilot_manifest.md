# Archived Result Manifest: GSE164378 Pilot

## Archive Status

- Original path: `results/gse164378_pilot/`
- Disposition: Deleted after recording this manifest
- Deletion date: 2026-07-27
- Logical size before deletion: 51.66 GiB
- File count before deletion: 30
- Output date range: 2026-02-17 to 2026-02-19

## Historical Purpose

**Inferred:** this directory was an early real-data pilot used to validate the
large-dataset workflow, scVI integration, CellTypist annotation, reporting,
and initial age-analysis modules before the full one-million-cell run.

It is not a current scientific result set.

## Final Annotated Object

- File: `06_annotated.h5ad`
- Shape: 192,317 cells x 27,678 genes
- Cell-type labels: 70
- Stored representations: `X_scVI`, `X_umap`
- Source-label fields: `donor_id_x`, `donor_id_y`
- Canonical `donor_id`: absent
- `biological_replicate_id`: absent

## Reason For Removal

- The stored 192,317-cell run does not match the current
  `config/gse164378_pilot.yaml` setting of `max_cells: 50000`, so the exact run is
  not reproducible from the current profile.
- Donor-level outputs used non-canonical donor labels and predate the validated
  `Tube_id` to `biological_replicate_id` remediation.
- No current workflow, documentation asset, or corrected analysis depends on
  this result directory.
- The full annotated checkpoint and corrected donor-level results remain
  available under `results/gse164378_full/` and
  `results/gse164378_full_replicate_corrected/`.

## Removed Contents

- Seven workflow H5AD files from `01_raw.h5ad` through
  `06_annotated.h5ad`.
- One scVI model (`models/scvi_model/model.pt`, approximately 70 MiB).
- Pilot report, composition, signature, and age-prediction figures and tables.
- `reports.done`.
