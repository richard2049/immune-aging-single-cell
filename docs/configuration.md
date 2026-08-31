# Configuration Profiles

Tracked profiles are reproducible execution contracts. Do not edit them in
place for exploratory runs.

| Profile | Role | Workflow | Output |
| --- | --- | --- | --- |
| `config/demo.yaml` | Small fixture and CI DAG check | `workflows/Snakefile` | `results/` |
| `config/blood_age_atlas_pilot.yaml` | Bounded 50,000-cell study run | `workflows/Snakefile` | `results/blood_age_atlas_pilot/` |
| `config/blood_age_atlas_1m.yaml` | One-million-cell checkpoint | `workflows/Snakefile` | `results/blood_age_atlas_1m/` |
| `config/blood_age_atlas_longitudinal.yaml` | Subject-aware downstream inference | `workflows/Snakefile.longitudinal` | `results/blood_age_atlas_longitudinal/` |
| `config/custom.example.yaml` | Template for another standardized H5AD | `workflows/Snakefile` | `results/custom/` |

The Blood Age Atlas is the Terekhova et al. cohort distributed as
`syn49637038`. Paths beginning with `results/gse164378` are historical and must
not be used as output targets.

## Longitudinal Unit Contract

The maintained downstream profile defines:

- `subject_id` from `Donor_id`: dependency and cross-validation group;
- `sample_unit_id` from `Tube_id`: longitudinal observation and aggregation unit;
- `technical_library_id` from `File_name`: technical input identifier.

The mapping must contain 166 subjects and 317 sample units in the complete
study metadata, with one subject, age, sex, and batch assignment per sample
unit. These are executable acceptance checks, not filename assumptions.

## scVI Execution Contract

Every supported profile explicitly records the scVI architecture and training
arguments. The one-million-cell profile fixes the eight epochs produced by the
accepted `scvi-tools 1.3.3` heuristic instead of allowing future package
versions to change the run silently. Hardware overrides may change only
`accelerator` and `devices`.

After reconstruction, `src.validate_scvi_checkpoint` must pass before
clustering. It checks cell and gene identity and order, latent dimensions and
finite values, exact configured parameters, package provenance, and persisted
model files.

## Clustering Execution Contract

Supported profiles explicitly fix neighbours, UMAP, and Leiden arguments.
Clustering uses `leidenalg`, resolution 0.8, seed 42, and two optimization
iterations. The two-iteration bound was approved after the historical
convergence setting remained active for more than four hours on the complete
one-million-cell graph. It preserves the backend and resolution but may change
cluster assignments; candidate clusters therefore require technical audit and
comparison with the historical partition before interpretation.

`src.validate_cluster_checkpoint` is the executable gate before annotation. It
requires unchanged ordered cell and gene identities, finite UMAP and sparse
graph outputs, complete Leiden labels, and exact agreement between stored and
configured clustering parameters. Historical ARI and NMI are reported only as
diagnostics; they are not acceptance thresholds.

## Annotation And Design Qualification

The longitudinal profile preserves CellTypist output in `cell_type` and uses
the versioned mapping at
`config/annotation_mappings/blood_age_atlas_celltypist_v1.csv` to propose
`cell_type_analysis`. Each raw label has a primary, exploratory, or unresolved
disposition plus rationale and review status. Mapping decisions are made
without inspecting age associations.

Run `annotation_design_qualification` to generate annotation support,
raw-count marker summaries, cluster concordance, and method-specific design
eligibility. This target may complete while the mapping is pending so that a
reviewer can inspect the evidence. Downstream rules require the separate
`annotation_design_approved.json` sentinel; it is generated only when the
fingerprinted mapping and qualification report are both approved. Passing a
design check establishes support and estimability, not biological validity.

The review package includes raw-label marker evidence in
`annotation_raw_label_marker_summary.csv` and a one-row-per-label index in
`annotation_mapping_review_table.csv`. Generated evidence is read-only. Follow
`docs/annotation_mapping_review_guide.md` and record human judgements in
`docs/reviews/annotation_mapping_review.yml`. The review must also approve the
composition denominator and the method-specific eligibility contract before
downstream regeneration.

## Inheritance

`config/blood_age_atlas_longitudinal.yaml` extends
`config/blood_age_atlas_1m.yaml`. `src/utils.py` resolves the base path relative
to the child profile and deep-merges mappings. Each tracked profile records its
own filename in `cfg_path`; tests require this to match.

## Custom Data

Create an ignored local profile:

```powershell
Copy-Item config/custom.example.yaml config/custom.local.yaml
```

Set a distinct `cfg_path` and `project.out_dir`. Define subject, sample,
technical-library, batch, and cohort fields from authoritative metadata. A
parsed filename or donor-like label is not sufficient evidence of independence.
Validate join keys, cardinality, missingness, genome/annotation provenance,
covariates, and model support as described in [the real-data guide](real_data.md).

## CPU And GPU

The tracked profiles use CPU for portability. `environment-gpu.yml` provides a
separately qualified NVIDIA CUDA environment. A local `.local.yaml` profile may
override only `scvi.accelerator`, `devices`, and its output directory after a
bounded CPU/GPU smoke test. Hardware selection must not change the seed, input
inclusion, covariates, layer, model settings, or statistical units.

Changing sample identity, inclusion, covariates, model formulas, or inferential
thresholds is an analysis-design change and requires scientific review plus a
new output namespace.
