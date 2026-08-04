# Pseudobulk DE Feasibility And Proposed Analysis Contract

## Status

Feasibility audit, biological-replicate remediation, analysis-contract
approval, and bounded runtime qualification are complete. The donor-aware
pseudobulk implementation is available as a targeted corrected-profile stage.
The one-million-cell GSE164378 pseudobulk analysis and automated technical
validation are complete. Biological interpretation remains pending.

## Audit Scope

The audit used:

- profile: `config/gse164378_1m.yaml`;
- checkpoint: `results/gse164378_full/06_annotated.h5ad`;
- source cell metadata:
  `data/raw/raw_counts_h5ad/all_pbmcs/all_pbmcs_metadata.csv`;
- command:

```powershell
python -u -m src.pseudobulk_de_feasibility `
  --config config/gse164378_1m.yaml `
  --inp results/gse164378_full/06_annotated.h5ad `
  --outdir results/gse164378_full/pseudobulk_de_feasibility
```

The command scans the sparse count values and gene indices by row blocks. It
does not modify the checkpoint or convert the cell-by-gene matrix to a dense
array.

## Observed Findings

### Count matrix

- `adata.X` is a CSR matrix with 1,000,000 cells, 36,601 genes, and
  1,252,929,694 stored nonzero values.
- A complete scan found no non-finite, negative, fractional, or out-of-range
  values.
- The matrix therefore passes the raw-count integrity checks required for
  pseudobulk aggregation.
- The scan used an estimated peak of about 94 MB for count/index chunks and
  row-library-size state.

### Biological replicate identity

- `donor_id` is not globally unique: 166 donor labels are reused across the
  source pools.
- Grouping by `donor_id` mixes different ages for 102 labels and different
  batches for 99 labels.
- `sample_id` identifies a multiplexed technical library, not an individual.
  One library contains up to six biological replicates, and one biological
  replicate appears in up to three libraries.
- `Tube_id` identifies 317 biological replicates. It maps to one donor label,
  age, sex, and batch without conflicts and is present for all 1,000,000
  audited cells through the source metadata join.

At the time of this audit, donor-level analyses configured with
`donor_col: donor_id` merged distinct biological replicates. Those analyses
were subsequently regenerated with `biological_replicate_id`, compared with
the provisional outputs, and reviewed before the corrected README results were
approved. The original finding remains relevant to the pseudobulk design
because the same corrected experimental unit must be used.

### Candidate cell-type support

Using candidate thresholds of at least 50 cells per replicate x cell type,
at least 12 biological replicates, at least 12 years of age coverage, and at
least five residual model degrees of freedom:

| Cell type | Replicates | Age span | Residual df | Candidate status |
| --- | ---: | ---: | ---: | --- |
| Tcm/Naive helper T cells | 317 | 56 | 301 | Eligible |
| Classical monocytes | 312 | 56 | 296 | Eligible |
| CD16+ NK cells | 312 | 56 | 296 | Eligible |
| Tem/Effector helper T cells | 312 | 56 | 296 | Eligible |
| Tcm/Naive cytotoxic T cells | 281 | 56 | 265 | Eligible |
| Tem/Trm cytotoxic T cells | 275 | 56 | 259 | Eligible |
| Regulatory T cells | 247 | 56 | 231 | Eligible |
| Tem/Temra cytotoxic T cells | 246 | 56 | 230 | Eligible |
| Non-classical monocytes | 231 | 56 | 215 | Eligible |
| Naive B cells | 175 | 56 | 159 | Eligible |
| MAIT cells | 113 | 48 | 97 | Eligible |
| Memory B cells | 105 | 55 | 89 | Eligible |
| NK cells | 18 | 54 | 8 | Eligible, but borderline |

The other 71 labels fail at least one candidate support or design criterion.
The full exclusion reasons are recorded in
`celltype_support_summary.csv`. Candidate eligibility is not equivalent to
approval for final inference; NK cells in particular have limited residual
degrees of freedom after covariate adjustment.

## Approved Analysis Contract

The following contract was approved on 2026-07-30:

1. Define `biological_replicate_id` from `Tube_id` during metadata integration.
2. Sum raw counts over all technical libraries within each
   `biological_replicate_id x cell_type` group.
3. Preserve age as a continuous variable and report the age coefficient as
   log2 fold change per 10 years.
4. Fit one model per approved cell type using `age + sex + batch`, but exclude
   or redesign any cell type whose adjusted design is rank-deficient or has
   inadequate residual degrees of freedom.
5. Use `edgeR` quasi-likelihood with TMM normalization and
   design-aware low-expression filtering.
6. Report both within-cell-type BH FDR and a global BH FDR over all tested
   gene-cell-type combinations. Reserve cross-cell-type headline claims for
   the global correction.
7. Store pseudobulk counts, replicate metadata, gene-level results, exclusion
   reasons, package versions, and the exact design formula.

The 12 well-supported populations are primary. NK cells are processed under
the same model but remain exploratory because they have only 18 qualifying
replicates and eight residual degrees of freedom. If the adjusted model is not
estimable, the cell type is excluded from primary adjusted inference and an
age-only result is reported as sensitivity analysis.

## Required Remediation Sequence

Before implementing DE:

1. Add `tube_id` to the integrated metadata and expose a canonical
   `biological_replicate_id`.
2. Change all donor-level analysis configurations to use that canonical key.
3. Regenerate composition, signature, age-prediction, sensitivity, and
   supplementary outputs from the corrected checkpoint.
4. Re-audit the corrected checkpoint and compare old versus corrected
   conclusions explicitly.
5. Add and validate the selected count-model runtime (completed).
6. Implement pseudobulk aggregation and DE only after the revised design is
   approved (completed).
7. Execute the complete corrected-profile target and review technical diagnostics
   before any biological interpretation (automated acceptance completed;
   manual diagnostic-plot review and biological interpretation pending).

## Generated Audit Outputs

Under `results/gse164378_full/pseudobulk_de_feasibility/`:

- `feasibility_summary.json`
- `matrix_audit.json`
- `replicate_metadata_audit.csv`
- `replicate_celltype_support.csv`

The implemented one-million-cell stage writes under
`results/gse164378_full_replicate_corrected/pseudobulk_de/`. Its output
contract and targeted commands are documented in `docs/pseudobulk_de.md`.
- `celltype_support_summary.csv`
- `design_diagnostics.csv`
- `cell_library_sizes.npy`
- `audit.log`

These are generated artifacts and remain outside version control.
