# Pseudobulk Robustness Audit And Evidence Preparation

> [!IMPORTANT]
> This audit applies to historical outputs whose annotated checkpoint is now
> pending provenance-complete requalification. Its gene-review queue must not be
> interpreted or promoted until pseudobulk outputs are regenerated and the
> technical audit is repeated.

## Scope

This stage evaluates technical and statistical stability of the frozen
donor-aware pseudobulk results. It does not approve genes, pathways,
biological interpretations, or public claims.

## Audit Status

- Automated contract: passed
- Global-FDR candidate rows preserved: 7,970 (7,961 primary; 9 exploratory)
- Compact manual-review queue: 249
- Completed sensitivity models: 212
- Sensitivity execution record: validated outputs; wrapper log incomplete after Docker client hang
- Diagnostic plot review: complete; follow-up flags recorded
- Biological interpretation: not performed

## Sensitivity Model Inventory

| Scenario type | Status | Models |
| --- | --- | ---: |
| covariate | completed | 26 |
| leave one batch out | completed | 175 |
| leave one batch out | not estimable | 1 |
| support | completed | 11 |
| support | not estimable | 2 |

## Non-Estimable Scenarios

| Cell type | Tier | Scenario | Replicates | Residual df | Reason |
| --- | --- | --- | ---: | ---: | --- |
| Memory B cells | primary | high cell support | 13 | 2 | insufficient residual df |
| NK cells | exploratory | high cell support | 3 | NA | insufficient replicates |
| NK cells | exploratory | leave one batch out: AS053 | 15 | NA | model matrix error: contrasts can be applied only to factors with 2 or more levels |

## Cell-Type Audit Summary

| Cell type | Tier | Replicates | Residual df | Global-FDR candidates | Completed scenarios | Median direction concordance | MD plot review |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| NK cells | exploratory | 18 | 8 | 9 | 9 | 1.000 | reviewed; follow-up required |
| CD16+ NK cells | primary | 312 | 296 | 666 | 17 | 1.000 | reviewed; no blocking issue |
| Classical monocytes | primary | 312 | 296 | 1000 | 17 | 1.000 | reviewed; no blocking issue |
| MAIT cells | primary | 113 | 97 | 110 | 17 | 1.000 | reviewed; follow-up required |
| Memory B cells | primary | 105 | 89 | 151 | 16 | 1.000 | reviewed; follow-up required |
| Naive B cells | primary | 175 | 159 | 52 | 17 | 1.000 | reviewed; follow-up required |
| Non-classical monocytes | primary | 231 | 215 | 353 | 17 | 1.000 | reviewed; no blocking issue |
| Regulatory T cells | primary | 247 | 231 | 236 | 17 | 1.000 | reviewed; no blocking issue |
| Tcm/Naive cytotoxic T cells | primary | 281 | 265 | 1446 | 17 | 1.000 | reviewed; no blocking issue |
| Tcm/Naive helper T cells | primary | 317 | 301 | 2592 | 17 | 1.000 | reviewed; no blocking issue |
| Tem/Effector helper T cells | primary | 312 | 296 | 695 | 17 | 1.000 | reviewed; no blocking issue |
| Tem/Temra cytotoxic T cells | primary | 246 | 230 | 331 | 17 | 1.000 | reviewed; no blocking issue |
| Tem/Trm cytotoxic T cells | primary | 275 | 259 | 329 | 17 | 1.000 | reviewed; no blocking issue |

## Runtime Provenance

| Component | Expected | Observed |
| --- | --- | --- |
| R | 4.6 | 4.6.1 |
| Bioconductor | 3.23 | 3.23 |
| edgeR | 4.10.1 | 4.10.1 |
| Matrix | not pinned | 1.7.5 |

## Interpretation Boundary

Sensitivity models reuse the same cohort and are not independent
replications. Direction concordance and effect ranges are review
diagnostics, not automated evidence of biological validity. Omission
of sex or batch deliberately weakens the adjustment set and is used
only to expose adjustment dependence. Leave-one-batch-out models can
identify sensitivity to a processing batch but cannot remove all
residual confounding.

The exploratory NK-cell population remains separate because its
reference model used only 18 qualifying replicates and eight residual
degrees of freedom. Gene annotation, pathway enrichment, candidate
disposition, external replication, and public wording require a separate human
scientific review.

## Operational Notes

- The Docker container wrote the complete sensitivity outputs, but the Windows Docker client did not close its wrapper log. The two identified host processes were stopped after the container exited; acceptance relies on independent schema, key, count, scenario, and model-status checks.
