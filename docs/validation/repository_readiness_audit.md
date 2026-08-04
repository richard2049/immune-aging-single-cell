# Repository Readiness Audit

## Audit Decision

**Ready to resume bounded scientific development, with interpretation and
external-validation gates still open.**

This decision applies to repository engineering, reproducible execution, and
continued analysis development. It does not certify the study as biologically
complete, externally validated, causal, longitudinal, diagnostic, or
clinically useful.

- Audit date: 2026-08-03
- Audited baseline: `main` after PR 5 (`78dfc61`)
- Repository: [richard2049/immune-aging-single-cell](https://github.com/richard2049/immune-aging-single-cell)

## Scope

The audit evaluated the public Git tree, documentation, environment contract,
workflow entry points, automated tests, CI configuration, scientific review
boundaries, citation and license metadata, and the documented maintenance
route. It did not rerun the one-million-cell analysis, retrain scVI, or repeat
the complete edgeR and robustness workloads. Those expensive stages retain their
existing technical acceptance evidence and remain available through the
documented maintenance commands.

## Readiness Evidence

| Area | Outcome | Evidence and boundary |
| --- | --- | --- |
| Local and remote state | Pass | Local `main` and `origin/main` matched at the audited baseline before this report was prepared. |
| Public repository metadata | Pass | Repository is public, uses `main`, exposes a clear scientific description, and GitHub detects the BSD 3-Clause license. |
| Repository hygiene | Pass | No raw data, generated `results/`, environments, caches, Snakemake state, private working records, credentials, or absolute local paths were found in the tracked tree. |
| Citation and release | Pass | `CITATION.cff` parses with version `0.1.0`, BSD-3-Clause, author metadata, and ORCID. Release `v0.1.0` remains the Phase 1 checkpoint. |
| Documentation | Pass | Sixteen tracked Markdown files, including this report, had valid relative links. The README now distinguishes demo execution, incremental maintenance, one-million-cell reconstruction, and checkpoint-based corrected regeneration. |
| Runtime environment | Pass with residual risk | `environment.yml` is version-constrained and was created successfully on GitHub-hosted Linux. It is not a platform lockfile, so future dependency resolution can still vary within allowed patch versions. |
| Automated verification | Pass | [PR CI run 30819119516](https://github.com/richard2049/immune-aging-single-cell/actions/runs/30819119516) and [main CI run 30819522816](https://github.com/richard2049/immune-aging-single-cell/actions/runs/30819522816) passed Ruff, formatting, compilation, 42 unit tests, and demo-DAG construction. Two Docker tests were skipped by the routine CI contract. |
| CI security and scope | Pass | Workflow permissions are read-only, actions are pinned by commit SHA, superseded runs are cancelled, and heavy or external qualifications are explicitly excluded from routine CI. |
| Demo workflow graph | Pass | The 13-rule demo DAG constructs without executing jobs or requiring real data. |
| One-million-cell maintenance route | Pass, not re-executed | The README runs the one-million-cell profile first and then explicitly targets both corrected default outputs and `pseudobulk_robustness_evidence`. Complete recomputation remains an intentional long-running operation. |
| Replicate correction | Pass | The corrected workflow uses `Tube_id` as `biological_replicate_id`; its mapping, grouped analyses, validation, and reviewed presentation scope are documented in the [corrected-results review](../corrected_results_review.md). |
| Pseudobulk implementation | Technically accepted | The donor-aware edgeR run passed its [one-million-cell acceptance contract](pseudobulk_de_1m_run_acceptance.md), and the predefined [robustness audit](../pseudobulk_robustness_audit.md) completed without model failures. |
| Gene and pathway interpretation | Open | The bounded pseudobulk queue supports human review but does not approve genes, pathways, or mechanisms. |
| External validation | Open | Corrected composition, signature, and predictive findings have not been independently replicated in a compatible donor cohort. |
| Branch policy | Recommendation | `main` is not protected. CI reports status but does not currently prevent direct pushes or merging with failed checks. |

## Maintenance Contract

The authoritative GSE164378 rebuild is intentionally sequential:

1. Run the one-million-cell profile to produce
   `results/gse164378_full/06_annotated.h5ad` and provisional outputs.
2. Run the replicate-corrected workflow with both `all` and
   `pseudobulk_robustness_evidence` as targets.
3. Inspect technical validation files, logs, and tracked-document diffs before
   accepting regenerated outputs.

The copy-pasteable PowerShell procedure, preflight checks, failure gates, and
logging location are provided in the
[README maintenance section](../../README.md#maintenance-rebuilds).

Routine maintenance should rely on Snakemake's incomplete and outdated-job
detection. `--forceall` is appropriate only when a complete recomputation is
deliberately required. Validated raw inputs and large checkpoints should remain
immutable and outside Git.

## Scientific Boundaries

The following remain mandatory when scientific work resumes:

- use biological replicate or donor-level inference for age-related questions;
- preserve the approved metadata mapping and grouped cross-validation rules;
- treat the current age associations as cross-sectional;
- describe age prediction as internal and transductive until the complete
  representation-learning procedure is tested in unseen donors;
- keep automated pseudobulk ranking separate from biological disposition;
- require human review before adding gene-, pathway-, causal-, diagnostic-, or
  clinical claims;
- use independent data for external validation rather than selecting results
  for agreement with prior literature.

## Remaining Actions

1. Consider protecting `main` and requiring the two CI jobs for pull requests.
2. Review and disposition the bounded pseudobulk evidence queue before any
   gene- or pathway-level presentation.
3. Plan independent-cohort validation with compatible donor, age, sex, batch,
   and cell-type metadata.
4. Consider a platform lockfile if exact long-term environment reconstruction
   becomes a release requirement.

The first action is repository governance. The remaining actions are new
scientific or reproducibility stages and require their own predefined scope and
acceptance criteria.
