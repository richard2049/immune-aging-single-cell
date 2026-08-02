# Decision Question Bank

This file is the compact index of questions recorded in decision cards. It is
not a general brainstorming list: add a question only when it exists in a
decision card.

Keep one canonical row per distinct question. If the same question appears in
more than one card, keep one row and include all relevant card links in the
final column. Use the stage and priority from the source card metadata and the
status recorded for the question.

## Interview Questions

| Question | Stage | Priority | Status | Skill tested | Decision card |
| --- | --- | --- | --- | --- | --- |
| Why should age-associated expression in this single-cell dataset be tested with donor-level pseudobulk profiles rather than by treating cells as independent replicates? | Donor-aware pseudobulk DE design and implementation | High | Ready | Experimental-unit identification, pseudoreplication, pseudobulk aggregation, and count-based differential-expression design | [Executed and technically validated pseudobulk DE](decisions/001_pseudobulk_de_feasibility_and_design.md) |
| How would you correct a reused donor identifier in an established single-cell workflow without unnecessarily retraining upstream models? | Biological replicate remediation and donor-level result regeneration | High | Ready | Experimental-unit definition, metadata provenance, checkpoint reuse, leakage prevention, and reproducible result migration | [Completed biological replicate remediation](decisions/002_biological_replicate_remediation.md) |
| How would you decide whether corrected single-cell age associations are ready for public presentation after changing the experimental unit? | Human review of corrected immune-aging results and public promotion | High | Ready | Scientific result triage, effect-size interpretation, multiple-testing awareness, sensitivity analysis, confounding assessment, and cautious communication | [Completed corrected results review and promotion](decisions/003_corrected_results_interpretation_and_promotion.md) |
| Which workflow checks prevent metadata and environment differences from silently changing a single-cell analysis across machines? | Metadata and model reproducibility remediation | High | Ready | Metadata cardinality, experimental-unit validation, dependency contracts, stochastic reproducibility, and interpretation boundaries | [Implemented and verified scientific reproducibility remediation](decisions/004_scientific_reproducibility_remediation.md) |
| How would you audit robustness of donor-level pseudobulk differential-expression results before interpreting the top genes? | Pseudobulk robustness audit and evidence preparation | High | Ready | Sensitivity design, confounding assessment, experimental-unit preservation, multiple-testing awareness, and evidence triage | [Implemented pseudobulk robustness audit](decisions/005_pseudobulk_robustness_and_evidence.md) |

## Follow-up Biological Questions

| Question | Stage | Priority | Status | Biological theme | Decision card |
| --- | --- | --- | --- | --- | --- |
| Which genes and pathways show donor-aware age associations within sufficiently supported immune cell populations, beyond composition shifts and predefined signature scores? | Donor-aware pseudobulk DE design and implementation | High | Open | Cell-type-specific transcriptional immune aging | [Executed and technically validated pseudobulk DE](decisions/001_pseudobulk_de_feasibility_and_design.md) |
| Which provisional immune-aging findings remain stable after regrouping cells by the correct biological replicate? | Biological replicate remediation and donor-level result regeneration | High | Open | Robustness of immune-aging associations to replicate definition | [Biological replicate remediation and result regeneration](decisions/002_biological_replicate_remediation.md) |
| Which corrected immune-aging associations reproduce in an independent donor cohort with comparable cell populations and metadata? | Human review of corrected immune-aging results and public promotion | High | Open | External validity and reproducibility of immune-aging associations | [Completed corrected results review and promotion](decisions/003_corrected_results_interpretation_and_promotion.md) |
| Do the retained immune-aging findings remain stable when the full representation-learning and prediction procedure is evaluated on unseen donors? | Metadata and model reproducibility remediation | High | Open | External validity of latent immune-aging representations | [Implemented and verified scientific reproducibility remediation](decisions/004_scientific_reproducibility_remediation.md) |
| Which cell-type-specific age-associated genes retain comparable effect direction and magnitude across support and batch-sensitivity analyses? | Pseudobulk robustness audit and evidence preparation | High | Open | Robustness of cell-type-specific transcriptional immune-aging associations | [Implemented pseudobulk robustness audit](decisions/005_pseudobulk_robustness_and_evidence.md) |
