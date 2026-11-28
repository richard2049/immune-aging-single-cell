# Documentation

This index separates operational guidance, analysis design, and reviewed
evidence so that the workflow's methods and limitations can be inspected
without relying on generated result directories.

## Running the workflow

- [Project overview and quickstart](../README.md)
- [Reproduction runbook and validation evidence](reproducibility.md)
- [Configuration profiles](configuration.md)
- [Real-data ingestion and metadata](real_data.md)
- [Troubleshooting](troubleshooting.md)

## Analysis design

- [Analysis decisions](analysis_decisions.md)
- [Scientific verification contract](verification_contract.md)
- [Annotation mapping review](annotation_mapping_review_guide.md)
- [Completed annotation mapping review](reviews/annotation_mapping_review.yml)
- [Pseudobulk differential-expression design](pseudobulk_de_design.md)
- [Pseudobulk differential-expression workflow](pseudobulk_de.md)

## Reviewed results

The public results cover composition, internal age prediction, and one CD8
RNA signature. Each includes the estimand, uncertainty, and interpretation
limits. [Results and aggregate provenance](../RESULTS.md)

## Technical qualification

The reconstructed scVI, clustering, and CellTypist checkpoints passed their
configured identity, provenance, and output gates. The approved annotation
mapping then qualified ten broad primary populations for method-specific use.
Subject-aware composition, signature, prediction, sensitivity, and pseudobulk
stages have completed, and the bounded and full `dream` runs pass technical
validation. Pseudobulk gene and pathway interpretation and the remaining signature and
sensitivity outputs await scientific review.

## Superseded Historical Evidence

- [Former edgeR runtime qualification](validation/pseudobulk_de_runtime_qualification.md)
- [Former one-million-cell pseudobulk acceptance](validation/pseudobulk_de_1m_run_acceptance.md)

These two records document technical execution of the superseded analysis.
That analysis treated sample units as independent and used an incorrect
study identifier. The [analysis decisions](analysis_decisions.md) explain
the corrected study identity and subject-level design used for current results.

Raw data, generated results, trained models, and local execution evidence are
not versioned. Selected figures used by tracked documentation are stored in
[`docs/assets/`](assets/).
