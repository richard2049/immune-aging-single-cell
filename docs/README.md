# Documentation

This index separates operational guidance, analysis design, and reviewed
evidence so that the workflow's methods and limitations can be inspected
without relying on generated result directories.

## Running the workflow

- [Project overview and quickstart](../README.md)
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

## Current Validation

The reconstructed scVI, clustering, and CellTypist checkpoints passed their
configured identity, provenance, and output gates. The approved annotation
mapping then qualified ten broad primary populations for method-specific use.
Subject-aware composition, signature, prediction, sensitivity, and pseudobulk
stages have completed, and the bounded and full `dream` runs pass technical
validation. These checks establish traceability and technical validity, but
biological interpretation of the regenerated results is still pending.

## Superseded Historical Evidence

- [Historical GSE164378 pilot manifest](archive/gse164378_pilot_manifest.md)
- [Former corrected-results review](corrected_results_review.md)
- [Former pseudobulk robustness audit](pseudobulk_robustness_audit.md)
- [Former edgeR runtime qualification](validation/pseudobulk_de_runtime_qualification.md)
- [Former one-million-cell pseudobulk acceptance](validation/pseudobulk_de_1m_run_acceptance.md)
- [Former result dispositions](reviews/corrected_results_dispositions.yml)
- [Former public-promotion record](reviews/corrected_public_promotion_record.md)

These records preserve provenance for the superseded analysis that treated
sample units as independent. They are historical only and should not be used to
interpret results from the current workflow.

Raw data, generated results, trained models, and local execution evidence are
not versioned. Selected figures used by tracked documentation are stored in
[`docs/assets/`](assets/).
