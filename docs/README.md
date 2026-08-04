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
- [Pseudobulk differential-expression design](pseudobulk_de_design.md)
- [Pseudobulk differential-expression workflow](pseudobulk_de.md)

## Validation and reviewed evidence

- [Corrected-results review](corrected_results_review.md)
- [Pseudobulk robustness audit](pseudobulk_robustness_audit.md)
- [edgeR runtime qualification](validation/pseudobulk_de_runtime_qualification.md)
- [One-million-cell pseudobulk acceptance](validation/pseudobulk_de_1m_run_acceptance.md)
- [Demo workflow end-to-end qualification](validation/demo_end_to_end.md)
- [Repository readiness audit](validation/repository_readiness_audit.md)
- [Corrected-result dispositions](reviews/corrected_results_dispositions.yml)
- [Pseudobulk diagnostic review](reviews/pseudobulk_diagnostic_plot_review.yml)
- [Corrected public-promotion record](reviews/corrected_public_promotion_record.md)

## Provenance and archived material

- [Historical GSE164378 pilot manifest](archive/gse164378_pilot_manifest.md)

Raw data, generated results, trained models, and local execution evidence are
not versioned. Selected figures used by tracked documentation are stored in
[`docs/assets/`](assets/).
