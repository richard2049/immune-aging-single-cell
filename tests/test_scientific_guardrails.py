from __future__ import annotations

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd
import yaml

from src.audit_scientific_checkpoint import audit
from src.composition_age import (
    ALL_QC_DENOMINATOR,
    _build_donor_fraction_table,
    _compute_trends,
    _prepare_obs,
    _resolve_composition_cell_types,
)
from src.longitudinal_stats import (
    GEENonConvergenceError,
    fit_independent_age_model,
    fit_subject_aware_gee,
    select_one_sample_per_subject,
)
from src.scientific_guardrails import (
    build_complete_nuisance_design,
    require_placeholder_permission,
    resolve_covariates,
    validate_replicate_covariates,
)
from src.scvi_train import _index_sha256, _resolved_config_sha256
from src.signature_age import _fit_signature_gee, _pick_columns
from src.utils import load_config
from src.validate_scvi_checkpoint import validate as validate_scvi_checkpoint


class ScientificGuardrailTests(unittest.TestCase):
    @staticmethod
    def _temporary_directory() -> TemporaryDirectory:
        root = Path(".tmp/tests")
        root.mkdir(parents=True, exist_ok=True)
        return TemporaryDirectory(dir=root)

    @staticmethod
    def _checkpoint_fixture(root: Path, *, complete_provenance: bool) -> tuple[Path, Path, Path]:
        neighbors_args = {
            "n_neighbors": 15,
            "n_pcs": "null",
            "knn": True,
            "method": "umap",
            "transformer": "null",
            "metric": "euclidean",
            "metric_kwds": {},
        }
        umap_args = {
            "min_dist": 0.5,
            "spread": 1.0,
            "n_components": 2,
            "maxiter": "null",
            "alpha": 1.0,
            "gamma": 1.0,
            "negative_sample_rate": 5,
            "init_pos": "spectral",
            "a": "null",
            "b": "null",
            "method": "umap",
        }
        leiden_args = {
            "directed": "null",
            "use_weights": True,
            "n_iterations": 2,
            "neighbors_key": "null",
            "flavor": "leidenalg",
        }
        obs = pd.DataFrame(
            {
                "age": [30, 30, 50, 50],
                "sex": ["F", "F", "M", "M"],
                "batch": ["B1", "B1", "B2", "B2"],
                "cell_type": ["T", "B", "T", "B"],
                "cell_type_confidence": [0.9, 0.8, 0.7, 0.6],
                "biological_replicate_id": ["r1", "r1", "r2", "r2"],
            },
            index=pd.Index(["c1", "c2", "c3", "c4"], name="cell_id"),
        )
        checkpoint = ad.AnnData(
            X=np.ones((4, 3), dtype=np.float32),
            obs=obs,
            var=pd.DataFrame(index=["g1", "g2", "g3"]),
            obsm={"X_scVI": np.ones((4, 2), dtype=np.float32)},
        )
        checkpoint.uns["clustering_provenance"] = {
            "started_at_utc": "2026-08-16T00:01:00+00:00",
            "completed_at_utc": "2026-08-16T00:02:00+00:00",
            "input_path": "fixture-scvi.h5ad",
            "input_shape": [4, 3],
            "input_obs_names_sha256": "fixture-obs",
            "input_var_names_sha256": "fixture-var",
            "config_path": "fixture.yaml",
            "resolved_config_sha256": "fixture-config",
            "seed": 42,
            "n_jobs": 1,
            "use_rep": "X_scVI",
            "leiden_resolution": 0.8,
            "neighbors_args": neighbors_args,
            "umap_args": umap_args,
            "leiden_args": leiden_args,
            "scanpy_version": "fixture",
        }
        checkpoint.uns["celltypist_provenance"] = {
            "model_identifier": "fixture.pkl",
            "resolved_model_path": "fixture.pkl",
            "model_sha256": "fixture",
            "model_source": "test fixture",
            "celltypist_version": "fixture",
            "majority_voting": False,
            "chunk_size": 10,
            "normalize_target_sum": 10_000,
        }
        if complete_provenance:
            checkpoint.uns["scvi_training_provenance"] = {
                "started_at_utc": "2026-08-16T00:00:00+00:00",
                "completed_at_utc": "2026-08-16T00:01:00+00:00",
                "input_path": "fixture.h5ad",
                "input_shape": [4, 3],
                "input_obs_names_sha256": "fixture-obs",
                "input_var_names_sha256": "fixture-var",
                "config_path": "fixture.yaml",
                "resolved_config_sha256": "fixture-config",
                "seed": 42,
                "layer": "X",
                "categorical_covariates": ["batch"],
                "continuous_covariates": [],
                "allow_missing_covariates": False,
                "model_args": {"n_latent": 2},
                "training_args": {"max_epochs": 1},
                "accelerator": "cpu",
                "devices": 1,
                "scvi_tools_version": "fixture",
                "torch_version": "fixture",
                "anndata_version": "fixture",
                "scanpy_version": "fixture",
            }
        checkpoint_path = root / "checkpoint.h5ad"
        checkpoint.write_h5ad(checkpoint_path)

        source = pd.DataFrame(
            {
                "cell_id": obs.index,
                "tube_id": obs["biological_replicate_id"].to_numpy(),
                "age": obs["age"].to_numpy(),
                "sex": obs["sex"].to_numpy(),
                "batch": obs["batch"].to_numpy(),
            }
        )
        source_path = root / "source.csv"
        source.to_csv(source_path, index=False)
        mapping = source.rename(columns={"tube_id": "biological_replicate_id"})
        mapping_path = root / "mapping.csv"
        mapping.to_csv(mapping_path, index=False)

        config = {
            "run": {"seed": 42, "allow_placeholder_outputs": False},
            "clustering": {
                "n_jobs": 1,
                "use_rep": "X_scVI",
                "leiden_resolution": 0.8,
                "neighbors_args": {
                    key: None if value == "null" else value for key, value in neighbors_args.items()
                },
                "umap_args": {
                    key: None if value == "null" else value for key, value in umap_args.items()
                },
                "leiden_args": {
                    key: None if value == "null" else value for key, value in leiden_args.items()
                },
            },
            "composition_age": {
                "age_col": "age",
                "celltype_col": "cell_type",
                "covariate_cols": ["sex", "batch"],
            },
            "signature_age": {
                "age_col": "age",
                "celltype_col": "cell_type",
                "covariate_cols": ["sex", "batch"],
            },
            "age_prediction": {
                "age_col": "age",
                "celltype_col": "cell_type",
                "latent_key": "X_scVI",
            },
            "biological_replicates": {
                "canonical_col": "biological_replicate_id",
                "source_col": "tube_id",
                "source_table_path": str(source_path),
                "source_table_join_key": "cell_id",
                "age_col": "age",
                "sex_col": "sex",
                "batch_col": "batch",
                "mapping_path": str(mapping_path),
                "strict": True,
            },
            "scientific_audit": {
                "annotation_confidence_col": "cell_type_confidence",
                "stage_h5ads": {"annotated_checkpoint": str(checkpoint_path)},
            },
        }
        config_path = root / "config.yaml"
        config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
        return config_path, checkpoint_path, mapping_path

    def test_checkpoint_audit_accepts_attributable_fixture(self) -> None:
        with self._temporary_directory() as directory:
            config_path, checkpoint_path, mapping_path = self._checkpoint_fixture(
                Path(directory), complete_provenance=True
            )
            report, retention, confidence = audit(
                config_path=config_path,
                checkpoint_path=checkpoint_path,
                mapping_path=mapping_path,
            )

        self.assertTrue(report["passed"])
        self.assertFalse(retention.empty)
        self.assertEqual(len(confidence), 2)

    def test_checkpoint_audit_rejects_incomplete_provenance(self) -> None:
        with self._temporary_directory() as directory:
            config_path, checkpoint_path, mapping_path = self._checkpoint_fixture(
                Path(directory), complete_provenance=False
            )
            report, _, _ = audit(
                config_path=config_path,
                checkpoint_path=checkpoint_path,
                mapping_path=mapping_path,
            )

        self.assertFalse(report["passed"])
        self.assertIn(
            "provenance__scvi_training_provenance",
            report["failed_checks"],
        )

    def test_checkpoint_audit_rejects_clustering_parameter_mismatch(self) -> None:
        with self._temporary_directory() as directory:
            config_path, checkpoint_path, mapping_path = self._checkpoint_fixture(
                Path(directory), complete_provenance=True
            )
            checkpoint = ad.read_h5ad(checkpoint_path)
            provenance = dict(checkpoint.uns["clustering_provenance"])
            leiden_args = dict(provenance["leiden_args"])
            leiden_args["n_iterations"] = 3
            provenance["leiden_args"] = leiden_args
            checkpoint.uns["clustering_provenance"] = provenance
            checkpoint.write_h5ad(checkpoint_path)

            report, _, _ = audit(
                config_path=config_path,
                checkpoint_path=checkpoint_path,
                mapping_path=mapping_path,
            )

        self.assertFalse(report["passed"])
        self.assertIn("clustering_parameters_match_config", report["failed_checks"])

    def test_scvi_checkpoint_gate_preserves_identity_and_contract(self) -> None:
        with self._temporary_directory() as directory:
            root = Path(directory)
            source = ad.AnnData(
                X=np.ones((4, 3), dtype=np.float32),
                obs=pd.DataFrame(
                    {"batch": ["A", "A", "B", "B"]},
                    index=["c1", "c2", "c3", "c4"],
                ),
                var=pd.DataFrame(index=["g1", "g2", "g3"]),
            )
            source_path = root / "03_nodoublets.h5ad"
            source.write_h5ad(source_path)

            model_args = {
                "n_hidden": 128,
                "n_latent": 2,
                "n_layers": 1,
                "dropout_rate": 0.1,
                "dispersion": "gene",
                "gene_likelihood": "zinb",
                "use_observed_lib_size": True,
                "latent_distribution": "normal",
            }
            training_args = {
                "max_epochs": 1,
                "batch_size": 2,
                "early_stopping": False,
                "train_size": 0.75,
                "validation_size": None,
                "shuffle_set_split": True,
                "load_sparse_tensor": False,
            }
            out_dir = root / "candidate"
            config = {
                "cfg_path": str(root / "config.yaml"),
                "project": {"out_dir": str(out_dir)},
                "run": {"seed": 42},
                "scvi": {
                    "layer": None,
                    "categorical_covariates": ["batch"],
                    "continuous_covariates": [],
                    "allow_missing_covariates": False,
                    "accelerator": "cpu",
                    "devices": 1,
                    "model_args": model_args,
                    "training_args": training_args,
                },
            }
            config_path = root / "config.yaml"
            config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
            resolved_config = load_config(config_path)

            candidate = source.copy()
            candidate.obsm["X_scVI"] = np.ones((4, 2), dtype=np.float32)
            candidate.uns["scvi_training_provenance"] = {
                "started_at_utc": "2026-08-16T00:00:00+00:00",
                "completed_at_utc": "2026-08-16T00:01:00+00:00",
                "input_path": str(source_path.resolve()),
                "input_shape": [4, 3],
                "input_obs_names_sha256": _index_sha256(source.obs_names),
                "input_var_names_sha256": _index_sha256(source.var_names),
                "config_path": str(config_path.resolve()),
                "resolved_config_sha256": _resolved_config_sha256(resolved_config),
                "seed": 42,
                "layer": "X",
                "categorical_covariates": ["batch"],
                "continuous_covariates": [],
                "allow_missing_covariates": False,
                "model_args": model_args,
                "training_args": training_args,
                "accelerator": "cpu",
                "devices": 1,
                "scvi_tools_version": "fixture",
                "torch_version": "fixture",
                "anndata_version": "fixture",
                "scanpy_version": "fixture",
            }
            checkpoint_path = out_dir / "04_scvi.h5ad"
            checkpoint_path.parent.mkdir(parents=True)
            candidate.write_h5ad(checkpoint_path)
            model_dir = out_dir / "models" / "scvi_model"
            model_dir.mkdir(parents=True)
            (model_dir / "model.pt").write_bytes(b"fixture model")

            report = validate_scvi_checkpoint(
                config_path=config_path,
                input_path=source_path,
                checkpoint_path=checkpoint_path,
                report_path=out_dir / "validation" / "04_scvi.json",
            )

            self.assertEqual(report["status"], "passed")
            self.assertFalse(report["failed_checks"])

            reordered = candidate[["c2", "c1", "c3", "c4"]].copy()
            reordered.write_h5ad(checkpoint_path)
            failed_report = out_dir / "validation" / "04_scvi_reordered.json"
            with self.assertRaisesRegex(ValueError, "cell_identity_and_order_preserved"):
                validate_scvi_checkpoint(
                    config_path=config_path,
                    input_path=source_path,
                    checkpoint_path=checkpoint_path,
                    report_path=failed_report,
                )
            self.assertEqual(
                yaml.safe_load(failed_report.read_text(encoding="utf-8"))["status"],
                "failed",
            )

    def test_composition_restores_structural_zeros(self) -> None:
        obs = pd.DataFrame(
            {
                "donor_id": ["d1", "d1", "d2", "d2", "d2"],
                "age": [30, 30, 50, 50, 50],
                "cell_type": ["A", "A", "A", "B", "B"],
            }
        )

        result = _build_donor_fraction_table(obs, min_cells_per_donor=1, covariate_cols=[])

        self.assertEqual(len(result), 4)
        missing_population = result[result["donor_id"].eq("d1") & result["cell_type"].eq("B")].iloc[
            0
        ]
        self.assertEqual(int(missing_population["n_cells"]), 0)
        self.assertEqual(float(missing_population["fraction"]), 0.0)
        sums = result.groupby("donor_id")["fraction"].sum()
        np.testing.assert_allclose(sums.to_numpy(), np.ones(len(sums)))

    def test_composition_denominator_retains_nonprimary_cells_without_testing_labels(self) -> None:
        obs = pd.DataFrame(
            {
                "cell_type": ["CD4 T cells", "NKT cells", "Fibroblasts"],
                "cell_type_analysis": ["CD4 T cells", pd.NA, pd.NA],
                "cell_type_analysis_disposition": ["primary", "exploratory", "unresolved"],
            }
        )
        config = {
            "annotation_qualification": {
                "enabled": True,
                "disposition_col": "cell_type_analysis_disposition",
            },
            "composition_age": {
                "denominator": ALL_QC_DENOMINATOR,
                "other_unresolved_label": "Other/unresolved",
            },
        }

        labels, metadata = _resolve_composition_cell_types(
            obs,
            config,
            celltype_col="cell_type_analysis",
        )

        self.assertEqual(labels.tolist(), ["CD4 T cells", "Other/unresolved", "Other/unresolved"])
        self.assertEqual(obs["cell_type"].tolist(), ["CD4 T cells", "NKT cells", "Fibroblasts"])
        self.assertEqual(metadata["nonprimary_cells_in_denominator"], 2)
        fractions = _build_donor_fraction_table(
            pd.DataFrame(
                {
                    "donor_id": ["sample_1"] * 3,
                    "age": [40] * 3,
                    "cell_type": labels,
                }
            ),
            min_cells_per_donor=1,
            covariate_cols=[],
        )
        self.assertEqual(int(fractions["total_cells"].iloc[0]), 3)
        observed = dict(zip(fractions["cell_type"], fractions["fraction"]))
        self.assertAlmostEqual(observed["CD4 T cells"], 1 / 3)
        self.assertAlmostEqual(observed["Other/unresolved"], 2 / 3)

    def test_complete_composition_denominator_requires_dispositions(self) -> None:
        obs = pd.DataFrame({"cell_type_analysis": ["CD4 T cells", pd.NA]})
        config = {
            "annotation_qualification": {"enabled": True},
            "composition_age": {"denominator": ALL_QC_DENOMINATOR},
        }

        with self.assertRaisesRegex(KeyError, "disposition column"):
            _resolve_composition_cell_types(
                obs,
                config,
                celltype_col="cell_type_analysis",
            )

    def test_trends_separate_tested_and_detected_replicates(self) -> None:
        fractions = pd.DataFrame(
            {
                "donor_id": [f"d{i}" for i in range(6)],
                "age": [20, 30, 40, 50, 60, 70],
                "cell_type": ["rare"] * 6,
                "n_cells": [0, 0, 1, 0, 2, 3],
                "total_cells": [10] * 6,
                "fraction": [0.0, 0.0, 0.1, 0.0, 0.2, 0.3],
            }
        )

        trends = _compute_trends(
            fractions,
            min_donors_per_celltype=3,
            covariate_cols=[],
            adjust_covariates=False,
            bootstrap_iterations=0,
            bootstrap_ci=0.95,
            seed=42,
        )

        self.assertEqual(int(trends.iloc[0]["n_donors"]), 6)
        self.assertEqual(int(trends.iloc[0]["n_donors_detected"]), 3)

    def test_longitudinal_composition_support_is_counted_by_subject(self) -> None:
        fractions = pd.DataFrame(
            {
                "donor_id": [f"sample_{index}" for index in range(6)],
                "subject_id": np.repeat(["subject_1", "subject_2", "subject_3"], 2),
                "age": [30, 31, 40, 41, 50, 51],
                "cell_type": ["T"] * 6,
                "n_cells": [2, 3, 4, 5, 6, 7],
                "total_cells": [10] * 6,
                "fraction": [0.2, 0.3, 0.4, 0.5, 0.6, 0.7],
            }
        )

        trends = _compute_trends(
            fractions,
            min_donors_per_celltype=4,
            covariate_cols=[],
            adjust_covariates=False,
            bootstrap_iterations=0,
            bootstrap_ci=0.95,
            seed=42,
            subject_aware=True,
        )

        self.assertTrue(trends.empty)

    def test_subject_aware_composition_requires_explicit_subject_id(self) -> None:
        obs = pd.DataFrame(
            {
                "donor_id": ["sample_1", "sample_2"],
                "age": [30, 50],
                "cell_type": ["A", "A"],
            }
        )

        with self.assertRaisesRegex(KeyError, "requires an explicit 'subject_id'"):
            _build_donor_fraction_table(
                obs,
                min_cells_per_donor=1,
                covariate_cols=[],
                require_subject_id=True,
            )

    def test_longitudinal_analyses_do_not_fallback_to_sample_as_subject(self) -> None:
        with self._temporary_directory() as directory:
            checkpoint = Path(directory) / "checkpoint.h5ad"
            ad.AnnData(
                X=np.ones((2, 1)),
                obs=pd.DataFrame(
                    {
                        "age": [30, 50],
                        "sample_unit_id": ["sample_1", "sample_2"],
                        "cell_type": ["T", "T"],
                    },
                    index=["cell_1", "cell_2"],
                ),
            ).write_h5ad(checkpoint)
            config = {
                "composition_age": {
                    "age_col": "age",
                    "sample_unit_col": "sample_unit_id",
                    "subject_col": "subject_id",
                    "celltype_col": "cell_type",
                    "covariate_cols": [],
                }
            }

            prepared, columns = _prepare_obs(str(checkpoint), config)

        self.assertTrue(prepared.empty)
        self.assertIsNone(columns["subject_col"])
        signature_columns = _pick_columns(
            ["age", "sample_unit_id", "cell_type"],
            {
                "signature_age": {
                    "age_col": "age",
                    "sample_unit_col": "sample_unit_id",
                    "subject_col": "subject_id",
                    "celltype_col": "cell_type",
                }
            },
        )
        self.assertIsNone(signature_columns["subject"])

    def test_subject_aware_gee_recovers_positive_age_direction(self) -> None:
        subject_ids = np.repeat([f"subject_{index}" for index in range(8)], 2)
        age = np.repeat(np.arange(25, 65, 5, dtype=float), 2) + np.tile([0.0, 1.0], 8)
        subject_offset = np.repeat(np.linspace(-0.15, 0.15, 8), 2)
        fraction = 1.0 / (1.0 + np.exp(-(-2.0 + 0.03 * age + subject_offset)))
        frame = pd.DataFrame(
            {
                "subject_id": subject_ids,
                "age": age,
                "fraction": fraction,
                "total_cells": 100,
            }
        )

        result = fit_subject_aware_gee(
            frame,
            outcome_col="fraction",
            subject_col="subject_id",
            covariate_cols=[],
            family="binomial",
            weights_col="total_cells",
            context="bounded composition fixture",
        )

        self.assertTrue(result["converged"])
        self.assertEqual(result["n_subjects"], 8)
        self.assertGreater(result["effect_per_10y"], 0.0)
        self.assertTrue(np.isfinite(result["pvalue"]))
        self.assertEqual(result["working_correlation_structure"], "exchangeable")
        self.assertEqual(result["covariance_type"], "robust_subject_clustered")

    def test_subject_aware_independence_gee_retains_clustered_covariance(self) -> None:
        subject_ids = np.repeat([f"subject_{index}" for index in range(12)], 2)
        age = np.repeat(np.arange(25, 85, 5, dtype=float), 2) + np.tile([0.0, 1.0], 12)
        subject_offset = np.repeat(np.linspace(-0.3, 0.3, 12), 2)
        score = 0.02 * age + subject_offset
        frame = pd.DataFrame({"subject_id": subject_ids, "age": age, "score": score})

        result = fit_subject_aware_gee(
            frame,
            outcome_col="score",
            subject_col="subject_id",
            covariate_cols=[],
            family="gaussian",
            working_correlation="independence",
            context="bounded signature fallback fixture",
        )

        self.assertTrue(result["converged"])
        self.assertEqual(result["n_subjects"], 12)
        self.assertEqual(result["working_correlation_structure"], "independence")
        self.assertEqual(result["covariance_type"], "robust_subject_clustered")
        self.assertTrue(np.isnan(result["working_correlation"]))
        self.assertGreater(result["effect_per_10y"], 0.0)
        self.assertTrue(np.isfinite(result["standard_error"]))

    def test_signature_gee_uses_fallback_only_after_specific_nonconvergence(self) -> None:
        frame = pd.DataFrame(
            {
                "subject_id": ["S1", "S1", "S2", "S2", "S3", "S3"],
                "age": [30, 31, 40, 41, 50, 51],
                "score": [0.1, 0.2, 0.2, 0.3, 0.4, 0.5],
            }
        )
        fallback_result = {
            "model": "GEE-gaussian-independence",
            "working_correlation_structure": "independence",
            "working_correlation": np.nan,
            "covariance_type": "robust_subject_clustered",
        }

        with patch(
            "src.signature_age.fit_subject_aware_gee",
            side_effect=[GEENonConvergenceError("primary failed"), fallback_result],
        ) as fit:
            result = _fit_signature_gee(
                frame,
                outcome_col="score",
                subject_col="subject_id",
                covariate_cols=[],
                primary_working_correlation="exchangeable",
                nonconvergence_fallback="independence",
                context="bounded fallback trigger fixture",
            )

        self.assertEqual(fit.call_count, 2)
        self.assertIs(fit.call_args_list[0].args[0], frame)
        self.assertIs(fit.call_args_list[1].args[0], frame)
        self.assertEqual(fit.call_args_list[0].kwargs["working_correlation"], "exchangeable")
        self.assertEqual(fit.call_args_list[1].kwargs["working_correlation"], "independence")
        self.assertTrue(result["fallback_used"])
        self.assertEqual(result["fallback_reason"], "exchangeable_gee_nonconvergence")
        self.assertTrue(result["fallback_requires_manual_review"])

    def test_signature_gee_does_not_mask_nonconvergence_unrelated_errors(self) -> None:
        frame = pd.DataFrame(
            {"subject_id": ["S1", "S2", "S3"], "age": [30, 40, 50], "score": [1, 2, 3]}
        )
        with patch(
            "src.signature_age.fit_subject_aware_gee",
            side_effect=RuntimeError("non-finite inference"),
        ) as fit:
            with self.assertRaisesRegex(RuntimeError, "non-finite inference"):
                _fit_signature_gee(
                    frame,
                    outcome_col="score",
                    subject_col="subject_id",
                    covariate_cols=[],
                    primary_working_correlation="exchangeable",
                    nonconvergence_fallback="independence",
                    context="bounded fallback negative fixture",
                )

        self.assertEqual(fit.call_count, 1)

    def test_signature_gee_stops_when_prespecified_fallback_does_not_converge(self) -> None:
        frame = pd.DataFrame(
            {"subject_id": ["S1", "S2", "S3"], "age": [30, 40, 50], "score": [1, 2, 3]}
        )
        with patch(
            "src.signature_age.fit_subject_aware_gee",
            side_effect=[
                GEENonConvergenceError("primary failed"),
                GEENonConvergenceError("fallback failed"),
            ],
        ) as fit:
            with self.assertRaisesRegex(GEENonConvergenceError, "fallback failed"):
                _fit_signature_gee(
                    frame,
                    outcome_col="score",
                    subject_col="subject_id",
                    covariate_cols=[],
                    primary_working_correlation="exchangeable",
                    nonconvergence_fallback="independence",
                    context="bounded double nonconvergence fixture",
                )

        self.assertEqual(fit.call_count, 2)

    def test_one_sample_sensitivity_is_deterministic_and_estimable(self) -> None:
        metadata = pd.DataFrame(
            {
                "subject_id": ["S1", "S1", "S2", "S2", "S3", "S4"],
                "sample_unit_id": ["B", "A", "D", "C", "E", "F"],
                "age": [30, 30, 51, 50, 60, 70],
            }
        )
        selected = select_one_sample_per_subject(
            metadata,
            subject_col="subject_id",
            sample_col="sample_unit_id",
        )
        self.assertEqual(selected, {"A", "C", "E", "F"})

        frame = pd.DataFrame(
            {
                "age": [30, 40, 50, 60, 70, 80],
                "score": [0.1, 0.2, 0.25, 0.4, 0.5, 0.55],
            }
        )
        result = fit_independent_age_model(
            frame,
            outcome_col="score",
            covariate_cols=[],
            family="gaussian",
            context="bounded one-sample fixture",
        )
        self.assertTrue(result["converged"])
        self.assertGreater(result["effect_per_10y"], 0.0)

    def test_missing_explicit_covariate_fails(self) -> None:
        with self.assertRaisesRegex(KeyError, "covariates are missing"):
            resolve_covariates(
                ["age", "sex"],
                {"covariate_cols": ["sex", "batch"]},
                section_name="composition_age",
            )

    def test_missing_covariate_values_are_not_imputed(self) -> None:
        frame = pd.DataFrame({"sex": ["F", None], "batch": ["A", "B"]})
        with self.assertRaisesRegex(ValueError, "missing values"):
            build_complete_nuisance_design(frame, context="test model")

    def test_rank_deficient_design_fails(self) -> None:
        frame = pd.DataFrame({"covariate_a": [1, 2, 3], "covariate_b": [2, 4, 6]})
        with self.assertRaisesRegex(ValueError, "rank deficient"):
            build_complete_nuisance_design(frame, context="test model")

    def test_numeric_looking_category_is_not_treated_as_continuous(self) -> None:
        frame = pd.DataFrame({"batch": ["1", "2", "3"]}, dtype="string")
        design = build_complete_nuisance_design(frame, context="test model")
        self.assertEqual(design.shape, (3, 3))

    def test_constant_configured_covariate_fails(self) -> None:
        frame = pd.DataFrame({"sex": ["F", "F", "F"]})
        with self.assertRaisesRegex(ValueError, "has no variation"):
            build_complete_nuisance_design(frame, context="test model")

    def test_inconsistent_replicate_covariate_fails(self) -> None:
        frame = pd.DataFrame({"replicate": ["r1", "r1"], "batch": ["A", "B"]})
        with self.assertRaisesRegex(ValueError, "inconsistent"):
            validate_replicate_covariates(
                frame,
                replicate_col="replicate",
                covariate_cols=["batch"],
                context="test model",
            )

    def test_real_profile_cannot_write_placeholder(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "disabled"):
            require_placeholder_permission(
                {"run": {"allow_placeholder_outputs": False}},
                "Analysis failed.",
            )

    def test_demo_profile_can_write_placeholder(self) -> None:
        require_placeholder_permission(
            {"run": {"allow_placeholder_outputs": True}},
            "Analysis failed.",
        )


if __name__ == "__main__":
    unittest.main()
