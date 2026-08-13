from __future__ import annotations

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import anndata as ad
import numpy as np
import pandas as pd
import yaml

from src.audit_scientific_checkpoint import audit
from src.composition_age import (
    _build_donor_fraction_table,
    _compute_trends,
)
from src.scientific_guardrails import (
    build_complete_nuisance_design,
    require_placeholder_permission,
    resolve_covariates,
    validate_replicate_covariates,
)


class ScientificGuardrailTests(unittest.TestCase):
    @staticmethod
    def _temporary_directory() -> TemporaryDirectory:
        root = Path(".tmp/tests")
        root.mkdir(parents=True, exist_ok=True)
        return TemporaryDirectory(dir=root)

    @staticmethod
    def _checkpoint_fixture(root: Path, *, complete_provenance: bool) -> tuple[Path, Path, Path]:
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
            "seed": 42,
            "use_rep": "X_scVI",
            "leiden_resolution": 0.8,
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
                "seed": 42,
                "layer": "X",
                "categorical_covariates": ["batch"],
                "continuous_covariates": [],
                "allow_missing_covariates": False,
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
