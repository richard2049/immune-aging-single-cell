from __future__ import annotations

import builtins
import hashlib
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd
import yaml

from src.age_prediction import _load_model_builders
from src.annotate_celltypist import _load_model_with_provenance
from src.build_donor_metadata import _aggregate_replicates, _is_raw_output
from src.metadata_integrate import _merge_external_metadata
from src.scvi_train import _resolve_covariates


class ReproducibilityContractTests(unittest.TestCase):
    @staticmethod
    def _temporary_directory() -> tempfile.TemporaryDirectory:
        root = Path(".tmp") / "tests"
        root.mkdir(parents=True, exist_ok=True)
        return tempfile.TemporaryDirectory(dir=root)

    def test_metadata_join_preserves_order_and_accepts_exact_duplicates(self) -> None:
        with self._temporary_directory() as tmp:
            table = Path(tmp) / "metadata.csv"
            pd.DataFrame(
                {
                    "cell_id": ["cell_b", "cell_a", "cell_a"],
                    "age": [60, 30, 30],
                }
            ).to_csv(table, index=False)
            adata = ad.AnnData(
                X=np.ones((2, 1)),
                obs=pd.DataFrame(index=["cell_a", "cell_b"]),
            )

            _merge_external_metadata(
                adata,
                {
                    "table_path": str(table),
                    "obs_join_key": "cell_id",
                    "table_join_key": "cell_id",
                },
            )

            self.assertEqual(adata.obs_names.tolist(), ["cell_a", "cell_b"])
            self.assertEqual(adata.obs["age"].tolist(), [30, 60])
            self.assertTrue(adata.obs["metadata_matched"].all())

    def test_metadata_join_rejects_conflicting_duplicate_keys(self) -> None:
        with self._temporary_directory() as tmp:
            table = Path(tmp) / "metadata.csv"
            pd.DataFrame({"cell_id": ["cell_a", "cell_a"], "age": [30, 31]}).to_csv(
                table, index=False
            )
            adata = ad.AnnData(
                X=np.ones((1, 1)),
                obs=pd.DataFrame(index=["cell_a"]),
            )

            with self.assertRaisesRegex(ValueError, "conflicting duplicate"):
                _merge_external_metadata(
                    adata,
                    {
                        "table_path": str(table),
                        "obs_join_key": "cell_id",
                        "table_join_key": "cell_id",
                    },
                )

    def test_metadata_join_rejects_unmatched_cells_by_default(self) -> None:
        with self._temporary_directory() as tmp:
            table = Path(tmp) / "metadata.csv"
            pd.DataFrame({"cell_id": ["cell_a"], "age": [30]}).to_csv(table, index=False)
            adata = ad.AnnData(
                X=np.ones((2, 1)),
                obs=pd.DataFrame(index=["cell_a", "cell_b"]),
            )

            with self.assertRaisesRegex(ValueError, "observations unmatched"):
                _merge_external_metadata(
                    adata,
                    {
                        "table_path": str(table),
                        "obs_join_key": "cell_id",
                        "table_join_key": "cell_id",
                    },
                )

    def test_replicate_aggregation_does_not_merge_reused_donor_labels(self) -> None:
        cells = pd.DataFrame(
            {
                "biological_replicate_id": ["tube_1", "tube_1", "tube_2"],
                "donor_id": ["label_1", "label_1", "label_1"],
                "sample_id": ["lib_1", "lib_1", "lib_2"],
                "batch": ["batch_1", "batch_1", "batch_2"],
                "age": [30, 30, 60],
            }
        )

        observed = _aggregate_replicates(cells)

        self.assertEqual(observed["biological_replicate_id"].tolist(), ["tube_1", "tube_2"])
        self.assertEqual(observed["n_cells"].tolist(), [2, 1])

    def test_replicate_aggregation_rejects_age_conflicts(self) -> None:
        cells = pd.DataFrame(
            {
                "biological_replicate_id": ["tube_1", "tube_1"],
                "sample_id": ["lib_1", "lib_1"],
                "batch": ["batch_1", "batch_1"],
                "age": [30, 31],
            }
        )

        with self.assertRaisesRegex(ValueError, "conflicting age"):
            _aggregate_replicates(cells)

    def test_derived_metadata_cannot_be_written_under_raw(self) -> None:
        self.assertTrue(_is_raw_output(Path("data/raw/donor_metadata.csv")))
        self.assertFalse(_is_raw_output(Path("data/derived/donor_metadata.csv")))

    def test_configured_prediction_dependencies_are_declared(self) -> None:
        environment = Path("environment.yml").read_text(encoding="utf-8").lower()
        for path in (
            Path("config/gse164378_pilot.yaml"),
            Path("config/gse164378_1m.yaml"),
        ):
            config = yaml.safe_load(path.read_text(encoding="utf-8"))
            requested = config["age_prediction"]["model_order"]
            if "xgboost" in requested:
                self.assertIn("xgboost=", environment)

    def test_requested_prediction_dependency_cannot_be_silently_skipped(self) -> None:
        original_import = builtins.__import__

        def block_xgboost(name, *args, **kwargs):
            if name == "xgboost":
                raise ImportError("bounded missing-dependency fixture")
            return original_import(name, *args, **kwargs)

        config = {"age_prediction": {"model_order": ["xgboost"]}}
        with patch("builtins.__import__", side_effect=block_xgboost):
            with self.assertRaisesRegex(ImportError, "requests 'xgboost'"):
                _load_model_builders(config, latent_cols=["z0"], seed=42)

    def test_celltypist_provenance_records_and_checks_model_hash(self) -> None:
        with self._temporary_directory() as tmp:
            model_path = Path(tmp) / "model.pkl"
            model_path.write_bytes(b"bounded-celltypist-model")
            expected_hash = hashlib.sha256(model_path.read_bytes()).hexdigest()
            loaded_model = SimpleNamespace(
                description={"date": "2024-01-01", "version": "fixture-v1"}
            )
            fake_models = SimpleNamespace(
                Model=SimpleNamespace(load=lambda model: loaded_model),
                get_model_path=lambda model: str(model_path),
            )
            fake_celltypist = SimpleNamespace(models=fake_models)

            with patch.dict(sys.modules, {"celltypist": fake_celltypist}):
                with patch("src.annotate_celltypist.version", return_value="1.7.1"):
                    model, provenance = _load_model_with_provenance(
                        {
                            "model": "fixture.pkl",
                            "model_sha256": expected_hash,
                            "model_source": "bounded fixture",
                        }
                    )
                    self.assertIs(model, loaded_model)
                    self.assertEqual(provenance["model_sha256"], expected_hash)
                    self.assertEqual(provenance["model_version"], "fixture-v1")

                    with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                        _load_model_with_provenance(
                            {
                                "model": "fixture.pkl",
                                "model_sha256": "0" * 64,
                            }
                        )

    def test_missing_scvi_covariate_fails_before_training(self) -> None:
        with self.assertRaisesRegex(KeyError, "Configured scVI covariates"):
            _resolve_covariates(
                ["batch"],
                {
                    "categorical_covariates": ["batch"],
                    "continuous_covariates": ["pct_counts_mt"],
                    "allow_missing_covariates": False,
                },
            )

    def test_explicit_alternative_scvi_model_can_allow_missing_covariate(self) -> None:
        categorical, continuous = _resolve_covariates(
            ["batch"],
            {
                "categorical_covariates": ["batch"],
                "continuous_covariates": ["pct_counts_mt"],
                "allow_missing_covariates": True,
            },
        )
        self.assertEqual(categorical, ["batch"])
        self.assertEqual(continuous, [])


if __name__ == "__main__":
    unittest.main()
