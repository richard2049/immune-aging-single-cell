from __future__ import annotations

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import anndata as ad
import numpy as np
import pandas as pd
import yaml
from scipy import sparse

from src.cluster import _h5ad_safe_arguments
from src.scvi_train import _index_sha256
from src.utils import load_config
from src.validate_cluster_checkpoint import validate


class ClusterCheckpointTests(unittest.TestCase):
    @staticmethod
    def _temporary_directory() -> TemporaryDirectory:
        root = Path(".tmp/tests")
        root.mkdir(parents=True, exist_ok=True)
        return TemporaryDirectory(dir=root)

    @staticmethod
    def _fixture(root: Path) -> tuple[Path, Path, Path]:
        config = load_config("config/demo.yaml")
        config["project"]["out_dir"] = str(root / "results")
        config_path = root / "config.yaml"
        config_path.write_text(yaml.safe_dump(config), encoding="utf-8")

        obs = pd.DataFrame(index=[f"cell_{index}" for index in range(6)])
        var = pd.DataFrame(index=[f"gene_{index}" for index in range(3)])
        source = ad.AnnData(X=sparse.csr_matrix(np.ones((6, 3))), obs=obs, var=var)
        source.obsm["X_scVI"] = np.arange(12, dtype=np.float32).reshape(6, 2)
        source_path = root / "04_scvi.h5ad"
        source.write_h5ad(source_path)

        candidate = source.copy()
        candidate.obsm["X_umap"] = np.arange(12, dtype=np.float32).reshape(6, 2)
        candidate.obsp["distances"] = sparse.eye(6, format="csr")
        candidate.obsp["connectivities"] = sparse.eye(6, format="csr")
        candidate.obs["leiden"] = pd.Categorical(["0", "0", "0", "1", "1", "1"])
        clustering = config["clustering"]
        candidate.uns["clustering_provenance"] = {
            "started_at_utc": "2026-08-20T00:00:00+00:00",
            "completed_at_utc": "2026-08-20T00:01:00+00:00",
            "input_path": str(source_path.resolve()),
            "input_shape": [6, 3],
            "input_obs_names_sha256": _index_sha256(source.obs_names),
            "input_var_names_sha256": _index_sha256(source.var_names),
            "config_path": str(config_path.resolve()),
            "resolved_config_sha256": "fixture-config",
            "seed": int(config["run"]["seed"]),
            "n_jobs": int(clustering["n_jobs"]),
            "use_rep": clustering["use_rep"],
            "leiden_resolution": float(clustering["leiden_resolution"]),
            "neighbors_args": _h5ad_safe_arguments(clustering["neighbors_args"]),
            "umap_args": _h5ad_safe_arguments(clustering["umap_args"]),
            "leiden_args": _h5ad_safe_arguments(clustering["leiden_args"]),
            "scanpy_version": "fixture",
        }
        checkpoint_path = root / "05_clustered.h5ad"
        candidate.write_h5ad(checkpoint_path)
        return config_path, source_path, checkpoint_path

    def test_cluster_gate_accepts_complete_checkpoint(self) -> None:
        with self._temporary_directory() as directory:
            root = Path(directory)
            config_path, source_path, checkpoint_path = self._fixture(root)
            report = validate(
                config_path=config_path,
                input_path=source_path,
                checkpoint_path=checkpoint_path,
                report_path=root / "validation.json",
                historical_path=checkpoint_path,
            )

        self.assertEqual(report["status"], "passed")
        self.assertFalse(report["failed_checks"])
        self.assertEqual(report["summary"]["candidate_clusters"], 2)

    def test_cluster_gate_rejects_reordered_cells(self) -> None:
        with self._temporary_directory() as directory:
            root = Path(directory)
            config_path, source_path, checkpoint_path = self._fixture(root)
            candidate = ad.read_h5ad(checkpoint_path)
            candidate = candidate[
                ["cell_1", "cell_0", "cell_2", "cell_3", "cell_4", "cell_5"]
            ].copy()
            candidate.write_h5ad(checkpoint_path)

            with self.assertRaisesRegex(ValueError, "cell_identity_and_order_preserved"):
                validate(
                    config_path=config_path,
                    input_path=source_path,
                    checkpoint_path=checkpoint_path,
                    report_path=root / "validation.json",
                )


if __name__ == "__main__":
    unittest.main()
