from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from src.biological_replicates import (
    _build_mapping,
    attach_biological_replicates,
)
from src.signature_age import _materialize_signature_genes


def _config(table: Path, mapping: Path | None = None) -> dict:
    return {
        "biological_replicates": {
            "canonical_col": "biological_replicate_id",
            "source_col": "tube_id",
            "source_table_path": str(table),
            "source_table_join_key": "cell_id",
            "mapping_path": str(mapping) if mapping else "",
            "strict": True,
        }
    }


def _source_table(path: Path) -> None:
    pd.DataFrame(
        {
            "cell_id": ["cell_a", "cell_b", "cell_c"],
            "Tube_id": ["tube_1", "tube_1", "tube_2"],
            "Donor_id": ["label_1", "label_1", "label_2"],
            "Age": [30, 30, 60],
            "Sex": ["Female", "Female", "Male"],
            "Batch": ["batch_1", "batch_1", "batch_2"],
            "File_name": ["lib_1", "lib_1", "lib_2"],
        }
    ).to_csv(path, index=False)


class BiologicalReplicateTests(unittest.TestCase):
    @staticmethod
    def _temporary_directory() -> tempfile.TemporaryDirectory:
        root = Path(".tmp") / "tests"
        root.mkdir(parents=True, exist_ok=True)
        return tempfile.TemporaryDirectory(dir=root)

    def test_build_mapping_is_keyed_and_conflict_free(self) -> None:
        with self._temporary_directory() as tmp:
            table = Path(tmp) / "metadata.csv"
            _source_table(table)
            obs = pd.DataFrame(
                {
                    "donor_id": ["label_2", "label_1", "label_1"],
                    "age": [60, 30, 30],
                    "sex": ["Male", "Female", "Female"],
                    "batch": ["batch_2", "batch_1", "batch_1"],
                    "sample_id": ["lib_2", "lib_1", "lib_1"],
                },
                index=["cell_c", "cell_a", "cell_b"],
            )

            mapping, replicate_metadata, report = _build_mapping(
                obs,
                _config(table),
            )

            self.assertEqual(
                mapping["biological_replicate_id"].tolist(),
                ["tube_2", "tube_1", "tube_1"],
            )
            self.assertEqual(replicate_metadata.shape[0], 2)
            self.assertTrue(report["passed"])
            self.assertEqual(report["unique_biological_replicates"], 2)

    def test_attach_mapping_rejects_missing_cells(self) -> None:
        with self._temporary_directory() as tmp:
            tmp_path = Path(tmp)
            mapping = tmp_path / "mapping.csv"
            pd.DataFrame(
                {
                    "cell_id": ["cell_a"],
                    "biological_replicate_id": ["tube_1"],
                }
            ).to_csv(mapping, index=False)
            adata = ad.AnnData(
                X=np.ones((2, 1)),
                obs=pd.DataFrame(index=["cell_a", "cell_b"]),
            )

            with self.assertRaisesRegex(ValueError, "missed 1"):
                attach_biological_replicates(
                    adata,
                    _config(tmp_path / "unused.csv", mapping),
                )

    def test_build_mapping_rejects_replicate_age_conflict(self) -> None:
        with self._temporary_directory() as tmp:
            table = Path(tmp) / "metadata.csv"
            _source_table(table)
            source = pd.read_csv(table)
            source.loc[source["cell_id"].eq("cell_b"), "Age"] = 31
            source.to_csv(table, index=False)
            obs = pd.DataFrame(index=["cell_a", "cell_b", "cell_c"])

            _, _, report = _build_mapping(obs, _config(table))

            self.assertFalse(report["passed"])
            self.assertEqual(
                report["replicate_conflict_counts"]["age_nunique"],
                1,
            )

    def test_signature_gene_materialization_preserves_normalization(
        self,
    ) -> None:
        matrix = np.array(
            [
                [1, 3, 6],
                [0, 4, 6],
                [5, 0, 5],
            ],
            dtype=np.float32,
        )
        adata = ad.AnnData(
            X=matrix,
            obs=pd.DataFrame(index=["a", "b", "c"]),
            var=pd.DataFrame(index=["G1", "G2", "G3"]),
        )
        observed = _materialize_signature_genes(
            adata,
            selected_indices=np.array([2, 0]),
            signatures={"test": ["G1", "G3"]},
            assume_log1p=False,
            target_sum=100.0,
            chunk_size=1,
        )
        expected = np.log1p(
            matrix[[2, 0]][:, [0, 2]] * (100.0 / matrix[[2, 0]].sum(axis=1))[:, None]
        )

        np.testing.assert_allclose(
            observed.X.toarray(),
            expected,
            rtol=1e-6,
        )
        self.assertEqual(observed.obs_names.tolist(), ["c", "a"])


if __name__ == "__main__":
    unittest.main()
