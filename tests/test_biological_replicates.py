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
        "biological_units": {
            "subject_col": "subject_id",
            "subject_source_col": "donor_id",
            "sample_unit_col": "sample_unit_id",
            "sample_unit_source_col": "tube_id",
            "technical_library_col": "technical_library_id",
            "technical_library_source_col": "file_name",
            "source_table_path": str(table),
            "source_table_join_key": "cell_id",
            "mapping_path": str(mapping) if mapping else "",
            "expected_subjects": 2,
            "expected_sample_units": 3,
            "strict": True,
        }
    }


def _source_table(path: Path) -> None:
    pd.DataFrame(
        {
            "cell_id": ["cell_a", "cell_b", "cell_c", "cell_d"],
            "Tube_id": ["tube_1", "tube_1", "tube_2", "tube_3"],
            "Donor_id": ["donor_1", "donor_1", "donor_1", "donor_2"],
            "Age": [30, 30, 31, 60],
            "Sex": ["Female", "Female", "Female", "Male"],
            "Batch": ["batch_1", "batch_1", "batch_2", "batch_3"],
            "File_name": ["lib_1", "lib_1", "lib_2", "lib_3"],
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
                    "donor_id": ["donor_2", "donor_1", "donor_1", "donor_1"],
                    "age": [60, 31, 30, 30],
                    "sex": ["Male", "Female", "Female", "Female"],
                    "batch": ["batch_3", "batch_2", "batch_1", "batch_1"],
                    "sample_id": ["lib_3", "lib_2", "lib_1", "lib_1"],
                },
                index=["cell_d", "cell_c", "cell_a", "cell_b"],
            )

            mapping, replicate_metadata, report = _build_mapping(
                obs,
                _config(table),
            )

            self.assertEqual(
                mapping["sample_unit_id"].tolist(),
                ["tube_3", "tube_2", "tube_1", "tube_1"],
            )
            self.assertEqual(mapping["subject_id"].nunique(), 2)
            self.assertEqual(replicate_metadata.shape[0], 3)
            self.assertTrue(report["passed"])
            self.assertEqual(report["unique_subjects"], 2)
            self.assertEqual(report["unique_sample_units"], 3)
            self.assertEqual(report["longitudinal_design"]["subjects_with_repeated_samples"], 1)
            self.assertEqual(report["longitudinal_design"]["subjects_with_age_variation"], 1)

    def test_attach_mapping_rejects_missing_cells(self) -> None:
        with self._temporary_directory() as tmp:
            tmp_path = Path(tmp)
            mapping = tmp_path / "mapping.csv"
            pd.DataFrame(
                {
                    "cell_id": ["cell_a"],
                    "subject_id": ["donor_1"],
                    "sample_unit_id": ["tube_1"],
                    "technical_library_id": ["lib_1"],
                }
            ).to_csv(mapping, index=False)
            adata = ad.AnnData(
                X=np.ones((2, 1)),
                obs=pd.DataFrame(index=["cell_a", "cell_b"]),
            )

            with self.assertRaisesRegex(ValueError, "sample_unit_id=1"):
                attach_biological_replicates(
                    adata,
                    _config(tmp_path / "unused.csv", mapping),
                )

    def test_build_mapping_rejects_replicate_age_conflict(self) -> None:
        with self._temporary_directory() as tmp:
            table = Path(tmp) / "metadata.csv"
            _source_table(table)
            source = pd.read_csv(table)
            source.loc[source["cell_id"].eq("cell_b"), "Age"] = 32
            source.to_csv(table, index=False)
            obs = pd.DataFrame(index=["cell_a", "cell_b", "cell_c", "cell_d"])

            _, _, report = _build_mapping(obs, _config(table))

            self.assertFalse(report["passed"])
            self.assertEqual(
                report["sample_unit_conflict_counts"]["age_nunique"],
                1,
            )

    def test_build_mapping_rejects_subject_sex_conflict_across_samples(self) -> None:
        with self._temporary_directory() as tmp:
            table = Path(tmp) / "metadata.csv"
            _source_table(table)
            source = pd.read_csv(table)
            source.loc[source["cell_id"].eq("cell_c"), "Sex"] = "Male"
            source.to_csv(table, index=False)
            obs = pd.DataFrame(index=["cell_a", "cell_b", "cell_c", "cell_d"])

            _, _, report = _build_mapping(obs, _config(table))

            self.assertFalse(report["passed"])
            self.assertEqual(report["subject_conflict_counts"]["sex_nunique"], 1)

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
