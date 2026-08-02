from __future__ import annotations

import gzip
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml
from scipy.io import mmread

from src.pseudobulk_aggregate import (
    _integer_sparse_chunk,
    _write_matrix_market_gzip_atomic,
    aggregate_sparse_counts,
    run,
)
from src.run_edger import build_command, run_command
from src.validate_pseudobulk_de import _bh_adjust, _resolve_manifest_path


class PseudobulkDifferentialExpressionTests(unittest.TestCase):
    @staticmethod
    def _temporary_directory() -> tempfile.TemporaryDirectory:
        root = Path(".tmp") / "tests"
        root.mkdir(parents=True, exist_ok=True)
        return tempfile.TemporaryDirectory(dir=root)

    def test_command_log_preserves_utf8_subprocess_output(self) -> None:
        with self._temporary_directory() as tmp:
            log_path = Path(tmp) / "command.log"
            run_command(
                [
                    sys.executable,
                    "-c",
                    ("import sys; sys.stdout.buffer.write('diagnóstico'.encode('utf-8'))"),
                ],
                repository=Path.cwd().resolve(),
                log_path=log_path,
                timeout_seconds=30,
            )

            log_text = log_path.read_text(encoding="utf-8")
            self.assertIn("diagnóstico", log_text)
            self.assertIn("return_code: 0", log_text)

    def test_sparse_aggregation_conserves_profile_counts(self) -> None:
        matrix = sp.csr_matrix(
            [
                [1, 0, 2],
                [0, 3, 1],
                [4, 0, 0],
                [1, 1, 1],
                [9, 9, 9],
            ],
            dtype=np.int64,
        )
        profile_codes = np.array([0, 0, 1, 1, -1], dtype=np.int64)

        observed, library_sizes, included_cells = aggregate_sparse_counts(
            matrix,
            profile_codes=profile_codes,
            n_profiles=2,
            chunk_size=2,
        )

        np.testing.assert_array_equal(
            observed.toarray(),
            np.array([[1, 3, 3], [5, 1, 1]], dtype=np.int64),
        )
        np.testing.assert_array_equal(library_sizes, [7, 7])
        self.assertEqual(included_cells, 4)

    def test_bh_adjustment_matches_known_values(self) -> None:
        observed = _bh_adjust(np.array([0.01, 0.04, 0.03, 0.002]))
        np.testing.assert_allclose(observed, [0.02, 0.04, 0.04, 0.008])

    def test_manifest_paths_must_be_relative_and_contained(self) -> None:
        with self._temporary_directory() as tmp:
            output_dir = Path(tmp).resolve()
            self.assertEqual(
                _resolve_manifest_path(output_dir, "by_cell_type/result.csv.gz"),
                output_dir / "by_cell_type" / "result.csv.gz",
            )
            with self.assertRaisesRegex(ValueError, "must be relative"):
                _resolve_manifest_path(output_dir, "/work/result.csv.gz")
            with self.assertRaisesRegex(ValueError, "escapes"):
                _resolve_manifest_path(output_dir, "../result.csv.gz")

    def test_count_validation_rejects_dense_and_fractional_values(self) -> None:
        with self.assertRaisesRegex(TypeError, "sparse count matrix"):
            _integer_sparse_chunk(np.ones((2, 2)))
        with self.assertRaisesRegex(ValueError, "fractional"):
            _integer_sparse_chunk(sp.csr_matrix([[1.0, 0.5]]))

    def test_bounded_aggregation_writes_auditable_matrix_market(self) -> None:
        with self._temporary_directory() as tmp:
            root = Path(tmp)
            obs = pd.DataFrame(
                {
                    "biological_replicate_id": [
                        "R1",
                        "R1",
                        "R2",
                        "R2",
                        "R1",
                        "R2",
                    ],
                    "cell_type": ["A", "A", "A", "A", "B", "B"],
                    "age": [30, 30, 50, 50, 30, 50],
                    "sex": ["F", "F", "M", "M", "F", "M"],
                    "batch": ["X", "X", "Y", "Y", "X", "Y"],
                    "sample_id": ["L1", "L1", "L2", "L2", "L1", "L2"],
                },
                index=[f"cell_{index}" for index in range(6)],
            )
            matrix = sp.csr_matrix(
                [
                    [1, 0, 2],
                    [0, 3, 0],
                    [4, 0, 1],
                    [1, 1, 0],
                    [2, 2, 2],
                    [3, 0, 1],
                ],
                dtype=np.int64,
            )
            h5ad_path = root / "input.h5ad"
            ad.AnnData(
                X=matrix,
                obs=obs,
                var=pd.DataFrame(index=["G1", "G2", "G3"]),
            ).write_h5ad(h5ad_path)

            config = {
                "biological_replicates": {
                    "canonical_col": "biological_replicate_id",
                    "strict": True,
                },
                "pseudobulk_de": {
                    "count_layer": None,
                    "replicate_col": "biological_replicate_id",
                    "celltype_col": "cell_type",
                    "age_col": "age",
                    "sex_col": "sex",
                    "batch_col": "batch",
                    "sample_col": "sample_id",
                    "min_cells_per_pseudobulk": 1,
                    "min_replicates_per_celltype": 2,
                    "min_age_span_years": 10,
                    "min_residual_df": 0,
                    "aggregation_chunk_size": 2,
                    "primary_formula": "~ sex + batch + age_decade",
                    "age_effect_scale": "log2_fold_change_per_10_years",
                    "primary_cell_types": ["A"],
                    "exploratory_cell_types": ["B"],
                },
            }
            config_path = root / "config.yml"
            config_path.write_text(
                yaml.safe_dump(config, sort_keys=False),
                encoding="utf-8",
            )
            outputs = {
                "matrix": root / "counts.mtx.gz",
                "profiles": root / "profiles.csv",
                "genes": root / "genes.csv",
                "eligibility": root / "eligibility.csv",
                "audit": root / "audit.json",
            }

            run(
                config_path=config_path,
                input_h5ad=h5ad_path,
                matrix_out=outputs["matrix"],
                profiles_out=outputs["profiles"],
                genes_out=outputs["genes"],
                eligibility_out=outputs["eligibility"],
                audit_out=outputs["audit"],
            )

            with gzip.open(outputs["matrix"], "rb") as handle:
                pseudobulk = mmread(handle).tocsr()
            profiles = pd.read_csv(outputs["profiles"])
            audit = json.loads(outputs["audit"].read_text(encoding="utf-8"))

            self.assertEqual(pseudobulk.shape, (4, 3))
            np.testing.assert_array_equal(
                np.asarray(pseudobulk.sum(axis=1)).ravel(),
                profiles["library_size"].to_numpy(),
            )
            self.assertEqual(int(pseudobulk.sum()), int(matrix.sum()))
            self.assertTrue(audit["count_conservation_passed"])
            self.assertEqual(audit["n_primary_profiles"], 2)
            self.assertEqual(audit["n_exploratory_profiles"], 2)

    def test_docker_command_contains_pinned_runtime_contract(self) -> None:
        repository = Path.cwd().resolve()
        paths = {
            "matrix": repository / ".tmp" / "counts.mtx.gz",
            "profiles": repository / ".tmp" / "profiles.csv",
            "genes": repository / ".tmp" / "genes.csv",
            "combined_out": repository / ".tmp" / "results.csv.gz",
            "celltype_dir": repository / ".tmp" / "by_cell_type",
            "manifest_out": repository / ".tmp" / "manifest.csv",
            "diagnostics_out": repository / ".tmp" / "diagnostics.csv",
            "plot_out": repository / ".tmp" / "diagnostics.pdf",
            "runtime_versions_out": repository / ".tmp" / "versions.csv",
            "session_info_out": repository / ".tmp" / "session.txt",
        }
        cfg = {
            "pseudobulk_de": {
                "primary_formula": "~ sex + batch + age_decade",
                "min_replicates_per_celltype": 12,
                "min_age_span_years": 12,
                "min_residual_df": 5,
                "runtime": {
                    "mode": "docker",
                    "docker_context": "desktop-linux",
                    "docker_image": ("immune-aging-edger:bioc-3.23-edger-4.10.1"),
                    "expected_r_version": "4.6",
                    "expected_bioconductor_version": "3.23",
                    "expected_edger_version": "4.10.1",
                },
            }
        }

        command = build_command(cfg, paths, repository)

        self.assertEqual(
            command[:5],
            ["docker", "--context", "desktop-linux", "run", "--rm"],
        )
        self.assertIn("immune-aging-edger:bioc-3.23-edger-4.10.1", command)
        self.assertEqual(
            command[command.index("--expected-edger-version") + 1],
            "4.10.1",
        )
        self.assertEqual(
            command[command.index("--primary-formula") + 1],
            "~ sex + batch + age_decade",
        )

    @unittest.skipUnless(
        os.environ.get("RUN_EDGER_DOCKER_TESTS") == "1",
        "set RUN_EDGER_DOCKER_TESTS=1 to qualify the pinned edgeR runtime",
    )
    def test_pinned_edger_runtime_on_bounded_fixture(self) -> None:
        with self._temporary_directory() as tmp:
            root = Path(tmp)
            random = np.random.default_rng(17)
            counts = random.poisson(30, size=(12, 40)).astype(np.int64)
            counts[:, 0] += np.arange(12, dtype=np.int64) * 4
            counts[:, 1] += np.arange(11, -1, -1, dtype=np.int64) * 3
            matrix_path = root / "counts.mtx.gz"
            _write_matrix_market_gzip_atomic(
                sp.csr_matrix(counts),
                matrix_path,
            )

            profiles = pd.DataFrame(
                {
                    "profile_id": [f"PB{index:05d}" for index in range(1, 13)],
                    "biological_replicate_id": [f"R{index:02d}" for index in range(1, 13)],
                    "cell_type": ["Fixture cells"] * 12,
                    "age": np.arange(20, 80, 5),
                    "age_decade": np.arange(20, 80, 5) / 10,
                    "sex": ["F", "M"] * 6,
                    "batch": ["B1", "B2", "B3"] * 4,
                    "analysis_tier": ["primary"] * 12,
                    "n_cells": [100] * 12,
                    "n_technical_libraries": [1] * 12,
                    "library_size": counts.sum(axis=1),
                }
            )
            profiles_path = root / "profiles.csv"
            profiles.to_csv(profiles_path, index=False)
            genes_path = root / "genes.csv"
            pd.DataFrame({"gene_id": [f"G{index:03d}" for index in range(1, 41)]}).to_csv(
                genes_path, index=False
            )

            config = {
                "pseudobulk_de": {
                    "primary_formula": "~ sex + batch + age_decade",
                    "min_replicates_per_celltype": 12,
                    "min_age_span_years": 12,
                    "min_residual_df": 5,
                    "runtime": {
                        "mode": "docker",
                        "docker_context": "desktop-linux",
                        "docker_image": ("immune-aging-edger:bioc-3.23-edger-4.10.1"),
                        "expected_r_version": "4.6",
                        "expected_bioconductor_version": "3.23",
                        "expected_edger_version": "4.10.1",
                        "timeout_seconds": 300,
                    },
                }
            }
            paths = {
                "matrix": matrix_path,
                "profiles": profiles_path,
                "genes": genes_path,
                "combined_out": root / "results.csv.gz",
                "celltype_dir": root / "by_cell_type",
                "manifest_out": root / "manifest.csv",
                "diagnostics_out": root / "diagnostics.csv",
                "plot_out": root / "diagnostics.pdf",
                "runtime_versions_out": root / "versions.csv",
                "session_info_out": root / "session.txt",
            }
            paths["celltype_dir"].mkdir()
            command = build_command(config, paths, Path.cwd().resolve())
            run_command(
                command,
                repository=Path.cwd().resolve(),
                log_path=root / "edgeR.log",
                timeout_seconds=300,
            )

            results = pd.read_csv(paths["combined_out"])
            diagnostics = pd.read_csv(paths["diagnostics_out"])
            versions = pd.read_csv(paths["runtime_versions_out"])
            self.assertEqual(len(results), 40)
            self.assertEqual(
                diagnostics.loc[0, "model_scope"],
                "adjusted_primary",
            )
            self.assertEqual(
                versions.loc[versions["component"].eq("edgeR"), "observed"].iloc[0],
                "4.10.1",
            )
            self.assertTrue(paths["plot_out"].stat().st_size > 0)


if __name__ == "__main__":
    unittest.main()
