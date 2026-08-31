from __future__ import annotations

import gzip
import json
import os
import sys
import tempfile
import unittest
import zlib
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml
from scipy.io import mmread

from src.pseudobulk_aggregate import (
    _build_profiles_and_eligibility,
    _integer_sparse_chunk,
    _replicate_metadata,
    _validate_pseudobulk_metadata,
    _write_matrix_market_gzip_atomic,
    aggregate_sparse_counts,
    run,
)
from src.run_pseudobulk_dream import build_command, run_command
from src.validate_pseudobulk_de import (
    _bh_adjust,
    _eligible_analysis_sets,
    _inspect_r_graphics_pdf,
    _resolve_manifest_path,
)


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

    def test_r_graphics_pdf_integrity_rejects_transcoded_stream(self) -> None:
        drawing = zlib.compress(b"0 0 m 10 10 l S\n")
        valid = (
            b"%PDF-1.4\n"
            b"1 0 obj\n<< /Type /Page /Contents 2 0 R >>\nendobj\n"
            + f"2 0 obj\n<< /Length {len(drawing)} /Filter /FlateDecode >>\nstream\n".encode()
            + drawing
            + b"\nendstream\nendobj\n%%EOF\n"
        )
        with self._temporary_directory() as tmp:
            valid_path = Path(tmp) / "valid.pdf"
            valid_path.write_bytes(valid)
            passed, detail = _inspect_r_graphics_pdf(valid_path, expected_pages=1)
            self.assertTrue(passed, detail)

            corrupt_path = Path(tmp) / "corrupt.pdf"
            corrupt_path.write_bytes(valid.replace(drawing[1:2], b"\xef\xbf\xbd", 1))
            passed, detail = _inspect_r_graphics_pdf(corrupt_path, expected_pages=1)
            self.assertFalse(passed)
            self.assertIn("invalid compressed stream", detail)

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

    def test_nonprimary_cells_are_excluded_without_weakening_metadata_checks(self) -> None:
        obs = pd.DataFrame(
            {
                "sample_unit_id": ["S1", "S1", "S2"],
                "cell_type_analysis": ["B cells", pd.NA, "CD4 T cells"],
                "age": [40, 40, 55],
                "sex": ["F", "F", "M"],
                "batch": ["B1", "B1", "B2"],
            }
        )
        required = ["sample_unit_id", "cell_type_analysis", "age", "sex", "batch"]

        report = _validate_pseudobulk_metadata(
            obs,
            required_columns=required,
            celltype_col="cell_type_analysis",
        )

        self.assertEqual(report["analysis_eligible_cells"], 2)
        self.assertEqual(report["excluded_nonprimary_cells"], 1)

        obs.loc[1, "batch"] = pd.NA
        with self.assertRaisesRegex(ValueError, "batch=1"):
            _validate_pseudobulk_metadata(
                obs,
                required_columns=required,
                celltype_col="cell_type_analysis",
            )

    def test_bh_adjustment_matches_known_values(self) -> None:
        observed = _bh_adjust(np.array([0.01, 0.04, 0.03, 0.002]))
        np.testing.assert_allclose(observed, [0.02, 0.04, 0.04, 0.008])

    def test_method_eligibility_preserves_configured_exclusions(self) -> None:
        contract = {
            "primary_cell_types": ["A", "B"],
            "exploratory_cell_types": ["C"],
        }
        eligibility = pd.DataFrame(
            {
                "cell_type": ["A", "B", "C"],
                "analysis_tier": ["primary", "primary", "exploratory"],
                "aggregation_status": ["included", "excluded", "included"],
                "exclusion_reasons": ["", "insufficient_subjects", ""],
            }
        )

        primary, exploratory, valid = _eligible_analysis_sets(contract, eligibility)

        self.assertTrue(valid)
        self.assertEqual(primary, {"A"})
        self.assertEqual(exploratory, {"C"})

        eligibility.loc[1, "exclusion_reasons"] = ""
        _, _, valid = _eligible_analysis_sets(contract, eligibility)
        self.assertFalse(valid)

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

    def test_longitudinal_eligibility_counts_independent_subjects(self) -> None:
        obs = pd.DataFrame(
            {
                "sample_unit_id": [f"S{index}" for index in range(1, 7)],
                "subject_id": ["D1"] * 6,
                "cell_type": ["A"] * 6,
                "age": np.arange(30, 36),
                "sex": ["F"] * 6,
                "batch": ["X", "Y"] * 3,
            }
        )
        metadata = _replicate_metadata(
            obs,
            replicate_col="sample_unit_id",
            age_col="age",
            sex_col="sex",
            batch_col="batch",
            subject_col="subject_id",
        )
        _, eligibility = _build_profiles_and_eligibility(
            obs,
            metadata,
            {
                "replicate_col": "sample_unit_id",
                "subject_col": "subject_id",
                "celltype_col": "cell_type",
                "min_cells_per_pseudobulk": 1,
                "min_replicates_per_celltype": 2,
                "min_subjects_per_celltype": 2,
                "min_age_span_years": 1,
                "min_residual_df": 0,
                "primary_cell_types": ["A"],
                "exploratory_cell_types": [],
            },
        )

        row = eligibility.iloc[0]
        self.assertEqual(row["n_sample_units_qualifying"], 6)
        self.assertEqual(row["n_subjects_qualifying"], 1)
        self.assertEqual(row["aggregation_status"], "excluded")
        self.assertIn("insufficient_qualifying_subjects", row["exclusion_reasons"])

    def test_subject_sex_conflict_is_rejected(self) -> None:
        obs = pd.DataFrame(
            {
                "sample_unit_id": ["S1", "S2"],
                "subject_id": ["D1", "D1"],
                "age": [30, 31],
                "sex": ["F", "M"],
                "batch": ["X", "Y"],
            }
        )

        with self.assertRaisesRegex(ValueError, "inconsistent across sample units"):
            _replicate_metadata(
                obs,
                replicate_col="sample_unit_id",
                age_col="age",
                sex_col="sex",
                batch_col="batch",
                subject_col="subject_id",
            )

    def test_bounded_aggregation_writes_auditable_matrix_market(self) -> None:
        with self._temporary_directory() as tmp:
            root = Path(tmp)
            obs = pd.DataFrame(
                {
                    "sample_unit_id": [f"S{i}" for i in range(1, 7)] * 2,
                    "subject_id": ["D1", "D1", "D2", "D3", "D4", "D5"] * 2,
                    "technical_library_id": [f"L{i}" for i in range(1, 7)] * 2,
                    "cell_type": ["A"] * 6 + ["B"] * 6,
                    "age": [30, 31, 40, 50, 60, 70] * 2,
                    "sex": ["F", "F", "M", "F", "M", "M"] * 2,
                    "batch": ["X", "Y", "X", "Y", "Y", "X"] * 2,
                },
                index=[f"cell_{index}" for index in range(12)],
            )
            matrix = sp.csr_matrix(
                [
                    [1, 0, 2],
                    [0, 3, 0],
                    [4, 0, 1],
                    [1, 1, 0],
                    [2, 2, 2],
                    [3, 0, 1],
                    [2, 0, 1],
                    [0, 2, 1],
                    [3, 1, 0],
                    [1, 2, 0],
                    [2, 1, 2],
                    [4, 0, 1],
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
                "biological_units": {
                    "subject_col": "subject_id",
                    "subject_source_col": "subject_id",
                    "sample_unit_col": "sample_unit_id",
                    "sample_unit_source_col": "sample_unit_id",
                    "technical_library_col": "technical_library_id",
                    "technical_library_source_col": "technical_library_id",
                    "strict": True,
                },
                "pseudobulk_de": {
                    "count_layer": None,
                    "replicate_col": "sample_unit_id",
                    "subject_col": "subject_id",
                    "celltype_col": "cell_type",
                    "age_col": "age",
                    "sex_col": "sex",
                    "batch_col": "batch",
                    "sample_col": "technical_library_id",
                    "min_cells_per_pseudobulk": 1,
                    "min_replicates_per_celltype": 2,
                    "min_age_span_years": 10,
                    "min_residual_df": 0,
                    "aggregation_chunk_size": 2,
                    "primary_formula": "~ sex + batch + age_decade + (1 | subject_id)",
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

            self.assertEqual(pseudobulk.shape, (12, 3))
            np.testing.assert_array_equal(
                np.asarray(pseudobulk.sum(axis=1)).ravel(),
                profiles["library_size"].to_numpy(),
            )
            self.assertEqual(int(pseudobulk.sum()), int(matrix.sum()))
            self.assertTrue(audit["count_conservation_passed"])
            self.assertEqual(audit["n_primary_profiles"], 6)
            self.assertEqual(audit["n_exploratory_profiles"], 6)
            self.assertEqual(profiles["subject_id"].nunique(), 5)

    def test_docker_command_contains_pinned_runtime_contract(self) -> None:
        repository = Path.cwd().resolve()
        paths = {
            "matrix": repository / ".tmp" / "counts.mtx.gz",
            "profiles": repository / ".tmp" / "profiles.csv",
            "genes": repository / ".tmp" / "genes.csv",
            "combined_out": repository / ".tmp" / "results.csv.gz",
            "sensitivity_out": repository / ".tmp" / "sensitivity.csv.gz",
            "celltype_dir": repository / ".tmp" / "by_cell_type",
            "manifest_out": repository / ".tmp" / "manifest.csv",
            "diagnostics_out": repository / ".tmp" / "diagnostics.csv",
            "plot_out": repository / ".tmp" / "diagnostics.pdf",
            "runtime_versions_out": repository / ".tmp" / "versions.csv",
            "session_info_out": repository / ".tmp" / "session.txt",
        }
        cfg = {
            "pseudobulk_de": {
                "primary_formula": "~ sex + batch + age_decade + (1 | subject_id)",
                "min_replicates_per_celltype": 12,
                "min_subjects_per_celltype": 12,
                "min_age_span_years": 12,
                "min_residual_df": 5,
                "runtime": {
                    "mode": "docker",
                    "docker_context": "desktop-linux",
                    "docker_image": ("immune-aging-dream:bioc-3.23-dream-1.42.0"),
                    "expected_r_version": "4.6",
                    "expected_bioconductor_version": "3.23",
                    "expected_edger_version": "4.10.1",
                    "expected_variance_partition_version": "1.42.0",
                },
            }
        }

        command = build_command(cfg, paths, repository)

        self.assertEqual(
            command[:5],
            ["docker", "--context", "desktop-linux", "run", "--rm"],
        )
        self.assertIn("immune-aging-dream:bioc-3.23-dream-1.42.0", command)
        self.assertEqual(
            command[command.index("--expected-edger-version") + 1],
            "4.10.1",
        )
        self.assertEqual(
            command[command.index("--expected-variance-partition-version") + 1],
            "1.42.0",
        )
        self.assertEqual(
            command[command.index("--primary-formula") + 1],
            "~ sex + batch + age_decade + (1 | subject_id)",
        )
        self.assertEqual(command[command.index("--min-subjects") + 1], "12")
        self.assertEqual(
            command[command.index("--sensitivity-sample-rule") + 1],
            "earliest_age_then_sample_id",
        )

    @unittest.skipUnless(
        os.environ.get("RUN_DREAM_DOCKER_TESTS") == "1",
        "set RUN_DREAM_DOCKER_TESTS=1 to qualify the pinned dream runtime",
    )
    def test_pinned_dream_runtime_on_bounded_fixture(self) -> None:
        with self._temporary_directory() as tmp:
            root = Path(tmp)
            random = np.random.default_rng(17)
            counts = random.poisson(30, size=(24, 40)).astype(np.int64)
            counts[:, 0] += np.arange(24, dtype=np.int64) * 4
            counts[:, 1] += np.arange(23, -1, -1, dtype=np.int64) * 3
            matrix_path = root / "counts.mtx.gz"
            _write_matrix_market_gzip_atomic(
                sp.csr_matrix(counts),
                matrix_path,
            )

            profiles = pd.DataFrame(
                {
                    "profile_id": [f"PB{index:05d}" for index in range(1, 25)],
                    "sample_unit_id": [f"S{index:02d}" for index in range(1, 25)],
                    "subject_id": [f"D{((index - 1) // 2) + 1:02d}" for index in range(1, 25)],
                    "cell_type": ["Fixture cells"] * 24,
                    "age": np.repeat(np.arange(25, 85, 5), 2) + np.tile([0, 1], 12),
                    "age_decade": (np.repeat(np.arange(25, 85, 5), 2) + np.tile([0, 1], 12)) / 10,
                    "sex": np.repeat(["F"] * 6 + ["M"] * 6, 2),
                    "batch": ["B1", "B2", "B2", "B3"] * 6,
                    "analysis_tier": ["primary"] * 24,
                    "n_cells": [100] * 24,
                    "n_technical_libraries": [1] * 24,
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
                    "primary_formula": "~ sex + batch + age_decade + (1 | subject_id)",
                    "min_replicates_per_celltype": 12,
                    "min_subjects_per_celltype": 12,
                    "min_age_span_years": 12,
                    "min_residual_df": 5,
                    "runtime": {
                        "mode": "docker",
                        "docker_context": "desktop-linux",
                        "docker_image": ("immune-aging-dream:bioc-3.23-dream-1.42.0"),
                        "expected_r_version": "4.6",
                        "expected_bioconductor_version": "3.23",
                        "expected_edger_version": "4.10.1",
                        "expected_variance_partition_version": "1.42.0",
                        "timeout_seconds": 300,
                    },
                }
            }
            paths = {
                "matrix": matrix_path,
                "profiles": profiles_path,
                "genes": genes_path,
                "combined_out": root / "results.csv.gz",
                "sensitivity_out": root / "sensitivity.csv.gz",
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
                log_path=root / "dream.log",
                timeout_seconds=300,
            )

            results = pd.read_csv(paths["combined_out"])
            sensitivity = pd.read_csv(paths["sensitivity_out"])
            diagnostics = pd.read_csv(paths["diagnostics_out"])
            versions = pd.read_csv(paths["runtime_versions_out"])
            self.assertEqual(len(results), 40)
            self.assertEqual(set(sensitivity["model_scope"]), {"one_sample_per_subject"})
            for column in (
                "log2_fc_per_10_years",
                "moderated_t_statistic",
                "z_standardized",
                "p_value",
                "fdr_within_celltype",
                "fdr_global",
            ):
                self.assertTrue(np.isfinite(results[column]).all(), column)
            for column in (
                "log2_fc_per_10_years",
                "moderated_t_statistic",
                "p_value",
                "fdr_within_celltype",
                "fdr_global",
            ):
                self.assertTrue(np.isfinite(sensitivity[column]).all(), column)
            self.assertEqual(diagnostics.loc[0, "model_scope"], "adjusted_repeated_measures")
            self.assertEqual(
                versions.loc[versions["component"].eq("edgeR"), "observed"].iloc[0],
                "4.10.1",
            )
            self.assertEqual(
                versions.loc[versions["component"].eq("variancePartition"), "observed"].iloc[0],
                "1.42.0",
            )
            pdf_ok, pdf_detail = _inspect_r_graphics_pdf(paths["plot_out"], expected_pages=1)
            self.assertTrue(pdf_ok, pdf_detail)


if __name__ == "__main__":
    unittest.main()
