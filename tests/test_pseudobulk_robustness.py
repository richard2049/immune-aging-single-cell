from __future__ import annotations

import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import scipy.sparse as sp

from src.prepare_pseudobulk_evidence import (
    _validate_sensitivity_contract,
    _write_csv,
    _write_json,
    build_candidate_evidence,
    build_review_queue,
)
from src.pseudobulk_aggregate import _write_matrix_market_gzip_atomic
from src.run_edger import run_command
from src.run_pseudobulk_robustness import build_command


class PseudobulkRobustnessTests(unittest.TestCase):
    @staticmethod
    def _temporary_directory() -> tempfile.TemporaryDirectory:
        root = Path(".tmp") / "tests"
        root.mkdir(parents=True, exist_ok=True)
        return tempfile.TemporaryDirectory(dir=root)

    @staticmethod
    def _primary_results() -> pd.DataFrame:
        common = {
            "model_scope": "adjusted_primary",
            "model_formula": "~ sex + batch + age_decade",
            "average_log_cpm": 5.0,
            "ql_f_statistic": 12.0,
            "p_value": 0.001,
            "fdr_within_celltype": 0.02,
            "n_replicates": 30,
            "age_span_years": 40.0,
        }
        return pd.DataFrame(
            [
                {
                    **common,
                    "gene_id": "GENE_UP",
                    "cell_type": "Cell A",
                    "analysis_tier": "primary",
                    "log2_fc_per_10_years": 1.0,
                    "fdr_global": 0.01,
                },
                {
                    **common,
                    "gene_id": "GENE_DOWN",
                    "cell_type": "Cell A",
                    "analysis_tier": "primary",
                    "log2_fc_per_10_years": -0.5,
                    "fdr_global": 0.04,
                },
                {
                    **common,
                    "gene_id": "NOT_SELECTED",
                    "cell_type": "Cell A",
                    "analysis_tier": "primary",
                    "log2_fc_per_10_years": 0.2,
                    "fdr_global": 0.20,
                },
                {
                    **common,
                    "gene_id": "NK_GENE",
                    "cell_type": "NK cells",
                    "analysis_tier": "exploratory",
                    "log2_fc_per_10_years": 2.0,
                    "fdr_global": 0.001,
                },
            ]
        )

    @staticmethod
    def _sensitivity_results() -> pd.DataFrame:
        effects = {
            "GENE_UP": {
                "high_cell_support": 0.8,
                "drop_sex": 1.2,
                "drop_batch": -0.2,
                "leave_one_batch_out__B1": 0.9,
            },
            "GENE_DOWN": {
                "high_cell_support": -0.4,
                "drop_sex": -0.6,
                "drop_batch": -0.3,
                "leave_one_batch_out__B1": -0.5,
            },
            "NK_GENE": {
                "high_cell_support": 1.5,
                "drop_sex": 2.2,
                "drop_batch": 1.8,
                "leave_one_batch_out__B1": 1.9,
            },
        }
        records = []
        for gene, scenarios in effects.items():
            cell_type = "NK cells" if gene == "NK_GENE" else "Cell A"
            tier = "exploratory" if gene == "NK_GENE" else "primary"
            for scenario, effect in scenarios.items():
                records.append(
                    {
                        "gene_id": gene,
                        "cell_type": cell_type,
                        "analysis_tier": tier,
                        "scenario_id": scenario,
                        "scenario_type": "leave_one_batch_out"
                        if scenario.startswith("leave_one_batch_out")
                        else "support"
                        if scenario == "high_cell_support"
                        else "covariate",
                        "excluded_batch": "B1"
                        if scenario.startswith("leave_one_batch_out")
                        else "",
                        "model_formula": "~ sex + batch + age_decade",
                        "log2_fc_per_10_years": effect,
                        "average_log_cpm": 5.0,
                        "ql_f_statistic": 8.0,
                        "p_value": 0.01,
                        "n_replicates": 25,
                        "age_span_years": 35.0,
                    }
                )
        return pd.DataFrame(records)

    @staticmethod
    def _diagnostics() -> pd.DataFrame:
        scenarios = [
            ("high_cell_support", "support"),
            ("drop_sex", "covariate"),
            ("drop_batch", "covariate"),
            ("leave_one_batch_out__B1", "leave_one_batch_out"),
        ]
        return pd.DataFrame(
            [
                {
                    "cell_type": cell_type,
                    "analysis_tier": tier,
                    "scenario_id": scenario,
                    "scenario_type": scenario_type,
                    "excluded_batch": "B1" if scenario_type == "leave_one_batch_out" else "",
                    "model_formula": "~ sex + batch + age_decade",
                    "n_replicates": 25,
                    "age_span_years": 35.0,
                    "design_columns": 5,
                    "design_rank": 5,
                    "residual_df": 20,
                    "n_genes_tested": 100,
                    "n_candidates_expected": 2 if cell_type == "Cell A" else 1,
                    "n_candidates_observed": 2 if cell_type == "Cell A" else 1,
                    "status": "completed",
                    "reason": "",
                }
                for cell_type, tier in [
                    ("Cell A", "primary"),
                    ("NK cells", "exploratory"),
                ]
                for scenario, scenario_type in scenarios
            ]
        )

    def test_candidate_evidence_preserves_global_fdr_set_and_known_effects(self) -> None:
        evidence = build_candidate_evidence(
            self._primary_results(),
            self._sensitivity_results(),
            self._diagnostics(),
            global_fdr_threshold=0.05,
        )

        self.assertEqual(set(evidence["gene_id"]), {"GENE_UP", "GENE_DOWN", "NK_GENE"})
        up = evidence.loc[evidence["gene_id"].eq("GENE_UP")].iloc[0]
        self.assertAlmostEqual(up["direction_concordance_fraction"], 0.75)
        self.assertAlmostEqual(up["leave_one_batch_out_direction_concordance_fraction"], 1.0)
        self.assertAlmostEqual(up["max_absolute_effect_change"], 1.2)
        self.assertTrue(evidence["requires_human_review"].all())
        self.assertEqual(set(evidence["automated_disposition"]), {"not_assigned"})
        self.assertEqual(
            evidence.loc[evidence["gene_id"].eq("NK_GENE"), "analysis_tier"].iloc[0],
            "exploratory",
        )

    def test_review_queue_is_bounded_by_celltype_and_direction(self) -> None:
        evidence = build_candidate_evidence(
            self._primary_results(),
            self._sensitivity_results(),
            self._diagnostics(),
            global_fdr_threshold=0.05,
        )
        queue = build_review_queue(evidence, limit=1)

        counts = queue.groupby(["cell_type", "effect_direction"]).size()
        self.assertTrue(counts.le(1).all())
        self.assertIn("NK_GENE", set(queue["gene_id"]))

    def test_non_candidate_sensitivity_row_is_rejected(self) -> None:
        sensitivity = self._sensitivity_results()
        extra = sensitivity.iloc[[0]].copy()
        extra["gene_id"] = "NOT_SELECTED"

        with self.assertRaisesRegex(ValueError, "non-candidate"):
            build_candidate_evidence(
                self._primary_results(),
                pd.concat([sensitivity, extra], ignore_index=True),
                self._diagnostics(),
                global_fdr_threshold=0.05,
            )

    def test_compressed_atomic_csv_is_readable(self) -> None:
        with self._temporary_directory() as tmp:
            path = Path(tmp) / "evidence.csv.gz"
            expected = pd.DataFrame({"gene_id": ["G1"], "effect": [0.25]})

            _write_csv(expected, path)

            self.assertEqual(path.read_bytes()[:2], b"\x1f\x8b")
            pd.testing.assert_frame_equal(pd.read_csv(path), expected)

    def test_json_writer_converts_numpy_scalars(self) -> None:
        with self._temporary_directory() as tmp:
            path = Path(tmp) / "summary.json"

            _write_json({"passed": np.bool_(True), "models": np.int64(3)}, path)

            self.assertEqual(
                path.read_text(encoding="utf-8"), '{\n  "models": 3,\n  "passed": true\n}\n'
            )

    def test_sensitivity_contract_rejects_diagnostic_count_mismatch(self) -> None:
        diagnostics = self._diagnostics()
        diagnostics.loc[0, "n_candidates_observed"] += 1
        manifest = diagnostics[
            [
                "scenario_id",
                "scenario_type",
                "model_formula",
                "excluded_batch",
            ]
        ].drop_duplicates()
        manifest["min_cells"] = 0

        with self.assertRaisesRegex(ValueError, "disagree"):
            _validate_sensitivity_contract(
                self._primary_results(),
                self._sensitivity_results(),
                diagnostics,
                manifest,
                global_fdr_threshold=0.05,
            )

    @patch("src.run_pseudobulk_robustness.shutil.which", return_value="docker")
    def test_docker_command_preserves_approved_contract(self, _: object) -> None:
        repository = Path.cwd().resolve()
        cfg = {
            "pseudobulk_de": {
                "primary_formula": "~ sex + batch + age_decade",
                "min_replicates_per_celltype": 12,
                "min_age_span_years": 12,
                "min_residual_df": 5,
                "runtime": {
                    "mode": "docker",
                    "docker_context": "desktop-linux",
                    "docker_image": "test-edger",
                    "expected_r_version": "4.6",
                    "expected_bioconductor_version": "3.23",
                    "expected_edger_version": "4.10.1",
                },
            },
            "pseudobulk_robustness": {
                "candidate_global_fdr_threshold": 0.05,
                "high_cell_min_cells": 100,
                "run_leave_one_batch_out": True,
            },
        }
        paths = {
            key: repository / "results" / value
            for key, value in {
                "matrix": "counts.mtx.gz",
                "profiles": "profiles.csv",
                "genes": "genes.csv",
                "primary_results": "primary.csv.gz",
                "sensitivity_out": "sensitivity.csv.gz",
                "diagnostics_out": "diagnostics.csv",
                "manifest_out": "manifest.csv",
            }.items()
        }

        command = build_command(cfg, paths, repository)

        self.assertIn("/work/src/pseudobulk_robustness.R", command)
        self.assertIn("~ sex + batch + age_decade", command)
        self.assertEqual(command[command.index("--candidate-fdr") + 1], "0.05")
        self.assertEqual(command[command.index("--high-cell-min") + 1], "100")
        self.assertEqual(command[command.index("--run-leave-one-batch-out") + 1], "true")

    @unittest.skipUnless(
        os.environ.get("RUN_EDGER_DOCKER_TESTS") == "1",
        "set RUN_EDGER_DOCKER_TESTS=1 to test robustness models in Docker",
    )
    def test_bounded_robustness_models_preserve_candidate_scope(self) -> None:
        with self._temporary_directory() as tmp:
            root = Path(tmp)
            random = np.random.default_rng(31)
            counts = random.poisson(30, size=(36, 50)).astype(np.int64)
            counts[:, 0] += np.arange(36, dtype=np.int64) * 3
            matrix_path = root / "counts.mtx.gz"
            _write_matrix_market_gzip_atomic(sp.csr_matrix(counts), matrix_path)

            profiles = pd.DataFrame(
                {
                    "profile_id": [f"PB{i:05d}" for i in range(1, 37)],
                    "biological_replicate_id": [f"R{i:03d}" for i in range(1, 37)],
                    "cell_type": ["Fixture cells"] * 36,
                    "age": np.linspace(20, 75, 36),
                    "age_decade": np.linspace(20, 75, 36) / 10,
                    "sex": ["F", "M"] * 18,
                    "batch": ["B1", "B2", "B3"] * 12,
                    "analysis_tier": ["primary"] * 36,
                    "n_cells": [120] * 36,
                    "n_technical_libraries": [1] * 36,
                    "library_size": counts.sum(axis=1),
                }
            )
            profiles_path = root / "profiles.csv"
            profiles.to_csv(profiles_path, index=False)
            genes_path = root / "genes.csv"
            pd.DataFrame({"gene_id": [f"G{i:03d}" for i in range(1, 51)]}).to_csv(
                genes_path, index=False
            )
            primary_path = root / "primary.csv.gz"
            pd.DataFrame(
                {
                    "gene_id": ["G001", "G002"],
                    "cell_type": ["Fixture cells", "Fixture cells"],
                    "analysis_tier": ["primary", "primary"],
                    "log2_fc_per_10_years": [0.5, 0.0],
                    "fdr_global": [0.01, 0.50],
                }
            ).to_csv(primary_path, index=False)

            config = {
                "pseudobulk_de": {
                    "primary_formula": "~ sex + batch + age_decade",
                    "min_replicates_per_celltype": 12,
                    "min_age_span_years": 12,
                    "min_residual_df": 5,
                    "runtime": {
                        "mode": "docker",
                        "docker_context": "desktop-linux",
                        "docker_image": "immune-aging-edger:bioc-3.23-edger-4.10.1",
                        "expected_r_version": "4.6",
                        "expected_bioconductor_version": "3.23",
                        "expected_edger_version": "4.10.1",
                    },
                },
                "pseudobulk_robustness": {
                    "candidate_global_fdr_threshold": 0.05,
                    "high_cell_min_cells": 100,
                    "run_leave_one_batch_out": True,
                },
            }
            paths = {
                "matrix": matrix_path,
                "profiles": profiles_path,
                "genes": genes_path,
                "primary_results": primary_path,
                "sensitivity_out": root / "sensitivity.csv.gz",
                "diagnostics_out": root / "diagnostics.csv",
                "manifest_out": root / "manifest.csv",
            }
            run_command(
                build_command(config, paths, Path.cwd().resolve()),
                repository=Path.cwd().resolve(),
                log_path=root / "robustness.log",
                timeout_seconds=300,
            )

            sensitivity = pd.read_csv(paths["sensitivity_out"])
            diagnostics = pd.read_csv(paths["diagnostics_out"])
            self.assertEqual(set(sensitivity["gene_id"]), {"G001"})
            self.assertEqual(len(diagnostics), 6)
            self.assertTrue(diagnostics["status"].eq("completed").all())
            self.assertEqual(sensitivity["scenario_id"].nunique(), 6)


if __name__ == "__main__":
    unittest.main()
