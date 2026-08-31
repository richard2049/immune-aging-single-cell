from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml

from src.annotation_mapping import (
    attach_analysis_cell_types,
    file_sha256,
    load_label_mapping,
)
from src.assert_annotation_qualification import assert_approved
from src.init_annotation_mapping_review import initialize_review
from src.qualify_annotation_design import qualify


class AnnotationQualificationTests(unittest.TestCase):
    @staticmethod
    def _temporary_directory() -> tempfile.TemporaryDirectory:
        root = Path(".tmp/tests")
        root.mkdir(parents=True, exist_ok=True)
        return tempfile.TemporaryDirectory(dir=root)

    @staticmethod
    def _write_label_mapping(path: Path, status: str = "pending_review") -> None:
        pd.DataFrame(
            [
                {
                    "raw_label": "A",
                    "analysis_label": "T cells",
                    "disposition": "primary",
                    "rationale": "Fixture broad lineage",
                    "review_status": status,
                    "mapping_version": "fixture-v1",
                },
                {
                    "raw_label": "B",
                    "analysis_label": "Rare subtype",
                    "disposition": "exploratory",
                    "rationale": "Fixture sparse subtype",
                    "review_status": status,
                    "mapping_version": "fixture-v1",
                },
            ]
        ).to_csv(path, index=False)

    @staticmethod
    def _write_fixture(root: Path) -> tuple[Path, Path, Path, Path]:
        sample_units = [f"s{index}" for index in range(1, 7)]
        subject_by_sample = {
            "s1": "d1",
            "s2": "d1",
            "s3": "d2",
            "s4": "d2",
            "s5": "d3",
            "s6": "d3",
        }
        age_by_sample = dict(zip(sample_units, [30, 31, 45, 46, 60, 61]))
        sex_by_sample = dict(zip(sample_units, ["F", "F", "M", "M", "F", "F"]))
        batch_by_sample = dict(zip(sample_units, ["B1", "B2"] * 3))

        rows: list[dict] = []
        for sample in sample_units:
            for raw_label in ("A", "A", "B"):
                rows.append(
                    {
                        "cell_id": f"cell_{len(rows):02d}",
                        "subject_id": subject_by_sample[sample],
                        "sample_unit_id": sample,
                        "technical_library_id": f"lib_{sample}",
                        "age": age_by_sample[sample],
                        "sex": sex_by_sample[sample],
                        "batch": batch_by_sample[sample],
                        "cell_type": raw_label,
                    }
                )
        cells = pd.DataFrame(rows)
        obs = cells[["cell_type"]].copy()
        obs["cell_type_confidence"] = np.where(obs["cell_type"].eq("A"), 0.9, 0.4)
        obs["leiden"] = np.where(obs["cell_type"].eq("A"), "0", "1")
        obs.index = pd.Index(cells["cell_id"], name="cell_id")
        matrix = sp.csr_matrix(np.tile([2, 1, 0], (len(obs), 1)), dtype=np.int32)
        checkpoint = ad.AnnData(
            X=matrix,
            obs=obs,
            var=pd.DataFrame(index=["CD3D", "TRAC", "MS4A1"]),
        )
        checkpoint_path = root / "checkpoint.h5ad"
        checkpoint.write_h5ad(checkpoint_path)

        unit_mapping_path = root / "units.csv"
        cells[
            [
                "cell_id",
                "subject_id",
                "sample_unit_id",
                "technical_library_id",
                "age",
                "sex",
                "batch",
            ]
        ].sample(frac=1.0, random_state=7).to_csv(unit_mapping_path, index=False)

        label_mapping_path = root / "labels.csv"
        AnnotationQualificationTests._write_label_mapping(label_mapping_path)
        config = {
            "annotation_qualification": {
                "enabled": True,
                "mapping_path": str(label_mapping_path),
                "raw_label_col": "cell_type",
                "confidence_col": "cell_type_confidence",
                "cluster_col": "leiden",
                "analysis_label_col": "cell_type_analysis",
                "require_approved_mapping": True,
                "low_confidence_reference": 0.5,
                "marker_chunk_size": 5,
                "marker_panels": {"T cells": ["CD3D", "TRAC"]},
            },
            "composition_age": {"min_donors_per_celltype": 2},
            "signature_age": {
                "min_cells_per_group": 2,
                "min_donors_per_celltype": 2,
                "min_age_span": 5,
            },
            "age_prediction": {
                "min_cells_per_group": 2,
                "min_groups_for_mae_plot": 3,
                "cv_folds": 2,
            },
            "pseudobulk_de": {
                "min_cells_per_pseudobulk": 2,
                "min_replicates_per_celltype": 3,
                "min_subjects_per_celltype": 2,
                "min_age_span_years": 5,
                "min_residual_df": 1,
            },
        }
        config_path = root / "config.yaml"
        config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
        return checkpoint_path, unit_mapping_path, label_mapping_path, config_path

    def test_mapping_requires_exactly_one_row_per_observed_label(self) -> None:
        with self._temporary_directory() as tmp:
            path = Path(tmp) / "labels.csv"
            self._write_label_mapping(path)
            with self.assertRaisesRegex(ValueError, "missing=.*C"):
                load_label_mapping(path, observed_labels={"A", "B", "C"})

            duplicated = pd.read_csv(path)
            pd.concat([duplicated, duplicated.iloc[[0]]]).to_csv(path, index=False)
            with self.assertRaisesRegex(ValueError, "duplicate raw labels"):
                load_label_mapping(path)

    def test_pending_mapping_blocks_downstream_attachment(self) -> None:
        with self._temporary_directory() as tmp:
            mapping_path = Path(tmp) / "labels.csv"
            self._write_label_mapping(mapping_path)
            adata = ad.AnnData(
                X=np.ones((2, 1)),
                obs=pd.DataFrame({"cell_type": ["A", "B"]}, index=["c1", "c2"]),
            )
            cfg = {
                "annotation_qualification": {
                    "enabled": True,
                    "mapping_path": str(mapping_path),
                    "require_approved_mapping": True,
                }
            }
            with self.assertRaisesRegex(RuntimeError, "not approved"):
                attach_analysis_cell_types(adata, cfg)

    def test_approved_mapping_attaches_primary_labels_only(self) -> None:
        with self._temporary_directory() as tmp:
            mapping_path = Path(tmp) / "labels.csv"
            self._write_label_mapping(mapping_path, status="approved")
            adata = ad.AnnData(
                X=np.ones((2, 1)),
                obs=pd.DataFrame({"cell_type": ["A", "B"]}, index=["c1", "c2"]),
            )
            report = attach_analysis_cell_types(
                adata,
                {
                    "annotation_qualification": {
                        "enabled": True,
                        "mapping_path": str(mapping_path),
                        "require_approved_mapping": True,
                    }
                },
            )
            self.assertEqual(adata.obs["cell_type_analysis"].tolist()[0], "T cells")
            self.assertTrue(pd.isna(adata.obs["cell_type_analysis"].tolist()[1]))
            self.assertEqual(report["primary_cells"], 1)

    def test_qualification_is_outcome_blind_and_pending_mapping_stays_blocked(self) -> None:
        with self._temporary_directory() as tmp:
            root = Path(tmp)
            checkpoint, units, labels, config = self._write_fixture(root)
            outputs = {
                name: root / filename
                for name, filename in {
                    "evidence": "evidence.csv",
                    "design": "design.csv",
                    "cluster": "cluster.csv",
                    "marker": "marker.csv",
                    "raw_marker": "raw_marker.csv",
                    "review_table": "review_table.csv",
                    "review": "review.md",
                    "review_record": "review_record.yml",
                    "report": "report.json",
                    "approval": "approved.json",
                }.items()
            }
            report = qualify(
                config_path=config,
                checkpoint_path=checkpoint,
                unit_mapping_path=units,
                label_mapping_path=labels,
                evidence_path=outputs["evidence"],
                design_path=outputs["design"],
                cluster_path=outputs["cluster"],
                marker_path=outputs["marker"],
                raw_marker_path=outputs["raw_marker"],
                review_table_path=outputs["review_table"],
                review_path=outputs["review"],
                report_path=outputs["report"],
            )
            self.assertEqual(report["status"], "pending_review")
            self.assertFalse(report["approved_for_downstream"])
            self.assertEqual(report["summary"]["n_cells"], 18)
            self.assertNotIn("age_association", json.dumps(report))
            pending_review_text = outputs["review"].read_text().lower()
            self.assertIn("no age associations", pending_review_text)
            self.assertIn("before changing all mapping rows", pending_review_text)

            evidence = pd.read_csv(outputs["evidence"])
            self.assertEqual(set(evidence["label_level"]), {"raw", "proposed_analysis"})
            design = pd.read_csv(outputs["design"])
            sparse_raw = design.loc[design["label_level"].eq("raw") & design["label"].eq("B")]
            self.assertFalse(sparse_raw["passes_current_contract"].all())
            marker = pd.read_csv(outputs["marker"])
            self.assertEqual(
                set(marker["expression_state"]), {"raw counts; descriptive marker evidence only"}
            )
            raw_marker = pd.read_csv(outputs["raw_marker"])
            self.assertEqual(set(raw_marker["raw_label"]), {"A", "B"})
            self.assertFalse(raw_marker.duplicated(["raw_label", "gene"]).any())
            self.assertTrue(
                raw_marker.loc[
                    raw_marker["raw_label"].eq("A") & raw_marker["gene"].eq("CD3D"),
                    "expected_marker",
                ].item()
            )
            review_table = pd.read_csv(outputs["review_table"])
            self.assertEqual(review_table["raw_label"].tolist(), ["B", "A"])
            self.assertIn("pseudobulk_passes_current_contract", review_table.columns)
            self.assertNotIn("age_association", " ".join(review_table.columns))

            review_record = initialize_review(
                mapping_path=labels,
                review_table_path=outputs["review_table"],
                qualification_report_path=outputs["report"],
                output_path=outputs["review_record"],
            )
            self.assertEqual(len(review_record["decisions"]), 2)
            self.assertFalse(review_record["review_metadata"]["outcome_blind_confirmation"])
            self.assertEqual(
                review_record["composition_denominator_decision"]["status"],
                "pending_review",
            )
            with self.assertRaises(FileExistsError):
                initialize_review(
                    mapping_path=labels,
                    review_table_path=outputs["review_table"],
                    qualification_report_path=outputs["report"],
                    output_path=outputs["review_record"],
                )
            with self.assertRaisesRegex(RuntimeError, "not approved"):
                assert_approved(
                    report_path=outputs["report"],
                    mapping_path=labels,
                    output_path=outputs["approval"],
                )
            self.assertFalse(outputs["approval"].exists())

            self._write_label_mapping(labels, status="approved")
            with self.assertRaisesRegex(RuntimeError, "changed after"):
                assert_approved(
                    report_path=outputs["report"],
                    mapping_path=labels,
                    output_path=outputs["approval"],
                )

            approved_report = qualify(
                config_path=config,
                checkpoint_path=checkpoint,
                unit_mapping_path=units,
                label_mapping_path=labels,
                evidence_path=outputs["evidence"],
                design_path=outputs["design"],
                cluster_path=outputs["cluster"],
                marker_path=outputs["marker"],
                raw_marker_path=outputs["raw_marker"],
                review_table_path=outputs["review_table"],
                review_path=outputs["review"],
                report_path=outputs["report"],
            )
            self.assertTrue(approved_report["approved_for_downstream"])
            approved_review_text = outputs["review"].read_text().lower()
            self.assertIn("completed human review", approved_review_text)
            self.assertNotIn("before changing all mapping rows", approved_review_text)
            approval = assert_approved(
                report_path=outputs["report"],
                mapping_path=labels,
                output_path=outputs["approval"],
            )
            self.assertEqual(approval["status"], "approved")
            self.assertEqual(
                approval["qualification_report_sha256"],
                file_sha256(outputs["report"]),
            )
            self.assertTrue(outputs["approval"].is_file())


if __name__ == "__main__":
    unittest.main()
