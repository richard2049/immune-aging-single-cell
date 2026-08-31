from __future__ import annotations

import unittest
from pathlib import Path

import yaml

from src.cluster import _required_clustering_arguments
from src.scvi_train import _required_arguments
from src.utils import load_config

CONFIG_DIR = Path("config")
SUPPORTED_PROFILES = (
    Path("config/blood_age_atlas_1m.yaml"),
    Path("config/blood_age_atlas_longitudinal.yaml"),
    Path("config/blood_age_atlas_pilot.yaml"),
    Path("config/custom.example.yaml"),
    Path("config/demo.yaml"),
)


def _read_raw(path: Path) -> dict:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict):
        raise TypeError(f"Configuration must be a mapping: {path}")
    return config


class ConfigProfileContractTests(unittest.TestCase):
    def test_supported_profile_inventory_is_explicit(self) -> None:
        observed = tuple(
            sorted(path for path in CONFIG_DIR.glob("*.y*ml") if ".local." not in path.name)
        )
        self.assertEqual(observed, SUPPORTED_PROFILES)
        self.assertFalse(Path("config/config.yaml").exists())
        self.assertFalse(Path("config/config.real.yml").exists())
        self.assertFalse(Path("config/config.real.full.yml").exists())

    def test_profiles_point_to_their_tracked_paths(self) -> None:
        for path in SUPPORTED_PROFILES:
            with self.subTest(path=path):
                self.assertEqual(_read_raw(path)["cfg_path"], path.as_posix())

    def test_resolved_profiles_define_seed_and_unique_output(self) -> None:
        outputs: dict[str, Path] = {}
        for path in SUPPORTED_PROFILES:
            with self.subTest(path=path):
                config = load_config(path)
                self.assertIsInstance(config["run"]["seed"], int)
                out_dir = config["project"]["out_dir"]
                self.assertNotIn(out_dir, outputs, msg=f"also used by {outputs.get(out_dir)}")
                outputs[out_dir] = path

    def test_only_demo_profile_allows_placeholder_outputs(self) -> None:
        for path in SUPPORTED_PROFILES:
            with self.subTest(path=path):
                config = load_config(path)
                expected = path == Path("config/demo.yaml")
                self.assertIs(config["run"]["allow_placeholder_outputs"], expected)

    def test_scvi_execution_parameters_are_explicit(self) -> None:
        for path in SUPPORTED_PROFILES:
            with self.subTest(path=path):
                config = load_config(path)
                model_args, training_args = _required_arguments(config["scvi"])
                self.assertGreater(int(model_args["n_latent"]), 0)
                self.assertGreater(int(training_args["max_epochs"]), 0)

    def test_clustering_execution_parameters_are_explicit(self) -> None:
        for path in SUPPORTED_PROFILES:
            with self.subTest(path=path):
                config = load_config(path)
                neighbors, umap, leiden = _required_clustering_arguments(config["clustering"])
                self.assertEqual(int(config["clustering"]["n_jobs"]), 1)
                self.assertGreater(int(neighbors["n_neighbors"]), 0)
                self.assertEqual(int(umap["n_components"]), 2)
                self.assertEqual(str(leiden["flavor"]), "leidenalg")
                self.assertEqual(int(leiden["n_iterations"]), 2)

    def test_incomplete_clustering_contract_fails(self) -> None:
        with self.assertRaisesRegex(KeyError, "clustering execution contract"):
            _required_clustering_arguments({"neighbors_args": {"n_neighbors": 15}})

    def test_blood_age_atlas_pilot_is_bounded_and_separate(self) -> None:
        config = load_config("config/blood_age_atlas_pilot.yaml")
        self.assertEqual(config["run"]["max_cells"], 50_000)
        self.assertEqual(config["project"]["out_dir"], "results/blood_age_atlas_pilot")
        self.assertEqual(config["study"]["synapse_accession"], "syn49637038")

    def test_blood_age_atlas_1m_preserves_accepted_cell_cap(self) -> None:
        config = load_config("config/blood_age_atlas_1m.yaml")
        self.assertEqual(config["run"]["max_cells"], 1_000_000)
        self.assertEqual(config["project"]["out_dir"], "results/blood_age_atlas_1m")
        self.assertEqual(config["biological_units"]["subject_col"], "donor_id")
        self.assertEqual(config["biological_units"]["sample_unit_col"], "tube_id")
        self.assertEqual(config["age_prediction"]["group_col"], "donor_id")

    def test_longitudinal_profile_preserves_unit_contract(self) -> None:
        raw = _read_raw(Path("config/blood_age_atlas_longitudinal.yaml"))
        self.assertEqual(raw["extends"], "blood_age_atlas_1m.yaml")

        config = load_config("config/blood_age_atlas_longitudinal.yaml")
        self.assertEqual(
            config["analysis_checkpoint"]["annotated_h5ad"],
            "results/blood_age_atlas_1m/06_annotated.h5ad",
        )
        self.assertEqual(
            config["project"]["out_dir"],
            "results/blood_age_atlas_longitudinal",
        )
        for section in ("composition_age", "signature_age"):
            self.assertEqual(config[section]["sample_unit_col"], "sample_unit_id")
            self.assertEqual(config[section]["subject_col"], "subject_id")
            self.assertEqual(config[section]["celltype_col"], "cell_type_analysis")
        self.assertEqual(
            config["composition_age"]["denominator"],
            "all_qc_passed_cells_with_other_unresolved",
        )
        self.assertEqual(
            config["composition_age"]["other_unresolved_label"],
            "Other/unresolved",
        )
        self.assertEqual(
            config["signature_age"]["gee_working_correlation"],
            "exchangeable",
        )
        self.assertEqual(
            config["signature_age"]["gee_nonconvergence_fallback"],
            "independence",
        )
        self.assertEqual(config["age_prediction"]["group_col"], "subject_id")
        self.assertEqual(
            config["age_prediction"]["celltype_col"],
            "cell_type_analysis",
        )
        self.assertEqual(
            config["pseudobulk_de"]["replicate_col"],
            "sample_unit_id",
        )
        self.assertEqual(config["pseudobulk_de"]["subject_col"], "subject_id")
        self.assertEqual(
            config["pseudobulk_de"]["celltype_col"],
            "cell_type_analysis",
        )
        self.assertEqual(config["pseudobulk_de"]["min_subjects_per_celltype"], 12)
        self.assertTrue(config["annotation_qualification"]["enabled"])
        self.assertTrue(config["annotation_qualification"]["require_approved_mapping"])

    def test_custom_example_uses_safe_generic_identity_and_paths(self) -> None:
        config = load_config("config/custom.example.yaml")
        self.assertEqual(config["run"]["dataset"], "custom_h5ad")
        self.assertEqual(config["paths"]["input_h5ad"], "data/raw/input.h5ad")
        self.assertEqual(config["project"]["out_dir"], "results/custom")
        for section in ("composition_age", "signature_age", "age_prediction"):
            self.assertEqual(config[section]["donor_col"], "biological_replicate_id")


if __name__ == "__main__":
    unittest.main()
