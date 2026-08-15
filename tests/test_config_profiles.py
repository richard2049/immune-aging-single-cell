from __future__ import annotations

import unittest
from pathlib import Path

import yaml

from src.utils import load_config

CONFIG_DIR = Path("config")
SUPPORTED_PROFILES = (
    Path("config/custom.example.yaml"),
    Path("config/demo.yaml"),
    Path("config/gse164378_1m.yaml"),
    Path("config/gse164378_corrected.yaml"),
    Path("config/gse164378_pilot.yaml"),
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

    def test_gse164378_pilot_is_bounded_and_separate(self) -> None:
        config = load_config("config/gse164378_pilot.yaml")
        self.assertEqual(config["run"]["max_cells"], 50_000)
        self.assertEqual(config["project"]["out_dir"], "results/gse164378_pilot")

    def test_gse164378_1m_preserves_accepted_cell_cap(self) -> None:
        config = load_config("config/gse164378_1m.yaml")
        self.assertEqual(config["run"]["max_cells"], 1_000_000)
        self.assertEqual(config["project"]["out_dir"], "results/gse164378_full")

    def test_corrected_profile_preserves_replicate_contract(self) -> None:
        raw = _read_raw(Path("config/gse164378_corrected.yaml"))
        self.assertEqual(raw["extends"], "gse164378_1m.yaml")

        config = load_config("config/gse164378_corrected.yaml")
        self.assertEqual(
            config["analysis_checkpoint"]["annotated_h5ad"],
            "results/gse164378_full/06_annotated.h5ad",
        )
        self.assertEqual(
            config["project"]["out_dir"],
            "results/gse164378_full_replicate_corrected",
        )
        for section in ("composition_age", "signature_age", "age_prediction"):
            self.assertEqual(config[section]["donor_col"], "biological_replicate_id")
        self.assertEqual(
            config["pseudobulk_de"]["replicate_col"],
            "biological_replicate_id",
        )

    def test_custom_example_uses_safe_generic_identity_and_paths(self) -> None:
        config = load_config("config/custom.example.yaml")
        self.assertEqual(config["run"]["dataset"], "custom_h5ad")
        self.assertEqual(config["paths"]["input_h5ad"], "data/raw/input.h5ad")
        self.assertEqual(config["project"]["out_dir"], "results/custom")
        for section in ("composition_age", "signature_age", "age_prediction"):
            self.assertEqual(config[section]["donor_col"], "biological_replicate_id")


if __name__ == "__main__":
    unittest.main()
