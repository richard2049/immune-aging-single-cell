from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest

from src.historical_output_paths import validate_historical_outputs
from src.render_reviewed_figures import (
    ASSETS,
    cd8_figure,
    cd8_inputs,
    composition_figure,
    composition_inputs,
    prediction_figure,
    prediction_inputs,
)


def bars(axis):
    return [tuple(item.lines[2][0].get_segments()[0][:, 0]) for item in axis.containers]


def data(name, delimiter=","):
    with (ASSETS / name).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter=delimiter))


def test_composition_uses_both_adjusted_models_for_all_populations():
    name = "blood_age_atlas_composition_summary.csv"
    fig = composition_figure(*composition_inputs(ASSETS / name))
    try:
        ax = fig.axes[0]
        source = data(name)
        assert [label.get_text() for label in ax.get_yticklabels()] == [
            r["cell_type"] for r in source
        ]
        series = [(r, p) for p in ("association", "sensitivity") for r in source]
        np.testing.assert_allclose(
            [float(container.lines[0].get_xdata()[0]) for container in ax.containers],
            [float(r[f"{p}_effect_per_10y"]) for r, p in series],
        )
        np.testing.assert_allclose(
            bars(ax),
            [[float(r[f"{p}_effect_ci_low"]), float(r[f"{p}_effect_ci_high"])] for r, p in series],
        )
    finally:
        plt.close(fig)


def test_prediction_plots_visit_mae_with_point_only_training_baseline():
    fig = prediction_figure(
        *prediction_inputs(ASSETS / "blood_age_atlas_age_prediction_summary.csv")
    )
    try:
        ax = fig.axes[0]
        assert "166 people, 316 visits, 5 person-grouped outer folds" in ax.get_title()
        assert len(ax.get_yticklabels()) == 3
        np.testing.assert_allclose(
            bars(ax),
            [(11.638430088973383, 13.626712926426078), (11.658658909393075, 13.876669307415586)],
        )
        np.testing.assert_allclose(
            [float(c.lines[0].get_xdata()[0]) for c in ax.containers],
            [12.640348584551813, 12.764120767699435],
        )
        assert ax.lines[-1].get_xdata()[0] == pytest.approx(15.160436692998319)
    finally:
        plt.close(fig)


def test_cd8_uses_adjusted_coefficients_and_model_specific_counts():
    fig = cd8_figure(cd8_inputs(ASSETS / "blood_age_atlas_cd8_sasp_summary.tsv"))
    try:
        ax = fig.axes[0]
        labels = [label.get_text() for label in ax.get_yticklabels()]
        assert labels[0].endswith("162 people, 309 visits")
        assert labels[1].endswith("161 people, 161 visits")
        np.testing.assert_allclose(
            [float(c.lines[0].get_xdata()[0]) for c in ax.containers],
            [0.0005871327345026147, 0.001077211107201668],
        )
        np.testing.assert_allclose(
            bars(ax),
            [
                (0.0002726933953898079, 0.0009015720736154216),
                (0.0006400393283442943, 0.0015143828860590416),
            ],
        )
    finally:
        plt.close(fig)


def test_cd8_rejects_invalid_interval(tmp_path: Path):
    path = tmp_path / "summary.tsv"
    original = (ASSETS / "blood_age_atlas_cd8_sasp_summary.tsv").read_text(encoding="utf-8")
    path.write_text(original.replace("0.0009015720736154216", "0.0001"), encoding="utf-8")
    with pytest.raises(ValueError, match="Invalid estimate or interval"):
        cd8_inputs(path)


@pytest.mark.parametrize(
    "path",
    [
        "docs/assets/old.png",
        "results/blood_age_atlas_1m/tables/old.csv",
        "results/blood_age_atlas_requalification_longitudinal/old.csv",
        "results/gse164378_1m/old.json",
        "scratch/../docs/assets/old.png",
    ],
)
def test_historical_outputs_reject_protected_paths(path):
    with pytest.raises(ValueError, match="separate scratch destination"):
        validate_historical_outputs(path)


def test_historical_outputs_allow_scratch_path():
    validate_historical_outputs("results/historical_replay/old.png")
