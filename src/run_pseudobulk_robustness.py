from __future__ import annotations

import argparse
import shutil
from pathlib import Path, PurePosixPath
from typing import Any

from .run_edger import run_command
from .utils import ensure_dir, load_config


def _container_path(path: Path, repository: Path) -> PurePosixPath:
    absolute = path.resolve()
    try:
        relative = absolute.relative_to(repository.resolve())
    except ValueError as error:
        raise ValueError(
            f"Docker runtime path must be inside the repository: {absolute}"
        ) from error
    return PurePosixPath("/work", *relative.parts)


def _required_sections(cfg: dict[str, Any]) -> tuple[dict[str, Any], ...]:
    de = cfg.get("pseudobulk_de", {})
    robustness = cfg.get("pseudobulk_robustness", {})
    runtime = de.get("runtime", {}) if isinstance(de, dict) else {}
    if not isinstance(de, dict) or not de:
        raise KeyError("pseudobulk_de configuration is required.")
    if not isinstance(robustness, dict) or not robustness:
        raise KeyError("pseudobulk_robustness configuration is required.")
    if not isinstance(runtime, dict) or not runtime:
        raise KeyError("pseudobulk_de.runtime configuration is required.")
    return de, robustness, runtime


def _script_arguments(
    paths: dict[str, Path | PurePosixPath],
    de: dict[str, Any],
    robustness: dict[str, Any],
    runtime: dict[str, Any],
) -> list[str]:
    return [
        "--matrix",
        str(paths["matrix"]),
        "--profiles",
        str(paths["profiles"]),
        "--genes",
        str(paths["genes"]),
        "--primary-results",
        str(paths["primary_results"]),
        "--sensitivity-out",
        str(paths["sensitivity_out"]),
        "--diagnostics-out",
        str(paths["diagnostics_out"]),
        "--manifest-out",
        str(paths["manifest_out"]),
        "--primary-formula",
        str(de["primary_formula"]),
        "--min-replicates",
        str(int(de["min_replicates_per_celltype"])),
        "--min-age-span",
        str(float(de["min_age_span_years"])),
        "--min-residual-df",
        str(int(de["min_residual_df"])),
        "--candidate-fdr",
        str(float(robustness["candidate_global_fdr_threshold"])),
        "--high-cell-min",
        str(int(robustness["high_cell_min_cells"])),
        "--run-leave-one-batch-out",
        str(bool(robustness["run_leave_one_batch_out"])).lower(),
        "--expected-r-version",
        str(runtime["expected_r_version"]),
        "--expected-bioconductor-version",
        str(runtime["expected_bioconductor_version"]),
        "--expected-edger-version",
        str(runtime["expected_edger_version"]),
    ]


def build_command(
    cfg: dict[str, Any],
    paths: dict[str, Path],
    repository: Path,
) -> list[str]:
    de, robustness, runtime = _required_sections(cfg)
    script = repository / "src" / "pseudobulk_robustness.R"
    if not script.is_file():
        raise FileNotFoundError(f"Robustness script not found: {script}")

    mode = str(runtime.get("mode", "")).lower()
    if mode == "native":
        executable = shutil.which(str(runtime.get("rscript", "Rscript")))
        if executable is None:
            raise FileNotFoundError("Configured Rscript executable was not found.")
        native_paths = {key: value.resolve() for key, value in paths.items()}
        return [
            executable,
            str(script.resolve()),
            *_script_arguments(native_paths, de, robustness, runtime),
        ]

    if mode != "docker":
        raise ValueError("pseudobulk_de.runtime.mode must be docker or native.")
    if shutil.which("docker") is None:
        raise FileNotFoundError("Docker executable was not found.")
    image = str(runtime.get("docker_image", "")).strip()
    context = str(runtime.get("docker_context", "")).strip()
    if not image or not context:
        raise KeyError("Docker context and image are required.")
    container_paths = {key: _container_path(value, repository) for key, value in paths.items()}
    return [
        "docker",
        "--context",
        context,
        "run",
        "--rm",
        "--mount",
        f"type=bind,source={repository.resolve()},target=/work",
        "--workdir",
        "/work",
        image,
        "Rscript",
        "/work/src/pseudobulk_robustness.R",
        *_script_arguments(container_paths, de, robustness, runtime),
    ]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run approved pseudobulk robustness sensitivity models."
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--profiles", required=True)
    parser.add_argument("--genes", required=True)
    parser.add_argument("--primary-results", required=True)
    parser.add_argument("--sensitivity-out", required=True)
    parser.add_argument("--diagnostics-out", required=True)
    parser.add_argument("--manifest-out", required=True)
    parser.add_argument("--log-out", required=True)
    args = parser.parse_args()

    repository = Path.cwd().resolve()
    cfg = load_config(args.config)
    paths = {
        "matrix": Path(args.matrix),
        "profiles": Path(args.profiles),
        "genes": Path(args.genes),
        "primary_results": Path(args.primary_results),
        "sensitivity_out": Path(args.sensitivity_out),
        "diagnostics_out": Path(args.diagnostics_out),
        "manifest_out": Path(args.manifest_out),
    }
    for key in ["matrix", "profiles", "genes", "primary_results"]:
        if not paths[key].is_file():
            raise FileNotFoundError(f"Required robustness input not found: {paths[key]}")
    for key in ["sensitivity_out", "diagnostics_out", "manifest_out"]:
        ensure_dir(paths[key].parent)

    _, robustness, runtime = _required_sections(cfg)
    timeout_seconds = int(robustness.get("timeout_seconds", runtime.get("timeout_seconds", 14400)))
    if timeout_seconds <= 0:
        raise ValueError("The pseudobulk robustness timeout must be positive.")
    run_command(
        build_command(cfg, paths, repository),
        repository=repository,
        log_path=Path(args.log_out),
        timeout_seconds=timeout_seconds,
    )


if __name__ == "__main__":
    main()
