from __future__ import annotations

import argparse
import shlex
import shutil
import subprocess
import threading
import time
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Any, TextIO

from .utils import ensure_dir, load_config


def _required_runtime_config(cfg: dict[str, Any]) -> dict[str, Any]:
    section = cfg.get("pseudobulk_de", {}).get("runtime", {})
    if not isinstance(section, dict) or not section:
        raise KeyError("pseudobulk_de.runtime configuration is required.")
    for key in [
        "mode",
        "expected_r_version",
        "expected_bioconductor_version",
        "expected_edger_version",
    ]:
        if key not in section:
            raise KeyError(f"Missing pseudobulk_de.runtime.{key}.")
    return section


def _contract_arguments(cfg: dict[str, Any]) -> list[str]:
    contract = cfg.get("pseudobulk_de", {})
    return [
        "--primary-formula",
        str(contract.get("primary_formula", "~ sex + batch + age_decade")),
        "--min-replicates",
        str(int(contract.get("min_replicates_per_celltype", 12))),
        "--min-age-span",
        str(float(contract.get("min_age_span_years", 12))),
        "--min-residual-df",
        str(int(contract.get("min_residual_df", 5))),
    ]


def _runtime_arguments(runtime: dict[str, Any]) -> list[str]:
    return [
        "--expected-r-version",
        str(runtime["expected_r_version"]),
        "--expected-bioconductor-version",
        str(runtime["expected_bioconductor_version"]),
        "--expected-edger-version",
        str(runtime["expected_edger_version"]),
    ]


def _analysis_arguments(paths: dict[str, Any]) -> list[str]:
    flags = {
        "matrix": "matrix",
        "profiles": "profiles",
        "genes": "genes",
        "combined_out": "combined-out",
        "celltype_dir": "celltype-dir",
        "manifest_out": "manifest-out",
        "diagnostics_out": "diagnostics-out",
        "plot_out": "plot-out",
        "runtime_versions_out": "runtime-versions-out",
        "session_info_out": "session-info-out",
    }
    values: list[str] = []
    for key, flag in flags.items():
        values.extend([f"--{flag}", str(paths[key])])
    return values


def _container_path(path: Path, repository: Path) -> PurePosixPath:
    absolute = path.resolve()
    try:
        relative = absolute.relative_to(repository.resolve())
    except ValueError as error:
        raise ValueError(
            f"Docker runtime path must be inside the repository: {absolute}"
        ) from error
    return PurePosixPath("/work", *relative.parts)


def build_command(
    cfg: dict[str, Any],
    paths: dict[str, Path],
    repository: Path,
) -> list[str]:
    runtime = _required_runtime_config(cfg)
    mode = str(runtime["mode"]).lower()
    script = repository / "src" / "pseudobulk_edger.R"
    if not script.exists():
        raise FileNotFoundError(f"edgeR script not found: {script}")

    if mode == "native":
        configured = str(runtime.get("rscript", "Rscript"))
        executable = shutil.which(configured)
        if executable is None:
            raise FileNotFoundError(
                f"Rscript executable not found: {configured}. "
                "Use the approved Docker runtime or install the pinned native "
                "R/Bioconductor environment."
            )
        runtime_paths = {key: value.resolve() for key, value in paths.items()}
        return (
            [executable, str(script.resolve())]
            + _analysis_arguments(runtime_paths)
            + _contract_arguments(cfg)
            + _runtime_arguments(runtime)
        )

    if mode != "docker":
        raise ValueError("pseudobulk_de.runtime.mode must be docker or native.")
    if shutil.which("docker") is None:
        raise FileNotFoundError(
            "Docker executable not found. Install Docker Desktop or select a "
            "validated native R runtime."
        )
    image = str(runtime.get("docker_image", "")).strip()
    if not image:
        raise KeyError("pseudobulk_de.runtime.docker_image is required.")
    context = str(runtime.get("docker_context", "desktop-linux")).strip()
    if not context:
        raise KeyError("pseudobulk_de.runtime.docker_context is required.")

    runtime_paths = {key: _container_path(value, repository) for key, value in paths.items()}
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
        "/work/src/pseudobulk_edger.R",
        *_analysis_arguments(runtime_paths),
        *_contract_arguments(cfg),
        *_runtime_arguments(runtime),
    ]


def _pump_output(stream: TextIO, log: TextIO) -> None:
    for line in iter(stream.readline, ""):
        print(line, end="", flush=True)
        log.write(line)
        log.flush()


def run_command(
    command: list[str],
    repository: Path,
    log_path: Path,
    timeout_seconds: int,
) -> None:
    ensure_dir(log_path.parent)
    started_at = datetime.now(timezone.utc)
    started_monotonic = time.monotonic()
    with log_path.open("w", encoding="utf-8") as log:
        log.write(f"started_at_utc: {started_at.isoformat()}\n")
        log.write(f"working_directory: {repository.resolve()}\n")
        log.write(f"command: {shlex.join(command)}\n\n")
        log.flush()
        process = subprocess.Popen(
            command,
            cwd=repository,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
            bufsize=1,
        )
        assert process.stdout is not None
        pump = threading.Thread(
            target=_pump_output,
            args=(process.stdout, log),
            daemon=True,
        )
        pump.start()
        try:
            return_code = process.wait(timeout=timeout_seconds)
        except subprocess.TimeoutExpired as error:
            process.terminate()
            try:
                process.wait(timeout=30)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
            pump.join(timeout=5)
            raise TimeoutError(
                f"edgeR exceeded the {timeout_seconds}-second timeout. Inspect {log_path}."
            ) from error
        finally:
            process.stdout.close()
            pump.join(timeout=5)

        elapsed = time.monotonic() - started_monotonic
        log.write(f"\nfinished_at_utc: {datetime.now(timezone.utc).isoformat()}\n")
        log.write(f"return_code: {return_code}\n")
        log.write(f"elapsed_seconds: {elapsed:.3f}\n")
        if return_code != 0:
            raise subprocess.CalledProcessError(return_code, command)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run the approved edgeR pseudobulk model in a pinned Docker or "
            "native R/Bioconductor runtime."
        )
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--profiles", required=True)
    parser.add_argument("--genes", required=True)
    parser.add_argument("--combined-out", required=True)
    parser.add_argument("--celltype-dir", required=True)
    parser.add_argument("--manifest-out", required=True)
    parser.add_argument("--diagnostics-out", required=True)
    parser.add_argument("--plot-out", required=True)
    parser.add_argument("--runtime-versions-out", required=True)
    parser.add_argument("--session-info-out", required=True)
    parser.add_argument("--log-out", required=True)
    args = parser.parse_args()

    repository = Path.cwd().resolve()
    cfg = load_config(args.config)
    paths = {
        "matrix": Path(args.matrix),
        "profiles": Path(args.profiles),
        "genes": Path(args.genes),
        "combined_out": Path(args.combined_out),
        "celltype_dir": Path(args.celltype_dir),
        "manifest_out": Path(args.manifest_out),
        "diagnostics_out": Path(args.diagnostics_out),
        "plot_out": Path(args.plot_out),
        "runtime_versions_out": Path(args.runtime_versions_out),
        "session_info_out": Path(args.session_info_out),
    }
    for key in ["matrix", "profiles", "genes"]:
        if not paths[key].exists():
            raise FileNotFoundError(f"Required edgeR input not found: {paths[key]}")
    for key, path in paths.items():
        if key not in {"matrix", "profiles", "genes", "celltype_dir"}:
            ensure_dir(path.parent)
    ensure_dir(paths["celltype_dir"])

    runtime = _required_runtime_config(cfg)
    timeout_seconds = int(runtime.get("timeout_seconds", 14400))
    if timeout_seconds <= 0:
        raise ValueError("pseudobulk_de.runtime.timeout_seconds must be positive.")
    command = build_command(cfg, paths, repository)
    run_command(
        command,
        repository=repository,
        log_path=Path(args.log_out),
        timeout_seconds=timeout_seconds,
    )


if __name__ == "__main__":
    main()
