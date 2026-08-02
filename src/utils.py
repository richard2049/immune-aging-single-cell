from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


def _deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    merged = dict(base)
    for key, value in override.items():
        if key in merged and isinstance(merged[key], dict) and isinstance(value, dict):
            merged[key] = _deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


def load_config(path: str | Path) -> dict:
    config_path = Path(path)
    with config_path.open("r", encoding="utf-8") as handle:
        current = yaml.safe_load(handle) or {}
    if not isinstance(current, dict):
        raise TypeError(f"Configuration must be a mapping: {config_path}")

    extends = current.pop("extends", None)
    if not extends:
        return current
    base_path = Path(str(extends))
    if not base_path.is_absolute():
        cwd_candidate = base_path
        relative_candidate = config_path.parent / base_path
        base_path = cwd_candidate if cwd_candidate.exists() else relative_candidate
    if base_path.resolve() == config_path.resolve():
        raise ValueError(f"Configuration cannot extend itself: {config_path}")
    return _deep_merge(load_config(base_path), current)


def ensure_dir(path: str | Path) -> None:
    Path(path).mkdir(parents=True, exist_ok=True)
