"""Configuration loading and project-relative path resolution."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


def load_config(path: str | Path) -> tuple[dict[str, Any], Path]:
    config_path = Path(path).expanduser().resolve()
    if not config_path.is_file():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    with config_path.open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    if not isinstance(config, dict):
        raise ValueError("Configuration root must be a mapping.")
    project_root = config_path.parent.parent
    return config, project_root


def resolve_path(project_root: Path, configured_path: str) -> Path:
    path = Path(configured_path).expanduser()
    return path.resolve() if path.is_absolute() else (project_root / path).resolve()


def configured_path(config: dict[str, Any], project_root: Path, key: str) -> Path:
    try:
        raw = config["paths"][key]
    except KeyError as exc:
        raise KeyError(f"Missing configuration path: paths.{key}") from exc
    return resolve_path(project_root, str(raw))
