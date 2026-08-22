"""Automated public-release checks; human disclosure review remains required."""

from __future__ import annotations

import csv
import re
from pathlib import Path

FORBIDDEN_IDENTIFIER_COLUMNS = {
    "rid",
    "ptid",
    "ptid_std",
    "med_id",
    "subject_id",
    "participant_id",
    "image_id",
    "scan_id",
    "series_uid",
}
CONTROLLED_DIRECTORIES = {
    Path("data/raw"),
    Path("data/derived"),
}
ALLOWED_CONTROLLED_FILES = {".gitkeep", "README.md"}
TEXT_SUFFIXES = {".py", ".r", ".md", ".yml", ".yaml", ".toml", ".cff", ".txt"}
TABULAR_SUFFIXES = {".csv", ".tsv"}
OPAQUE_TABULAR_SUFFIXES = {".xlsx", ".xls", ".parquet", ".feather"}


def _inside(path: Path, directory: Path) -> bool:
    try:
        path.relative_to(directory)
        return True
    except ValueError:
        return False


def audit_release_tree(root: Path) -> list[str]:
    root = root.resolve()
    findings: list[str] = []
    for path in sorted(root.rglob("*")):
        if not path.is_file() or ".git" in path.parts or "__pycache__" in path.parts:
            continue
        relative = path.relative_to(root)

        for controlled in CONTROLLED_DIRECTORIES:
            controlled_absolute = root / controlled
            if _inside(path, controlled_absolute) and path.name not in ALLOWED_CONTROLLED_FILES:
                findings.append(f"Controlled-data directory contains a file: {relative}")

        if path.suffix.lower() in TEXT_SUFFIXES:
            text = path.read_text(encoding="utf-8", errors="ignore")
            if re.search(r"(?<![A-Za-z])/(Users|home)/[^\\s\"']+", text):
                findings.append(f"Hard-coded home-directory path: {relative}")
            if re.search(
                r"(?i)(api[_-]?key|secret|token|password)\\s*[:=]\\s*[\"'][^\"']+[\"']",
                text,
            ):
                findings.append(f"Possible embedded secret: {relative}")

        if path.suffix.lower() in TABULAR_SUFFIXES:
            try:
                delimiter = "," if path.suffix.lower() == ".csv" else "\t"
                with path.open(encoding="utf-8-sig", errors="ignore", newline="") as handle:
                    columns = next(csv.reader(handle, delimiter=delimiter), [])
            except Exception as error:  # noqa: BLE001
                findings.append(f"Could not inspect tabular file {relative}: {error}")
                continue
            normalised = {str(column).strip().lower() for column in columns}
            exposed = sorted(normalised & FORBIDDEN_IDENTIFIER_COLUMNS)
            if exposed:
                findings.append(
                    f"Identifier columns {exposed} in public tabular file: {relative}"
                )
        if path.suffix.lower() in OPAQUE_TABULAR_SUFFIXES:
            findings.append(
                f"Opaque tabular format requires manual conversion/review: {relative}"
            )
    return findings
