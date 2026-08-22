#!/usr/bin/env python3
"""Auxiliary Python release audit for controlled data and identifiers."""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pved_ptau217.privacy import audit_release_tree  # noqa: E402


def main() -> int:
    findings = audit_release_tree(ROOT)
    if findings:
        print("RELEASE AUDIT FAILED")
        for finding in findings:
            print(f"- {finding}")
        return 1
    print("RELEASE AUDIT PASSED: no automated finding.")
    print("Human disclosure and data-governance review is still required.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
