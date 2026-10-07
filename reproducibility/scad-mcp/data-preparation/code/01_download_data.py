#!/usr/bin/env python3
"""Acquire ACS county values and boundaries without exposing the API key."""

from pathlib import Path
import os
import sys


HERE = Path(__file__).resolve()





key = os.environ.get("CENSUS_API_KEY", "")
if not key or key.strip() != key or any(character.isspace() for character in key):
    raise SystemExit(
        "CENSUS_API_KEY is required and must contain no whitespace; "
        "the key was not printed or written."
    )

ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parent))

from download_sources import main  # noqa: E402


if __name__ == "__main__":
    config = HERE.parents[1] / "config" / "data-sources-v1.json"
    raise SystemExit(main(["--strict", "--config", str(config)]))

