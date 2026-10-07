#!/usr/bin/env python3
"""Download only the official ACS 2024 metadata through the shared utility."""

from pathlib import Path
import sys


HERE = Path(__file__).resolve()





ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parent))

from download_sources import main  # noqa: E402


if __name__ == "__main__":
    config = HERE.parents[1] / "config" / "metadata-sources-v1.json"
    raise SystemExit(main(["--strict", "--config", str(config)]))

