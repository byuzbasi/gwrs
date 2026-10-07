"""Local foreground controller for the closure-corrected ACS study-v3."""

from __future__ import annotations

import importlib.util
import os
from pathlib import Path
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
APP = HERE.parent
V1_RUN = (APP.parent / "study-v1/code/run.py").resolve()

import workflow as study_v3_workflow

spec = importlib.util.spec_from_file_location("gwrs_acs_study_v1_run_for_v3", V1_RUN)
if spec is None or spec.loader is None:
    raise RuntimeError("Cannot load the hash-bound study-v1 controller")
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)

base.APP = APP
base.atomic_json = study_v3_workflow.atomic_json
base.environment = study_v3_workflow.environment
base.implementation_signature = study_v3_workflow.implementation_signature
base.read_json = study_v3_workflow.read_json
base.require = study_v3_workflow.require
base.sha256 = study_v3_workflow.sha256
base.verify = study_v3_workflow.verify


def local_full_gate():
    study_v3_workflow.require(
        os.environ.get("GWRS_STUDY_LOCAL_FULL") == "YES",
        "The approved local full-run launcher is required",
    )


base.require_calibration_approval = local_full_gate


def plan_priority(row):
    weights = {
        "local-en": 90, "local-lasso": 80, "local-scad": 70,
        "local-mcp": 60, "local-ridge": 30, "local-gwr": 20,
        "global-en": 15, "global-lasso": 14, "global-scad": 13,
        "global-mcp": 12, "global-ridge": 2, "global-ols": 1,
    }
    candidate = row.get("candidate_id", "") or ""
    q = 0
    if "-q" in candidate and "-qglobal" not in candidate:
        try:
            q = int(candidate.split("-q", 1)[1].split("-", 1)[0])
        except ValueError:
            q = 0
    return (-weights.get(row["method_key"], 0), -q, row["id"])


base.plan_priority = plan_priority
original_snapshot = base.Progress.snapshot


def snapshot(self):
    value = original_snapshot(self)
    value["schema"] = "gwrs-usa-counties-acs2024-study-v3-progress-v1"
    value["work_unit"] = "one study-v3 signature-validated model path shard"
    return value


base.Progress.snapshot = snapshot


if __name__ == "__main__":
    base.main()
