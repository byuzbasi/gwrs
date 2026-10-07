"""Foreground resumable controller for the ACS 2024 county study-v1."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import signal
import statistics
import subprocess
import sys
import tempfile
import time
from typing import Optional, Union

sys.dont_write_bytecode = True
from workflow import APP, atomic_json, environment, implementation_signature, read_json, require, sha256, verify


def verify_for_execution():
    if (APP / "VALIDATED.json").exists():
        verify()
    else:
        require(
            os.environ.get("GWRS_STUDY_DEVELOPMENT") == "YES",
            "Only an explicitly labelled local development smoke may use an unfrozen release",
        )
        verify(require_frozen_release=False)


def safe_output(path: Path) -> Path:
    path = Path(path)
    require(path.is_absolute(), "Output path must be absolute")
    require(path.resolve() == path, "Output path must be canonical and contain no symlink")
    require(path.is_relative_to(APP), "Output must stay inside the versioned study-v1 directory")
    current = path
    while current != APP.parent:
        if current.exists():
            require(not current.is_symlink(), f"Symlinked output component refused: {current}")
        if current == APP:
            break
        current = current.parent
    return path


def worker_command(
    action: str,
    output: Path,
    mode: str,
    subject: str = "-",
    event: Union[Path, str] = "-",
):
    return [
        "Rscript",
        "--vanilla",
        "code/worker.R",
        action,
        str(APP),
        str(output),
        mode,
        str(subject),
        implementation_signature(),
        str(event),
    ]


def atomic_text(path: Path, text: str, replace: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=".publish-", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as stream:
            stream.write(text)
            if text and not text.endswith("\n"):
                stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        if replace:
            os.replace(temporary, path)
        else:
            os.link(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


class ControllerLock:
    def __init__(self, output: Path):
        self.path = output / "controller.lock"
        self.stream = None

    def __enter__(self):
        self.stream = self.path.open("a+")
        try:
            fcntl.flock(self.stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError("Another controller holds the run lock") from error
        self.stream.seek(0)
        self.stream.truncate()
        json.dump(
            {
                "pid": os.getpid(),
                "host": os.uname().nodename,
                "slurm_job_id": os.environ.get("SLURM_JOB_ID", ""),
                "started_utc": datetime.now(timezone.utc).isoformat(),
            },
            self.stream,
            sort_keys=True,
        )
        self.stream.write("\n")
        self.stream.flush()
        os.fsync(self.stream.fileno())
        return self

    def __exit__(self, exc_type, exc, traceback):
        if self.stream is not None:
            fcntl.flock(self.stream.fileno(), fcntl.LOCK_UN)
            self.stream.close()


def read_plan(output: Path):
    path = output / "numeric/task-plan.csv"
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    require(rows and len({row["id"] for row in rows}) == len(rows), "Invalid task plan")
    return rows


def run_logged(command, log: Path, env: dict, check: bool = True):
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("x") as stream:
        result = subprocess.run(command, cwd=APP, env=env, stdout=stream, stderr=subprocess.STDOUT)
    if check and result.returncode != 0:
        raise RuntimeError(f"Command failed ({result.returncode}); inspect {log}")
    return result


def scan(output: Path, mode: str, stage: str, invocation: Path, env: dict):
    event = invocation / f"scan-{stage}-{time.time_ns()}.json"
    log = invocation / f"scan-{stage}-{time.time_ns()}.log"
    result = run_logged(worker_command("scan", output, mode, stage, event), log, env, check=False)
    require(event.exists(), f"Scan produced no receipt; inspect {log}")
    value = read_json(event)
    require(value["scientific_signature"], "Scan omitted signature")
    if value["invalid"]:
        detail = "; ".join(f"{row['id']}: {row['message']}" for row in value["invalid"][:3])
        raise RuntimeError("Invalid existing shard; no overwrite attempted: " + detail)
    require(result.returncode == 0, f"Scan failed; inspect {log}")
    return value


def plan_priority(row):
    method = row["method_key"]
    weights = {
        "local-en": 90,
        "local-lasso": 80,
        "local-scad": 70,
        "local-mcp": 60,
        "local-ridge": 30,
        "local-gwr": 20,
        "global-en": 15,
        "global-lasso": 14,
        "global-scad": 13,
        "global-mcp": 12,
        "global-ridge": 2,
        "global-ols": 1,
    }
    candidate = row.get("candidate_id", "") or ""
    k = 0
    if "-k" in candidate:
        try:
            k = int(candidate.split("-k", 1)[1].split("-", 1)[0])
        except ValueError:
            k = 0
    return (-weights.get(method, 0), -k, row["id"])


class Progress:
    fields = [
        "heartbeat_utc",
        "run_id",
        "slurm_job_id",
        "state",
        "phase",
        "total",
        "completed",
        "running",
        "failed",
        "pending",
        "percent_complete",
        "elapsed_seconds",
        "throughput_tasks_per_second",
        "eta_seconds",
        "estimated_completion_utc",
        "remaining_slurm_wall_seconds",
        "last_completed_task",
    ]

    def __init__(self, output: Path, plan, workers: int, heartbeat: float):
        self.output = output
        self.plan = {row["id"]: row for row in plan}
        self.total = len(plan)
        self.workers = workers
        self.heartbeat = heartbeat
        self.started = time.monotonic()
        self.started_utc = time.time()
        self.completed = set()
        self.running = set()
        self.failures = {}
        self.durations = {}
        self.current_fresh = set()
        self.phase = "initializing"
        self.state = "running"
        self.last_completed = None
        self.history = output / "progress.tsv"

    def seed(self, scan_value):
        for row in scan_value["valid"]:
            self.completed.add(row["id"])
            self.durations[row["id"]] = float(row["elapsed_seconds"])

    def done(self, event):
        task = event["id"]
        self.running.discard(task)
        self.failures.pop(task, None)
        self.completed.add(task)
        self.current_fresh.add(task)
        self.durations[task] = float(event["controller_elapsed_seconds"])
        self.last_completed = task

    def failed(self, task: str, message: str):
        self.running.discard(task)
        self.failures[task] = message

    def estimate(self):
        stage_ids = [task for task, row in self.plan.items() if row["stage"] == self.phase]
        if not stage_ids:
            return None, None
        remaining = [task for task in stage_ids if task not in self.completed]
        observed = [task for task in stage_ids if task in self.current_fresh]
        if len(observed) < 3 or not remaining:
            return (0.0 if not remaining else None), (datetime.now(timezone.utc).isoformat() if not remaining else None)
        by_class = {}
        for task in observed:
            by_class.setdefault(self.plan[task]["cost_class"], []).append(self.durations[task])
        if any(self.plan[task]["cost_class"] not in by_class for task in remaining):
            return None, None
        work = sum(
            statistics.median(by_class[self.plan[task]["cost_class"]])
            for task in remaining
        )
        eta = work / max(1, self.workers)
        completion = datetime.fromtimestamp(time.time() + eta, timezone.utc).isoformat()
        return eta, completion

    def snapshot(self):
        elapsed = max(0.0, time.monotonic() - self.started)
        eta, completion = self.estimate()
        if self.state == "complete":
            eta, completion = 0.0, datetime.now(timezone.utc).isoformat()
        failed = len(self.failures)
        pending = max(0, self.total - len(self.completed) - len(self.running) - failed)
        wall_limit = os.environ.get("GWRS_SLURM_WALLTIME_SECONDS", "")
        wall_remaining = None
        if os.environ.get("SLURM_JOB_END_TIME", "").isdigit():
            wall_remaining = max(0.0, float(os.environ["SLURM_JOB_END_TIME"]) - time.time())
        elif wall_limit:
            wall_remaining = max(0.0, float(wall_limit) - elapsed)
        return {
            "schema": "gwrs-usa-counties-acs2024-study-v1-progress-v1",
            "heartbeat_utc": datetime.now(timezone.utc).isoformat(),
            "heartbeat_stale_after_seconds": max(120.0, 3 * self.heartbeat),
            "run_id": self.output.name,
            "slurm_job_id": os.environ.get("SLURM_JOB_ID", "local"),
            "state": self.state,
            "phase": self.phase,
            "work_unit": "one signature-validated model path shard",
            "total": self.total,
            "completed": len(self.completed),
            "resumed_completed": len(self.completed - self.current_fresh),
            "current_invocation_completed": len(self.current_fresh),
            "running": len(self.running),
            "failed": failed,
            "pending": pending,
            "percent_complete": 100.0 * len(self.completed) / self.total,
            "elapsed_seconds": elapsed,
            "throughput_tasks_per_second": len(self.current_fresh) / max(elapsed, 1e-9),
            "eta_seconds": eta,
            "estimated_completion_utc": completion,
            "eta_scope": "current computational stage only; selection, later stages and final validation excluded",
            "remaining_slurm_wall_seconds": wall_remaining,
            "remaining_slurm_wall_scope": "allocation deadline, not computational ETA",
            "workers_requested": self.workers,
            "native_threads_per_worker": 1,
            "last_completed_task": self.last_completed,
            "failure_tasks": sorted(self.failures),
        }

    def publish(self):
        value = self.snapshot()
        atomic_json(self.output / "progress.json", value, replace=True)
        write_header = not self.history.exists()
        with self.history.open("a", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=self.fields, delimiter="\t", extrasaction="ignore")
            if write_header:
                writer.writeheader()
            writer.writerow(value)
            stream.flush()
            os.fsync(stream.fileno())
        eta = "estimating" if value["eta_seconds"] is None else f"{value['eta_seconds']:.1f}s"
        print(
            f"[{value['heartbeat_utc']}] run={value['run_id']} job={value['slurm_job_id']} "
            f"phase={value['phase']} state={value['state']} "
            f"done={value['completed']}/{value['total']} running={value['running']} "
            f"failed={value['failed']} pending={value['pending']} "
            f"percent={value['percent_complete']:.2f} ETA={eta} "
            f"last={value['last_completed_task'] or '-'}",
            flush=True,
        )


class Interrupted(Exception):
    pass


def terminate_active(active):
    for item in active.values():
        process = item["process"]
        if process.poll() is None:
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
    deadline = time.monotonic() + 10
    while time.monotonic() < deadline and any(item["process"].poll() is None for item in active.values()):
        time.sleep(0.1)
    for item in active.values():
        process = item["process"]
        if process.poll() is None:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        process.wait()
        item["stream"].close()


def run_stage(output, mode, stage, rows, invocation, env, progress, task_budget):
    snapshot = scan(output, mode, stage, invocation, env)
    valid = {row["id"] for row in snapshot["valid"]}
    queue = sorted([row for row in rows if row["id"] not in valid], key=plan_priority)
    active = {}
    failures = []
    launched = 0
    interrupted = False

    def signal_handler(signum, frame):
        nonlocal interrupted
        interrupted = True

    old_handlers = {sig: signal.signal(sig, signal_handler) for sig in (signal.SIGTERM, signal.SIGINT)}
    next_heartbeat = time.monotonic()
    try:
        while queue or active:
            if interrupted:
                raise Interrupted("Scheduler/user interruption")
            while queue and len(active) < progress.workers and launched < task_budget:
                row = queue.pop(0)
                task = row["id"]
                event = invocation / "events" / f"{task}.json"
                log = invocation / "tasks" / f"{task}.log"
                event.parent.mkdir(parents=True, exist_ok=True)
                log.parent.mkdir(parents=True, exist_ok=True)
                stream = log.open("x")
                process = subprocess.Popen(
                    worker_command("task", output, mode, task, event),
                    cwd=APP,
                    env=env,
                    stdout=stream,
                    stderr=subprocess.STDOUT,
                    start_new_session=True,
                )
                active[task] = {"process": process, "stream": stream, "event": event, "log": log, "started": time.monotonic()}
                progress.running.add(task)
                launched += 1
            finished = [task for task, item in active.items() if item["process"].poll() is not None]
            for task in finished:
                item = active.pop(task)
                item["stream"].close()
                code = item["process"].returncode
                if code == 0 and item["event"].exists():
                    event = read_json(item["event"])
                    require(event["id"] == task, "Task event ID mismatch")
                    event["controller_elapsed_seconds"] = time.monotonic() - item["started"]
                    progress.done(event)
                else:
                    message = f"exit={code}; log={item['log']}"
                    progress.failed(task, message)
                    failures.append((task, message))
            if time.monotonic() >= next_heartbeat or finished:
                progress.publish()
                next_heartbeat = time.monotonic() + progress.heartbeat
            if launched >= task_budget and not active:
                break
            time.sleep(0.2)
    except BaseException:
        terminate_active(active)
        raise
    finally:
        for sig, handler in old_handlers.items():
            signal.signal(sig, handler)
    verified = scan(output, mode, stage, invocation, env)
    missing = set(verified["missing"])
    if failures:
        raise RuntimeError(f"{len(failures)} {stage} task process(es) failed; valid shards retained")
    incomplete = bool(missing)
    return launched, incomplete


def invocation_manifest(invocation: Path):
    paths = sorted(path for path in invocation.rglob("*") if path.is_file() and path.name != "manifest-sha256.json")
    value = {
        "schema": "gwrs-usa-counties-acs2024-study-v1-invocation-manifest-v1",
        "files": [
            {"file": str(path.relative_to(invocation)), "bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in paths
        ],
    }
    atomic_json(invocation / "manifest-sha256.json", value)


def call_selection(output, mode, subject, label, invocation, env):
    event = invocation / "selection" / f"{label}.json"
    log = invocation / "selection" / f"{label}.log"
    event.parent.mkdir(parents=True, exist_ok=True)
    run_logged(worker_command("select", output, mode, str(subject), event), log, env)
    require(event.exists(), "Missing selection receipt")
    value = read_json(event)
    require(value["scope"] in ("outer", "full"), "Invalid selection receipt")


def verify_run(
    output: Path,
    mode: str,
    invocation: Optional[Path] = None,
    env: Optional[dict] = None,
):
    output = safe_output(output)
    env = environment() if env is None else env
    created_invocation = invocation is None
    if invocation is None:
        invocation = output / "invocations" / f"verify-{time.time_ns()}"
        invocation.mkdir(parents=True)
    log = invocation / "read-only-verification.log"
    try:
        result = run_logged(worker_command("verify", output, mode), log, env)
        require(result.returncode == 0, "Verification failed")
        print(f"Verified read-only: {output}")
    finally:
        if created_invocation:
            invocation_manifest(invocation)


def preflight(output: Path, mode: str, workers: int):
    output = safe_output(output)
    env = environment()
    env["GWRS_STUDY_WORKERS"] = str(workers)
    subprocess.run(worker_command("preflight", output, mode), cwd=APP, env=env, check=True)


def calibration_approval_path():
    return APP / "runtime/calibration-approval-v1.json"


def approve_calibration(output: Path, job_id: str):
    require(job_id.isdigit(), "A numeric completed SLURM job ID is required")
    verify()
    output = safe_output(output)
    verify_run(output, "calibration")
    completion = output / "COMPLETED"
    approval = calibration_approval_path()
    require(not approval.exists(), "Calibration approval already exists")
    atomic_json(
        approval,
        {
            "schema": "gwrs-usa-counties-acs2024-study-v1-calibration-approval-v1",
            "approved_utc": datetime.now(timezone.utc).isoformat(),
            "author": "Bahadir Yuzbasi",
            "slurm_job_id": job_id,
            "calibration_output": str(output.relative_to(APP)),
            "completion_sha256": sha256(completion),
            "implementation_signature": implementation_signature(),
            "meaning": "Author-reviewed Linux calibration gate; not a scientific result",
        },
    )
    print(approval)


def require_calibration_approval():
    path = Path(os.environ.get("GWRS_STUDY_CALIBRATION_APPROVAL", calibration_approval_path()))
    require(path.resolve() == path and path.is_relative_to(APP), "Unsafe calibration approval path")
    value = read_json(path)
    require(
        value["schema"] == "gwrs-usa-counties-acs2024-study-v1-calibration-approval-v1",
        "Wrong calibration approval",
    )
    require(value["implementation_signature"] == implementation_signature(), "Stale calibration approval")
    output = APP / value["calibration_output"]
    require(sha256(output / "COMPLETED") == value["completion_sha256"], "Calibration changed after approval")


def execute(output: Path, mode: str, workers: int, heartbeat: float, author_run: bool,
            max_new_tasks: Optional[int], allow_incomplete: bool):
    verify_for_execution()
    output = safe_output(output)
    require(workers >= 1, "workers must be positive")
    if mode in ("calibration", "full"):
        require(author_run, "Empirical execution requires --author-run")
    if mode == "full":
        require(os.environ.get("GWRS_STUDY_PRODUCTION") == "YES", "Production environment gate is closed")
        require_calibration_approval()
    output.mkdir(parents=True, exist_ok=True)
    with ControllerLock(output):
        if (output / "COMPLETED").exists():
            verify_run(output, mode)
            print("Completed run is immutable; no model was refitted.")
            return
        env = environment()
        env["GWRS_STUDY_WORKERS"] = str(workers)
        invocation = output / "invocations" / f"attempt-{time.strftime('%Y%m%dT%H%M%SZ', time.gmtime())}-{time.time_ns()}"
        invocation.mkdir(parents=True)
        atomic_json(
            invocation / "controller.json",
            {
                "schema": "gwrs-usa-counties-acs2024-study-v1-controller-v1",
                "python": sys.version,
                "executable": sys.executable,
                "mode": mode,
                "workers": workers,
                "implementation_signature": implementation_signature(),
                "slurm_job_id": os.environ.get("SLURM_JOB_ID", ""),
            },
        )
        try:
            run_logged(worker_command("preflight", output, mode), invocation / "preflight.log", env)
            run_logged(worker_command("initialize", output, mode), invocation / "initialize.log", env)
            plan = read_plan(output)
            all_scan = scan(output, mode, "all", invocation, env)
            progress = Progress(output, plan, workers, heartbeat)
            progress.seed(all_scan)
            progress.publish()
            budget = math.inf if max_new_tasks is None else max_new_tasks
            require(budget > 0, "max-new-tasks must be positive")
            used = 0
            incomplete = False
            if mode == "calibration":
                progress.phase = "calibration"
                launched, incomplete = run_stage(
                    output, mode, "calibration", plan, invocation, env, progress, budget
                )
                used += launched
            else:
                policy = read_json(APP / "config/policy.json")
                folds = (
                    policy["smoke"]["outer_folds"]
                    if mode == "smoke"
                    else policy["full"]["outer_folds"]
                )

                def available_budget():
                    return math.inf if math.isinf(budget) else max(0, budget - used)

                progress.phase = "nested_inner"
                nested_rows = [row for row in plan if row["stage"] == "nested_inner"]
                launched, incomplete = run_stage(
                    output, mode, "nested_inner", nested_rows,
                    invocation, env, progress, available_budget(),
                )
                used += launched

                if not incomplete:
                    progress.phase = "nested_selection"
                    progress.publish()
                    for outer in folds:
                        call_selection(
                            output, mode, outer, f"outer-{outer}", invocation, env
                        )

                    progress.phase = "outer_refit"
                    outer_rows = [row for row in plan if row["stage"] == "outer_refit"]
                    launched, incomplete = run_stage(
                        output, mode, "outer_refit", outer_rows,
                        invocation, env, progress, available_budget(),
                    )
                    used += launched

                if not incomplete:
                    progress.phase = "full_tune"
                    full_tune_rows = [row for row in plan if row["stage"] == "full_tune"]
                    launched, incomplete = run_stage(
                        output, mode, "full_tune", full_tune_rows,
                        invocation, env, progress, available_budget(),
                    )
                    used += launched

                if not incomplete:
                    progress.phase = "full_selection"
                    progress.publish()
                    call_selection(output, mode, "full", "full", invocation, env)

                    progress.phase = "full_refit"
                    full_refit_rows = [row for row in plan if row["stage"] == "full_refit"]
                    launched, incomplete = run_stage(
                        output, mode, "full_refit", full_refit_rows,
                        invocation, env, progress, available_budget(),
                    )
                    used += launched
            if incomplete:
                progress.state = "incomplete"
                progress.phase = "waiting_for_resume"
                progress.publish()
                if not allow_incomplete:
                    raise RuntimeError("Validated task budget reached; submit the same command to resume")
                print("Intentional incomplete stop; valid shards retained for resume.")
                return
            progress.phase = "aggregation"
            progress.publish()
            final_event = invocation / "finalize.json"
            run_logged(worker_command("finalize", output, mode, "-", final_event), invocation / "finalize.log", env)
            require(final_event.exists(), "Finalizer produced no receipt")
            verify_run(output, mode, invocation, env)
            final_scan = scan(output, mode, "all", invocation, env)
            progress.seed(final_scan)
            progress.state = "complete"
            progress.phase = "validated_completion"
            progress.publish()
            print(f"Run closed after task, OOF and checksum validation: {output}")
        except BaseException:
            if "progress" in locals():
                progress.state = "failed"
                progress.phase = "failure_requires_review"
                progress.publish()
            raise
        finally:
            invocation_manifest(invocation)


def status(output: Path):
    output = safe_output(output)
    value = read_json(output / "progress.json")
    heartbeat = datetime.fromisoformat(value["heartbeat_utc"])
    age = (datetime.now(timezone.utc) - heartbeat).total_seconds()
    value["heartbeat_age_seconds"] = age
    value["stale_running_heartbeat"] = value["state"] == "running" and age > value["heartbeat_stale_after_seconds"]
    print(json.dumps(value, indent=2, sort_keys=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("preflight", "run", "verify", "status", "approve-calibration"))
    parser.add_argument("--mode", choices=("smoke", "calibration", "full"), default="smoke")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--heartbeat", type=float, default=30.0)
    parser.add_argument("--author-run", action="store_true")
    parser.add_argument("--max-new-tasks", type=int)
    parser.add_argument("--allow-incomplete", action="store_true")
    parser.add_argument("--job-id", default="")
    args = parser.parse_args()
    require(args.heartbeat >= 1.0, "heartbeat must be at least one second")
    if args.action == "preflight":
        verify_for_execution()
        preflight(args.output, args.mode, args.workers)
    elif args.action == "run":
        execute(
            args.output, args.mode, args.workers, args.heartbeat,
            args.author_run, args.max_new_tasks, args.allow_incomplete,
        )
    elif args.action == "verify":
        verify_for_execution()
        verify_run(args.output, args.mode)
    elif args.action == "status":
        status(args.output)
    else:
        require(args.mode == "calibration", "Approval applies only to calibration")
        approve_calibration(args.output, args.job_id)


if __name__ == "__main__":
    main()
