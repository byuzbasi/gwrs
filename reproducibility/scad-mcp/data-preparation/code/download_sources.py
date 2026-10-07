#!/usr/bin/env python3
"""Immutable, manifest-driven acquisition for the project applications.

Only Python's standard library is used.  Files are streamed to ``.part``,
validated, hashed with SHA-256, and atomically moved into a dated raw snapshot.
Existing raw files are never overwritten.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import os
import shutil
import ssl
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path
from typing import Any, Iterable


USER_AGENT = "gwrs-real-data-acquisition/1.0 (+local reproducible research)"
CHUNK_BYTES = 1024 * 1024
FREE_SPACE_CHECK_BYTES = 64 * CHUNK_BYTES


def trusted_tls_context() -> ssl.SSLContext:
    """Use Python's configured CA file, then common verified system bundles."""
    candidates = [
        ssl.get_default_verify_paths().cafile,
        "/etc/ssl/cert.pem",
        "/private/etc/ssl/cert.pem",
        "/usr/local/etc/openssl@3/cert.pem",
    ]
    for candidate in candidates:
        if candidate and Path(candidate).is_file():
            return ssl.create_default_context(cafile=candidate)
    return ssl.create_default_context()


TLS_CONTEXT = trusted_tls_context()


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat().replace(
        "+00:00", "Z"
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(CHUNK_BYTES), b""):
            digest.update(chunk)
    return digest.hexdigest()


def bool_text(value: Any) -> str:
    return "true" if bool(value) else "false"


def sanitize(text: str, secrets: Iterable[str]) -> str:
    result = text
    for secret in secrets:
        if secret:
            result = result.replace(secret, "<REDACTED>")
    return result


def next_manifest_version(raw_dir: Path) -> int:
    versions: list[int] = []
    for path in raw_dir.glob("source_manifest_v*.csv"):
        suffix = path.stem.removeprefix("source_manifest_v")
        if suffix.isdigit():
            versions.append(int(suffix))
    return max(versions, default=0) + 1


def validate_download(path: Path, file_format: str, minimum_bytes: int) -> None:
    size = path.stat().st_size
    if size < minimum_bytes:
        raise ValueError(
            f"downloaded file has {size} bytes; expected at least {minimum_bytes}"
        )

    if file_format in {"gzip", "csv_gz", "tsv_gz"}:
        total_uncompressed = 0
        with gzip.open(path, "rb") as stream:
            for chunk in iter(lambda: stream.read(CHUNK_BYTES), b""):
                total_uncompressed += len(chunk)
        if total_uncompressed == 0:
            raise ValueError("gzip archive expands to an empty file")
    elif file_format in {"zip", "xlsx", "docx"}:
        with zipfile.ZipFile(path) as archive:
            if not archive.namelist():
                raise ValueError("zip archive contains no entries")
            damaged = archive.testzip()
            if damaged is not None:
                raise ValueError(f"zip integrity check failed at {damaged}")
    elif file_format == "xls":
        with path.open("rb") as stream:
            signature = stream.read(8)
        if signature != bytes.fromhex("d0cf11e0a1b11ae1"):
            raise ValueError("invalid legacy XLS/OLE signature")
    elif file_format == "json":
        with path.open("r", encoding="utf-8-sig") as stream:
            json.load(stream)
    elif file_format == "html":
        with path.open("rb") as stream:
            prefix = stream.read(4096).lower()
        if b"<html" not in prefix and b"<!doctype html" not in prefix:
            raise ValueError("invalid HTML signature")
    elif file_format in {"csv", "tsv"}:
        with path.open("r", encoding="utf-8-sig", errors="strict") as stream:
            sample = stream.read(65536)
        if not sample or sample.lstrip().lower().startswith(("<html", "<!doctype")):
            raise ValueError("invalid delimited text: empty or HTML response")
        dialect = csv.Sniffer().sniff(sample, delimiters=",;\t")
        first_row = next(csv.reader(sample.splitlines(), dialect=dialect), [])
        if len(first_row) < 2:
            raise ValueError("invalid delimited text: fewer than two columns")


def credentialized_url(source: dict[str, Any]) -> tuple[str | None, list[str]]:
    template = source.get("url", "")
    env_name = source.get("api_key_env")
    if not env_name:
        return template, []
    secret = os.environ.get(env_name, "")
    if not secret:
        return None, []
    placeholder = "{" + env_name + "}"
    return template.replace(placeholder, urllib.parse.quote(secret, safe="")), [secret]


def request(url: str, method: str = "GET") -> urllib.request.Request:
    return urllib.request.Request(
        url,
        method=method,
        headers={
            "User-Agent": USER_AGENT,
            "Accept": "*/*",
            "Accept-Encoding": "identity",
        },
    )


def probe(url: str, timeout_seconds: int) -> dict[str, Any]:
    try:
        with urllib.request.urlopen(
            request(url, method="HEAD"), timeout=timeout_seconds,
            context=TLS_CONTEXT,
        ) as response:
            length = response.headers.get("Content-Length")
            return {
                "status": getattr(response, "status", 200),
                "bytes": int(length) if length and length.isdigit() else None,
                "content_type": response.headers.get_content_type(),
            }
    except urllib.error.HTTPError as exc:
        if exc.code in {403, 405, 501}:
            return {"status": None, "bytes": None, "content_type": ""}
        raise


def probe_curl(url: str, timeout_seconds: int) -> dict[str, Any]:
    if not shutil.which("curl"):
        raise RuntimeError("curl transport was requested but curl is unavailable")
    command = [
        "curl", "--fail", "--location", "--silent", "--show-error",
        "--head", "--connect-timeout", "30", "--max-time",
        str(timeout_seconds), "--user-agent", USER_AGENT, url,
    ]
    try:
        result = subprocess.run(
            command, check=True, capture_output=True, text=True
        )
    except subprocess.CalledProcessError:
        return {"status": None, "bytes": None, "content_type": ""}
    lengths = []
    content_types = []
    for line in result.stdout.splitlines():
        key, separator, value = line.partition(":")
        if not separator:
            continue
        if key.strip().lower() == "content-length" and value.strip().isdigit():
            lengths.append(int(value.strip()))
        elif key.strip().lower() == "content-type":
            content_types.append(value.strip().split(";", 1)[0])
    return {
        "status": 200,
        "bytes": lengths[-1] if lengths else None,
        "content_type": content_types[-1] if content_types else "",
    }


def check_disk_space(
    raw_dir: Path,
    remote_bytes: int | None,
    minimum_free_bytes: int,
    unknown_size_reserve_bytes: int,
) -> None:
    free = shutil.disk_usage(raw_dir).free
    reserve = remote_bytes if remote_bytes is not None else unknown_size_reserve_bytes
    if free - reserve < minimum_free_bytes:
        raise OSError(
            "disk guard refused download: "
            f"free={free}, anticipated={reserve}, minimum_after={minimum_free_bytes}"
        )


def stream_download(
    url: str,
    partial: Path,
    timeout_seconds: int,
    minimum_free_bytes: int,
) -> tuple[int, str]:
    bytes_written = 0
    next_space_check = FREE_SPACE_CHECK_BYTES
    with urllib.request.urlopen(
        request(url), timeout=timeout_seconds, context=TLS_CONTEXT
    ) as response:
        content_type = response.headers.get_content_type()
        if content_type == "text/html" and partial.suffix.lower() not in {
            ".html",
            ".htm",
        }:
            raise ValueError("server returned HTML instead of the requested data file")
        with partial.open("xb") as output:
            while True:
                chunk = response.read(CHUNK_BYTES)
                if not chunk:
                    break
                output.write(chunk)
                bytes_written += len(chunk)
                if bytes_written >= next_space_check:
                    if shutil.disk_usage(partial.parent).free < minimum_free_bytes:
                        raise OSError("disk free space fell below the configured safety floor")
                    next_space_check += FREE_SPACE_CHECK_BYTES
            output.flush()
            os.fsync(output.fileno())
    return bytes_written, content_type


def stream_download_curl(
    url: str,
    partial: Path,
    timeout_seconds: int,
    minimum_free_bytes: int,
) -> tuple[int, str]:
    del minimum_free_bytes  # preflight disk guard remains authoritative
    if not shutil.which("curl"):
        raise RuntimeError("curl transport was requested but curl is unavailable")
    command = [
        "curl", "--fail", "--location", "--silent", "--show-error",
        "--connect-timeout", "30", "--max-time", str(timeout_seconds),
        "--retry", "2", "--retry-delay", "1", "--output", str(partial),
        "--user-agent", USER_AGENT, url,
    ]
    subprocess.run(command, check=True, capture_output=True, text=True)
    return partial.stat().st_size, ""


def base_manifest_row(config: dict[str, Any], source: dict[str, Any]) -> dict[str, str]:
    metadata: dict[str, Any] = {}
    for rule in config.get("source_role_rules", []):
        exact_match = "id" in rule and rule["id"] == source.get("id")
        prefix_match = "id_prefix" in rule and str(source.get("id", "")).startswith(
            str(rule["id_prefix"])
        )
        if exact_match or prefix_match:
            metadata.update(rule)
    metadata.update(source)
    role_defaults = config.get("source_role_defaults", {})
    if isinstance(role_defaults, dict):
        defaults = role_defaults.get(str(metadata.get("role", "")), {})
        if isinstance(defaults, dict):
            for key, value in defaults.items():
                metadata.setdefault(key, value)
    required_analysis_ids = set(
        config.get("analysis_requirements", {}).get(
            "required_non_download_source_ids", []
        )
    )
    return {
        "application_id": str(config["application_id"]),
        "source_id": str(source["id"]),
        "source_organization": str(source.get("source_organization", "")),
        "source_page": str(source.get("source_page", "")),
        "requested_url": str(source.get("url", "")),
        "filename": str(source.get("filename", "")),
        "acquisition_mode": str(source.get("mode", "download")),
        "required_for_proposal": bool_text(source.get("required_for_proposal", False)),
        "enabled": bool_text(source.get("enabled", True)),
        "snapshot_date": str(config["snapshot_date"]),
        "status": "",
        "downloaded_at_utc": "",
        "bytes": "",
        "sha256": "",
        "content_type": "",
        "relative_path": "",
        "format": str(source.get("format", "binary")),
        "source_role": str(metadata.get("role", "")),
        "temporal_role": str(metadata.get("temporal_role", "")),
        "spatial_resolution": str(metadata.get("spatial_resolution", "")),
        "allowed_use": str(metadata.get("allowed_use", "")),
        "required_for_analysis": bool_text(source.get("id") in required_analysis_ids),
        "forbid_downscaling": (
            bool_text(metadata["forbid_downscaling"])
            if "forbid_downscaling" in metadata
            else ""
        ),
        "duplicate_group": str(metadata.get("duplicate_group", "")),
        "license": str(source.get("license", "")),
        "notes": str(source.get("notes", "")),
        "error": "",
    }


def acquire_source(
    config: dict[str, Any], source: dict[str, Any], app_root: Path, raw_dir: Path
) -> dict[str, str]:
    row = base_manifest_row(config, source)
    mode = source.get("mode", "download")
    enabled = source.get("enabled", True)
    if not enabled or mode in {"catalog_only", "manual_request", "restricted"}:
        row["status"] = "catalog_only" if not enabled else mode
        return row

    url, secrets = credentialized_url(source)
    if url is None:
        row["status"] = "requires_credentials"
        row["error"] = f"environment variable {source['api_key_env']} is not set"
        return row

    filename = source.get("filename", "")
    if not filename or Path(filename).name != filename:
        row["status"] = "configuration_error"
        row["error"] = "filename must be a non-empty basename"
        return row

    destination = raw_dir / filename
    row["relative_path"] = destination.relative_to(app_root).as_posix()
    if destination.exists():
        try:
            validate_download(
                destination,
                source.get("format", "binary"),
                int(source.get("min_bytes", 1)),
            )
            row.update(
                status="preserved_existing",
                downloaded_at_utc=dt.datetime.fromtimestamp(
                    destination.stat().st_mtime, dt.timezone.utc
                )
                .replace(microsecond=0)
                .isoformat()
                .replace("+00:00", "Z"),
                bytes=str(destination.stat().st_size),
                sha256=sha256_file(destination),
            )
        except Exception as exc:  # existing raw data must never be replaced
            row["status"] = "existing_file_invalid"
            row["error"] = sanitize(str(exc), secrets)
        return row

    if mode == "local_file":
        row["status"] = "missing_local_file"
        row["error"] = (
            "place the independently acquired file at the declared immutable "
            "raw path, then rerun acquisition"
        )
        return row

    timeout_seconds = int(config.get("timeout_seconds", 90))
    minimum_free_bytes = int(
        float(config.get("minimum_free_gib_after_download", 10)) * 1024**3
    )
    unknown_reserve_bytes = int(
        float(config.get("unknown_size_reserve_gib", 1)) * 1024**3
    )
    retries = int(config.get("retries", 3))
    partial = destination.with_name(destination.name + ".part")
    if partial.exists():
        row["status"] = "partial_file_exists"
        row["error"] = "remove or inspect the stale .part file manually"
        return row

    last_error: Exception | None = None
    pre_download_delay_seconds = float(source.get("pre_download_delay_seconds", 0))
    if source.get("skip_probe", False):
        pre_download_delay_seconds = max(
            pre_download_delay_seconds,
            float(config.get("minimum_dynamic_download_delay_seconds", 0)),
        )
    if pre_download_delay_seconds < 0:
        row["status"] = "configuration_error"
        row["error"] = "pre_download_delay_seconds cannot be negative"
        return row
    if pre_download_delay_seconds:
        time.sleep(pre_download_delay_seconds)
    for attempt in range(1, retries + 1):
        try:
            transport = source.get("transport", "urllib")
            if source.get("skip_probe", False):
                remote = {"status": None, "bytes": None, "content_type": ""}
            else:
                remote = (
                    probe_curl(url, timeout_seconds)
                    if transport == "curl"
                    else probe(url, timeout_seconds)
                )
            check_disk_space(
                raw_dir,
                remote.get("bytes"),
                minimum_free_bytes,
                unknown_reserve_bytes,
            )
            if transport == "curl":
                bytes_written, content_type = stream_download_curl(
                    url, partial, timeout_seconds, minimum_free_bytes
                )
            else:
                bytes_written, content_type = stream_download(
                    url, partial, timeout_seconds, minimum_free_bytes
                )
            validate_download(
                partial,
                source.get("format", "binary"),
                int(source.get("min_bytes", 1)),
            )
            partial.replace(destination)
            row.update(
                status="downloaded",
                downloaded_at_utc=utc_now(),
                bytes=str(bytes_written),
                sha256=sha256_file(destination),
                content_type=content_type,
            )
            return row
        except urllib.error.HTTPError as exc:
            last_error = exc
            if exc.code in {400, 401, 403, 404, 410}:
                break
        except Exception as exc:
            last_error = exc
        finally:
            if partial.exists():
                partial.unlink()
        if attempt < retries:
            retry_delay_seconds = float(
                source.get("retry_delay_seconds", 2 ** (attempt - 1))
            )
            if source.get("skip_probe", False):
                retry_delay_seconds = max(
                    retry_delay_seconds,
                    float(config.get("minimum_dynamic_download_delay_seconds", 0)),
                )
            time.sleep(retry_delay_seconds)

    row["status"] = "download_failed"
    row["error"] = sanitize(str(last_error or "unknown download error"), secrets)
    return row


def write_csv_once(path: Path, rows: list[dict[str, str]]) -> None:
    if path.exists():
        raise FileExistsError(f"refusing to overwrite {path}")
    with path.open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run_config(config_path: Path) -> tuple[list[dict[str, str]], Path]:
    with config_path.open("r", encoding="utf-8") as stream:
        config = json.load(stream)
    app_root = config_path.parent.parent
    raw_dir = app_root / "data" / "raw" / str(config["snapshot_date"])
    log_dir = app_root / "logs" / "v1"
    raw_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    version = next_manifest_version(raw_dir)
    rows = [acquire_source(config, source, app_root, raw_dir) for source in config["sources"]]
    manifest = raw_dir / f"source_manifest_v{version}.csv"
    write_csv_once(manifest, rows)

    event_log = log_dir / f"download_log_{utc_now().replace(':', '')}_v{version}.jsonl"
    if event_log.exists():
        raise FileExistsError(f"refusing to overwrite {event_log}")
    with event_log.open("x", encoding="utf-8") as stream:
        for row in rows:
            stream.write(json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n")
    return rows, manifest


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        action="append",
        required=True,
        type=Path,
        help="Path to an application config/sources.json (repeatable).",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Return non-zero if an enabled required download is incomplete.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    incomplete_required = False
    for config_path in args.config:
        rows, manifest = run_config(config_path.resolve())
        print(f"Manifest: {manifest}")
        for row in rows:
            print(f"  {row['source_id']}: {row['status']}")
            if (
                row["required_for_proposal"] == "true"
                and row["enabled"] == "true"
                and row["acquisition_mode"] == "download"
                and row["status"] not in {"downloaded", "preserved_existing"}
            ):
                incomplete_required = True
    return 1 if args.strict and incomplete_required else 0


if __name__ == "__main__":
    raise SystemExit(main())
