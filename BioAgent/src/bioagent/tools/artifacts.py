from __future__ import annotations

import csv
import gzip
import hashlib
import json
import struct
import re
import subprocess
from pathlib import Path
from typing import Any

from bioagent.config import AgentConfig
from bioagent.skills import image_for_profile

from .workspace import resolve_run_path, run_relative_path


class ArtifactStore:
    def __init__(
        self,
        run_dir: Path,
        config: AgentConfig | None = None,
        *,
        prior_run_dirs: list[Path] | None = None,
    ) -> None:
        self.run_dir = run_dir
        self.config = config
        self.run_dirs = _artifact_run_dirs(run_dir, config, prior_run_dirs or [])
        self.path = run_dir / "state" / "artifact_evidence.json"

    def inspect(self, path: str, *, max_rows: int = 10) -> dict:
        target, source_run_dir = self._resolve_artifact_path(path)
        if not target.is_file():
            return {"ok": False, "error": f"Artifact not found: {target}"}
        digest = _sha256(target)
        facts = _inspect_file(
            target,
            max_rows=max(1, min(max_rows, 50)),
            run_dir=source_run_dir,
            config=self.config,
        )
        relative = run_relative_path(source_run_dir, target)
        record = {
            "path": str(target),
            "workspace_path": f"/work/{relative}" if source_run_dir == self.run_dir.resolve() else str(target),
            "source_run_id": source_run_dir.name,
            "sha256": digest,
            "facts": facts,
        }
        evidence_id = "evidence_" + hashlib.sha256(
            json.dumps(record, ensure_ascii=False, sort_keys=True, default=str).encode("utf-8")
        ).hexdigest()[:12]
        record["evidence_id"] = evidence_id
        records = self._load()
        records[evidence_id] = record
        self._save(records)
        return {"ok": True, **record}

    def _resolve_artifact_path(self, path: str) -> tuple[Path, Path]:
        raw = str(path or "").strip()
        if not raw or raw == "/work" or raw.startswith("/work/") or not Path(raw).is_absolute():
            return resolve_run_path(self.run_dir, raw), self.run_dir.resolve()
        target = Path(raw).resolve()
        for run_dir in self.run_dirs:
            if target == run_dir or run_dir in target.parents:
                return target, run_dir
        raise ValueError(f"Path is outside the session run workspaces: {path}")

    def finish_task(self, *, summary: str, evidence_ids: list[str]) -> dict:
        ids = list(dict.fromkeys(str(item).strip() for item in evidence_ids if str(item).strip()))
        if not ids:
            return {"ok": False, "error": "finish_task requires at least one inspect_artifact evidence_id."}
        records = self._load()
        selected: list[dict[str, Any]] = []
        for evidence_id in ids:
            record = records.get(evidence_id)
            if not isinstance(record, dict):
                return {"ok": False, "error": f"Unknown evidence_id: {evidence_id}"}
            target = Path(str(record.get("path") or ""))
            if not target.is_file():
                return {"ok": False, "error": f"Artifact no longer exists: {target}"}
            if _sha256(target) != record.get("sha256"):
                return {"ok": False, "error": f"Artifact changed since inspection: {target}. Inspect it again."}
            selected.append(record)
        unsupported_numbers = _unsupported_numbers(summary, selected)
        if unsupported_numbers:
            return {
                "ok": False,
                "error": (
                    "Summary contains numeric claims not present in the selected evidence: "
                    + ", ".join(unsupported_numbers)
                    + ". Inspect the supporting artifact or remove the unsupported claim."
                ),
            }
        lines = [summary.strip() or "Task completed with verified artifacts.", "", "Verified artifacts:"]
        for record in selected:
            facts = record.get("facts") or {}
            fact_text = ", ".join(f"{key}={value}" for key, value in facts.items() if key not in {"preview", "columns"})
            lines.append(f"- `{record['path']}` ({fact_text})")
        result = {
            "ok": True,
            "status": "verified",
            "summary": summary.strip(),
            "evidence_ids": ids,
            "evidence": selected,
            "final_answer": "\n".join(lines),
        }
        verification_path = self.run_dir / "state" / "final_verification.json"
        verification_path.parent.mkdir(parents=True, exist_ok=True)
        verification_path.write_text(json.dumps(result, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
        return result

    def _load(self) -> dict[str, dict[str, Any]]:
        try:
            value = json.loads(self.path.read_text(encoding="utf-8"))
        except (FileNotFoundError, json.JSONDecodeError, OSError):
            return {}
        return value if isinstance(value, dict) else {}

    def _save(self, records: dict[str, dict[str, Any]]) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        temporary = self.path.with_suffix(".tmp")
        temporary.write_text(json.dumps(records, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
        temporary.replace(self.path)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _inspect_file(
    path: Path,
    *,
    max_rows: int,
    run_dir: Path | None = None,
    config: AgentConfig | None = None,
) -> dict[str, Any]:
    suffixes = [item.lower() for item in path.suffixes]
    facts: dict[str, Any] = {"kind": "file", "size_bytes": path.stat().st_size}
    if path.suffix.lower() == ".json":
        value = json.loads(path.read_text(encoding="utf-8"))
        facts.update(kind="json", top_level_type=type(value).__name__)
        if isinstance(value, dict):
            facts["keys"] = list(value)[:50]
        elif isinstance(value, list):
            facts["item_count"] = len(value)
            facts["preview"] = value[:max_rows]
        return facts
    if path.suffix.lower() == ".h5ad":
        return {**facts, **_inspect_h5ad(path, run_dir=run_dir, config=config)}
    if path.suffix.lower() in {".png", ".jpg", ".jpeg"}:
        width, height = _image_dimensions(path)
        facts.update(kind="image", width=width, height=height)
        return facts
    if ".csv" in suffixes or ".tsv" in suffixes:
        delimiter = "\t" if ".tsv" in suffixes else ","
        opener = gzip.open if path.suffix.lower() == ".gz" else open
        with opener(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.reader(handle, delimiter=delimiter)
            header = next(reader, [])
            preview: list[list[str]] = []
            row_count = 0
            for row in reader:
                row_count += 1
                if len(preview) < max_rows:
                    preview.append(row)
        facts.update(
            kind="table",
            columns=header,
            row_count=row_count,
            column_count=len(header),
            preview=preview,
        )
        return facts
    if path.suffix.lower() == ".gz":
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            line_count = sum(1 for _ in handle)
        facts.update(kind="gzip_text", line_count=line_count)
        return facts
    if path.suffix.lower() in {".txt", ".md", ".log", ".py", ".r"}:
        text = path.read_text(encoding="utf-8", errors="replace")
        facts.update(kind="text", line_count=len(text.splitlines()), preview=text[:2000])
    return facts


def _image_dimensions(path: Path) -> tuple[int | None, int | None]:
    with path.open("rb") as handle:
        header = handle.read(24)
    if header.startswith(b"\x89PNG\r\n\x1a\n") and len(header) >= 24:
        return struct.unpack(">II", header[16:24])
    try:
        from PIL import Image

        with Image.open(path) as image:
            return image.size
    except Exception:
        return None, None


def _inspect_h5ad(
    path: Path,
    *,
    run_dir: Path | None = None,
    config: AgentConfig | None = None,
) -> dict[str, Any]:
    try:
        import h5py

        with h5py.File(path, "r") as handle:
            obs = handle.get("obs/_index")
            if obs is None:
                obs = handle.get("obs/index")
            var = handle.get("var/_index")
            if var is None:
                var = handle.get("var/index")
            return {
                "kind": "h5ad",
                "n_obs": len(obs) if obs is not None else None,
                "n_vars": len(var) if var is not None else None,
                "top_level_keys": list(handle.keys()),
            }
    except ImportError as exc:
        if run_dir is not None and config is not None:
            return _inspect_h5ad_in_docker(path, run_dir=run_dir, config=config)
        return {"kind": "h5ad", "inspection_error": str(exc)[:300]}
    except Exception as exc:
        return {"kind": "h5ad", "inspection_error": str(exc)[:300]}


def _inspect_h5ad_in_docker(path: Path, *, run_dir: Path, config: AgentConfig) -> dict[str, Any]:
    entry = image_for_profile(config.docker_root, "py-scverse-v1")
    image = str((entry or {}).get("image_tag") or "")
    if not image:
        return {"kind": "h5ad", "inspection_error": "py-scverse-v1 is not configured."}
    docker_path = f"/work/{run_relative_path(run_dir, path)}"
    code = (
        "import h5py,json,sys; "
        "f=h5py.File(sys.argv[1],'r'); "
        "obs=f.get('obs/_index') if f.get('obs/_index') is not None else f.get('obs/index'); "
        "var=f.get('var/_index') if f.get('var/_index') is not None else f.get('var/index'); "
        "print(json.dumps({'kind':'h5ad','n_obs':len(obs) if obs is not None else None,"
        "'n_vars':len(var) if var is not None else None,'top_level_keys':list(f.keys())})); f.close()"
    )
    try:
        result = subprocess.run(
            [
                "docker",
                "run",
                "--rm",
                "-v",
                f"{run_dir}:/work:ro",
                image,
                "python",
                "-c",
                code,
                docker_path,
            ],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=120,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {"kind": "h5ad", "inspection_error": str(exc)[:300]}
    if result.returncode != 0:
        return {"kind": "h5ad", "inspection_error": (result.stderr or result.stdout)[-500:]}
    try:
        return json.loads(result.stdout.strip().splitlines()[-1])
    except (IndexError, json.JSONDecodeError):
        return {"kind": "h5ad", "inspection_error": "Invalid py-scverse inspection output."}


def _unsupported_numbers(summary: str, evidence: list[dict[str, Any]]) -> list[str]:
    claims = set(re.findall(r"(?<![A-Za-z])\d+(?:\.\d+)?(?:[eE][+-]?\d+)?%?", summary or ""))
    if not claims:
        return []
    evidence_text = json.dumps(evidence, ensure_ascii=False, default=str)
    return sorted(value for value in claims if value not in evidence_text)


def _artifact_run_dirs(
    current_run: Path,
    config: AgentConfig | None,
    prior_run_dirs: list[Path],
) -> list[Path]:
    runs_root = config.runs_dir.resolve() if config is not None else current_run.resolve().parent
    result: list[Path] = []
    for candidate in [current_run, *prior_run_dirs]:
        resolved = Path(candidate).resolve()
        if resolved != runs_root and runs_root not in resolved.parents:
            continue
        if resolved not in result:
            result.append(resolved)
    return result
