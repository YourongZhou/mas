from __future__ import annotations

import sys
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any


def now_stamp() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


@dataclass
class TimedBlock:
    logger: "RunLogger"
    label: str
    start: float = 0.0

    def __enter__(self) -> "TimedBlock":
        self.start = time.perf_counter()
        self.logger.log(f"START {self.label}")
        return self

    def __exit__(self, exc_type: Any, exc: Any, tb: Any) -> None:
        elapsed = time.perf_counter() - self.start
        status = "ERROR" if exc else "END"
        self.logger.log(f"{status} {self.label} elapsed={elapsed:.2f}s")


class RunLogger:
    def __init__(self, logs_dir: Path, *, run_id: str | None = None) -> None:
        self.run_id = run_id or f"run_{now_stamp()}"
        self.path = logs_dir / f"{self.run_id}.log"
        self._handle = None

    def __enter__(self) -> "RunLogger":
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._handle = self.path.open("w", encoding="utf-8-sig")
        self.log(f"RUN START id={self.run_id}")
        return self

    def __exit__(self, exc_type: Any, exc: Any, tb: Any) -> None:
        self.log(f"RUN END id={self.run_id}")
        if self._handle:
            self._handle.close()
            self._handle = None

    def log(self, message: str) -> None:
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        line = f"[{ts}] {message}"
        encoding = sys.stdout.encoding or "utf-8"
        console_line = line.encode(encoding, errors="replace").decode(encoding, errors="replace")
        sys.stdout.write(console_line + "\n")
        sys.stdout.flush()
        if self._handle:
            self._handle.write(line + "\n")
            self._handle.flush()

    def section(self, title: str) -> None:
        self.log("=" * 88)
        self.log(title)
        self.log("=" * 88)

    def node(self, title: str, detail: str = "") -> None:
        self.section(f"节点：{title}")
        if detail:
            self.log(f"说明：{detail}")

    def progress(self, title: str, detail: str = "") -> None:
        suffix = f" | {detail}" if detail else ""
        self.log(f"进度：{title}{suffix}")

    def warning(self, message: str) -> None:
        self.log(f"警告：{message}")

    def error_reason(self, message: str) -> None:
        reason = message.strip() if message else "未能自动提取明确错误原因，请查看 stdout/stderr 详情。"
        self.log(f"错误原因：{reason}")

    def timed(self, label: str) -> TimedBlock:
        return TimedBlock(self, label)

    def preview(self, title: str, text: str, max_chars: int = 4000) -> None:
        body = str(text or "")
        if len(body) > max_chars:
            body = body[:max_chars] + "\n...[truncated]"
        self.log(f"{title}:\n{body}")


class TeeStdout:
    """Optional stdout tee for notebook cells that also need a complete log file."""

    def __init__(self, logger: RunLogger) -> None:
        self.logger = logger
        self.original = sys.stdout

    def write(self, text: str) -> int:
        written = self.original.write(text)
        if self.logger._handle:
            self.logger._handle.write(text)
            self.logger._handle.flush()
        return written

    def flush(self) -> None:
        self.original.flush()

