from __future__ import annotations

import json
import threading
from pathlib import Path
from typing import Protocol

from .schemas import MemoryEpisode


class EpisodeStore(Protocol):
    def add(self, episode: MemoryEpisode) -> None: ...

    def list_all(self) -> list[MemoryEpisode]: ...


_LOCKS: dict[Path, threading.Lock] = {}
_LOCKS_GUARD = threading.Lock()


class JsonEpisodeStore:
    def __init__(self, path: Path) -> None:
        self.path = path
        resolved = path.resolve()
        with _LOCKS_GUARD:
            self._lock = _LOCKS.setdefault(resolved, threading.Lock())

    def add(self, episode: MemoryEpisode) -> None:
        with self._lock:
            episodes = self._read_unlocked()
            episodes.append(episode)
            self.path.parent.mkdir(parents=True, exist_ok=True)
            temporary = self.path.with_suffix(self.path.suffix + ".tmp")
            temporary.write_text(
                json.dumps([item.to_dict() for item in episodes], ensure_ascii=True, indent=2),
                encoding="utf-8",
            )
            temporary.replace(self.path)

    def list_all(self) -> list[MemoryEpisode]:
        with self._lock:
            return self._read_unlocked()

    def _read_unlocked(self) -> list[MemoryEpisode]:
        if not self.path.is_file():
            return []
        try:
            payload = json.loads(self.path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return []
        if not isinstance(payload, list):
            return []
        episodes: list[MemoryEpisode] = []
        required = set(MemoryEpisode.field_names())
        for item in payload:
            if isinstance(item, dict) and set(item) == required:
                episodes.append(MemoryEpisode.from_dict(item))
        return episodes
