from __future__ import annotations

import re

from .long_term import EpisodeStore
from .schemas import MemoryEpisode


class EpisodeRetriever:
    def __init__(self, store: EpisodeStore) -> None:
        self.store = store

    def retrieve(self, task_signature: str, data_signature: str, *, limit: int = 3) -> list[MemoryEpisode]:
        query_tokens = _tokens(f"{task_signature} {data_signature}")
        ranked: list[tuple[float, str, MemoryEpisode]] = []
        for episode in self.store.list_all():
            episode_tokens = _tokens(
                " ".join(
                    [
                        episode.task_signature,
                        episode.data_signature,
                        episode.skill_id,
                        episode.runtime_environment,
                        episode.reusable_lesson,
                    ]
                )
            )
            overlap = len(query_tokens & episode_tokens)
            score = float(overlap)
            if data_signature and data_signature == episode.data_signature:
                score += 4.0
            if score > 0:
                ranked.append((score, episode.timestamp, episode))
        ranked.sort(key=lambda item: (item[0], item[1]), reverse=True)
        return [item[2] for item in ranked[: max(0, limit)]]

    def format_for_prompt(self, episodes: list[MemoryEpisode], *, max_chars: int = 1800) -> str:
        lines = [
            "Prior verified experience:",
            "Current observations always override historical experience.",
        ]
        if not episodes:
            lines.append("- No relevant verified episodes were found.")
        for episode in episodes:
            parts = [
                f"task={episode.task_signature}",
                f"data={episode.data_signature or '(unknown)'}",
                f"skill={episode.skill_id or '(none)'}",
                f"runtime={episode.runtime_environment or '(unknown)'}",
                f"outcome={episode.outcome}",
            ]
            if episode.verified_root_cause:
                parts.append(f"root_cause={episode.verified_root_cause}")
            if episode.verified_fix:
                parts.append(f"fix={episode.verified_fix}")
            if episode.reusable_lesson:
                parts.append(f"lesson={episode.reusable_lesson}")
            lines.append("- " + "; ".join(parts))
        return "\n".join(lines)[:max_chars]


def _tokens(value: str) -> set[str]:
    return set(re.findall(r"[a-z0-9_.-]+|[\u4e00-\u9fff]", value.lower()))
