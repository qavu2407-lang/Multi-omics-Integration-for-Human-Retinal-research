from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


@dataclass(slots=True)
class PipelineConfig:
    raw: dict[str, Any]
    path: Path

    @property
    def paths(self) -> dict[str, Any]:
        return self.raw["paths"]

    @property
    def data(self) -> dict[str, Any]:
        return self.raw["data"]

    @property
    def analysis(self) -> dict[str, Any]:
        return self.raw["analysis"]

    @property
    def mofa(self) -> dict[str, Any]:
        return self.raw["mofa"]

    @property
    def evaluation(self) -> dict[str, Any]:
        return self.raw["evaluation"]


def read_config(path: str | Path) -> PipelineConfig:
    config_path = Path(path).resolve()
    with config_path.open("r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle)

    required_sections = {"paths", "data", "analysis", "mofa", "evaluation"}
    missing = required_sections.difference(raw or {})
    if missing:
        raise ValueError(f"Missing required config sections: {sorted(missing)}")

    return PipelineConfig(raw=raw, path=config_path)


def resolve_path(config: PipelineConfig, value: str | Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return (config.path.parent / path).resolve()
