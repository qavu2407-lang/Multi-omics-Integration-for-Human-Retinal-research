from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd


@dataclass(slots=True)
class OmicsLayer:
    name: str
    matrix: pd.DataFrame
    sample_metadata: pd.DataFrame
    source_path: Path


@dataclass(slots=True)
class MultiOmicsDataset:
    views: dict[str, pd.DataFrame]
    sample_metadata: pd.DataFrame
    config_path: Path | None = None

