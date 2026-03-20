from __future__ import annotations

import json
import logging
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


def ensure_dir(path: str | Path) -> Path:
    out = Path(path)
    out.mkdir(parents=True, exist_ok=True)
    return out


def configure_logging(log_path: str | Path | None = None) -> None:
    handlers: list[logging.Handler] = [logging.StreamHandler(sys.stdout)]
    if log_path is not None:
        handlers.append(logging.FileHandler(log_path, encoding="utf-8"))

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=handlers,
        force=True,
    )


def save_json(data: dict[str, Any], path: str | Path) -> Path:
    output_path = Path(path)
    output_path.write_text(json.dumps(data, indent=2), encoding="utf-8")
    return output_path


def write_session_metadata(path: str | Path) -> Path:
    payload = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": sys.version,
        "platform": platform.platform(),
        "numpy_version": np.__version__,
    }
    return save_json(payload, path)


def safe_corrcoef(x: np.ndarray, y: np.ndarray) -> float:
    if x.ndim != 1 or y.ndim != 1:
        raise ValueError("safe_corrcoef expects 1D arrays.")
    if np.std(x) == 0 or np.std(y) == 0:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])

