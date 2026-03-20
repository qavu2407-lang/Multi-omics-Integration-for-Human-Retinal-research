from __future__ import annotations

import numpy as np
import pandas as pd

from .config import PipelineConfig, resolve_path
from .models import MultiOmicsDataset
from .utils import ensure_dir, write_session_metadata


def compute_missingness(matrix: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "feature_id": matrix.index,
            "missing_fraction": matrix.isna().mean(axis=1).to_numpy(),
            "variance": matrix.var(axis=1, skipna=True).to_numpy(),
        }
    )


def impute_matrix(matrix: pd.DataFrame, method: str = "median") -> pd.DataFrame:
    if method != "median":
        raise ValueError(f"Unsupported imputation method: {method}")

    medians = matrix.median(axis=1, skipna=True)
    imputed = matrix.T.fillna(medians).T
    return imputed


def filter_matrix_by_missingness(
    matrix: pd.DataFrame,
    feature_threshold: float,
    sample_threshold: float,
) -> pd.DataFrame:
    feature_keep = matrix.isna().mean(axis=1) <= feature_threshold
    sample_keep = matrix.isna().mean(axis=0) <= sample_threshold
    return matrix.loc[feature_keep, sample_keep]


def build_qc_summary(dataset: MultiOmicsDataset) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for view_name, matrix in dataset.views.items():
        records.append(
            {
                "view": view_name,
                "n_features": int(matrix.shape[0]),
                "n_samples": int(matrix.shape[1]),
                "missing_values": int(matrix.isna().sum().sum()),
                "duplicated_features": bool(matrix.index.duplicated().any()),
                "duplicated_samples": bool(matrix.columns.duplicated().any()),
                "median_feature_missingness": float(matrix.isna().mean(axis=1).median()),
                "median_sample_missingness": float(matrix.isna().mean(axis=0).median()),
            }
        )
    return pd.DataFrame.from_records(records)


def run_qc_workflow(dataset: MultiOmicsDataset, config: PipelineConfig) -> list[Path]:
    qc_dir = ensure_dir(resolve_path(config, config.paths["qc_dir"]))
    outputs: list[Path] = []

    qc_summary_path = qc_dir / "qc_summary.csv"
    build_qc_summary(dataset).to_csv(qc_summary_path, index=False)
    outputs.append(qc_summary_path)

    metadata_path = qc_dir / "sample_metadata.csv"
    dataset.sample_metadata.to_csv(metadata_path, index=False)
    outputs.append(metadata_path)

    for view_name, matrix in dataset.views.items():
        missingness_path = qc_dir / f"{view_name}_feature_missingness.csv"
        compute_missingness(matrix).to_csv(missingness_path, index=False)
        outputs.append(missingness_path)

    outputs.append(write_session_metadata(qc_dir / "session_info.json"))
    return outputs


def apply_variance_filter(
    views: dict[str, pd.DataFrame],
    config: PipelineConfig,
) -> dict[str, pd.DataFrame]:
    variance_filter = config.data.get("variance_filter", {})
    if not variance_filter.get("enabled", False):
        return views

    target_view = variance_filter["view"]
    if target_view not in views:
        raise ValueError(f"Variance filter target view '{target_view}' not found.")

    filtered = dict(views)
    target = filtered[target_view]
    ranked = target.var(axis=1, skipna=True).sort_values(ascending=False)
    keep_ids = ranked.head(min(int(variance_filter["top_n"]), len(ranked))).index
    filtered[target_view] = target.loc[keep_ids].copy()
    return filtered


def preprocess_multiomics(dataset: MultiOmicsDataset, config: PipelineConfig) -> MultiOmicsDataset:
    processed_views: dict[str, pd.DataFrame] = {}
    for view_name, matrix in dataset.views.items():
        filtered = filter_matrix_by_missingness(
            matrix,
            feature_threshold=float(config.data["feature_missingness_threshold"]),
            sample_threshold=float(config.data["sample_missingness_threshold"]),
        )
        if filtered.empty:
            raise ValueError(f"Filtering removed all data from view '{view_name}'.")
        imputed = impute_matrix(filtered, method=config.data["imputation_method"])
        if np.isnan(imputed.to_numpy()).all():
            raise ValueError(f"All values became NaN in view '{view_name}' after preprocessing.")
        processed_views[view_name] = imputed

    common_samples = sorted(set.intersection(*(set(df.columns) for df in processed_views.values())))
    if not common_samples:
        raise ValueError("No common samples remain after preprocessing.")

    processed_views = {
        view_name: df.loc[:, common_samples].copy() for view_name, df in processed_views.items()
    }
    processed_views = apply_variance_filter(processed_views, config)
    sample_metadata = (
        dataset.sample_metadata.loc[dataset.sample_metadata["sample_id"].isin(common_samples)]
        .set_index("sample_id")
        .loc[common_samples]
        .reset_index()
    )
    return MultiOmicsDataset(
        views=processed_views,
        sample_metadata=sample_metadata,
        config_path=dataset.config_path,
    )


def build_analysis_datasets(
    dataset: MultiOmicsDataset,
    config: PipelineConfig,
) -> tuple[MultiOmicsDataset, MultiOmicsDataset]:
    base_dataset = preprocess_multiomics(dataset, config)
    filtered_dataset = MultiOmicsDataset(
        views=apply_variance_filter(base_dataset.views, config),
        sample_metadata=base_dataset.sample_metadata.copy(),
        config_path=base_dataset.config_path,
    )
    return base_dataset, filtered_dataset
