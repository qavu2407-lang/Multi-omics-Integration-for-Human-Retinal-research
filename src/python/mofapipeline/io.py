from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from .config import PipelineConfig, resolve_path
from .models import MultiOmicsDataset, OmicsLayer


LOGGER = logging.getLogger(__name__)


def read_omics_layer(
    path: str | Path,
    view_name: str,
    id_column: str = "sample",
    label_token: str = "label",
) -> OmicsLayer:
    csv_path = Path(path)
    raw_df = pd.read_csv(csv_path)

    if id_column not in raw_df.columns:
        raise ValueError(f"Column '{id_column}' not found in {csv_path}.")

    label_rows = raw_df.loc[raw_df[id_column] == label_token]
    if label_rows.shape[0] != 1:
        raise ValueError(
            f"Expected exactly one label row in {csv_path}, found {label_rows.shape[0]}."
        )

    feature_df = raw_df.loc[raw_df[id_column] != label_token].copy()
    sample_ids = [col for col in feature_df.columns if col != id_column]

    if feature_df[id_column].duplicated().any():
        raise ValueError(f"Duplicated feature identifiers detected in {csv_path}.")
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError(f"Duplicated sample identifiers detected in {csv_path}.")

    for sample_id in sample_ids:
        feature_df[sample_id] = pd.to_numeric(feature_df[sample_id], errors="coerce")

    matrix = feature_df.set_index(id_column)
    metadata = pd.DataFrame(
        {
            "sample_id": sample_ids,
            "group": label_rows.iloc[0][sample_ids].astype(str).to_list(),
            "view": view_name,
        }
    )

    LOGGER.info(
        "Loaded %s from %s with %s features and %s samples.",
        view_name,
        csv_path,
        matrix.shape[0],
        matrix.shape[1],
    )

    return OmicsLayer(
        name=view_name,
        matrix=matrix,
        sample_metadata=metadata,
        source_path=csv_path.resolve(),
    )


def align_multiomics_layers(layers: list[OmicsLayer]) -> MultiOmicsDataset:
    sample_sets = [set(layer.matrix.columns) for layer in layers]
    common_samples = sorted(set.intersection(*sample_sets))
    if not common_samples:
        raise ValueError("No shared samples were found across omics layers.")

    aligned_views: dict[str, pd.DataFrame] = {}
    aligned_meta: list[pd.DataFrame] = []
    for layer in layers:
        aligned_views[layer.name] = layer.matrix.loc[:, common_samples].copy()
        meta = (
            layer.sample_metadata.loc[layer.sample_metadata["sample_id"].isin(common_samples)]
            .set_index("sample_id")
            .loc[common_samples]
            .reset_index()
        )
        aligned_meta.append(meta)

    reference = aligned_meta[0][["sample_id", "group"]].reset_index(drop=True)
    for current in aligned_meta[1:]:
        current_cmp = current[["sample_id", "group"]].reset_index(drop=True)
        if not reference.equals(current_cmp):
            raise ValueError("Sample group labels are not consistent across omics layers.")

    return MultiOmicsDataset(views=aligned_views, sample_metadata=reference)


def load_multiomics_data(config: PipelineConfig) -> MultiOmicsDataset:
    layer_specs = {
        "proteomics": config.paths["proteomics"],
        "lipidomics": config.paths["lipidomics"],
        "metabolomics": config.paths["metabolomics"],
    }
    layers = [
        read_omics_layer(
            resolve_path(config, relative_path),
            view_name=view_name,
            id_column=config.data["id_column"],
            label_token=config.data["label_token"],
        )
        for view_name, relative_path in layer_specs.items()
    ]
    dataset = align_multiomics_layers(layers)
    dataset.config_path = config.path
    return dataset

