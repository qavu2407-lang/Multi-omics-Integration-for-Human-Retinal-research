from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from .config import PipelineConfig, resolve_path
from .models import MultiOmicsDataset
from .utils import ensure_dir


def _save(fig: plt.Figure, path: Path) -> Path:
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_distributions(matrix: pd.DataFrame, omic_name: str, save_path: Path) -> Path:
    numeric = matrix.apply(pd.to_numeric, errors="coerce")
    fig, axes = plt.subplots(1, 4, figsize=(16, 3.5))
    flat_values = numeric.to_numpy().ravel()

    axes[0].hist(flat_values[~pd.isna(flat_values)], bins=50, edgecolor="black")
    axes[0].set_title(f"{omic_name}: Overall distribution")

    axes[1].boxplot([numeric[col].dropna().to_numpy() for col in numeric.columns], showfliers=False)
    axes[1].set_title(f"{omic_name}: Sample distributions")
    axes[1].tick_params(axis="x", rotation=90)

    for column in numeric.columns[: min(5, numeric.shape[1])]:
        numeric[column].dropna().plot(kind="density", ax=axes[2], alpha=0.5)
    axes[2].set_title(f"{omic_name}: Sample densities")

    feature_variance = numeric.var(axis=1, skipna=True).sort_values(ascending=False)
    axes[3].plot(feature_variance.to_numpy())
    axes[3].set_title(f"{omic_name}: Feature variances")
    axes[3].set_yscale("log")

    fig.tight_layout()
    return _save(fig, save_path)


def plot_top_variable_features(
    matrix: pd.DataFrame,
    omic_name: str,
    save_path: Path,
    n_top: int,
) -> Path:
    numeric = matrix.apply(pd.to_numeric, errors="coerce")
    feature_variance = numeric.var(axis=1, skipna=True).sort_values(ascending=False)
    top_features = feature_variance.head(min(n_top, len(feature_variance)))

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    axes[0].barh(range(len(top_features)), top_features.to_numpy())
    axes[0].set_yticks(range(len(top_features)))
    axes[0].set_yticklabels(top_features.index)
    axes[0].invert_yaxis()
    axes[0].set_title(f"{omic_name}: Top variable features")

    sns.heatmap(numeric.loc[top_features.index], center=0, ax=axes[1], cbar_kws={"label": "Expression"})
    axes[1].set_title(f"{omic_name}: Top feature heatmap")

    fig.tight_layout()
    return _save(fig, save_path)


def plot_pca(
    matrix: pd.DataFrame,
    labels: pd.Series,
    omic_name: str,
    save_path: Path,
    n_components: int = 2,
) -> Path:
    data = matrix.T.apply(pd.to_numeric, errors="coerce")
    scaled = StandardScaler().fit_transform(data)
    pcs = PCA(n_components=n_components).fit_transform(scaled)
    plot_df = pd.DataFrame(
        {
            "PC1": pcs[:, 0],
            "PC2": pcs[:, 1],
            "group": labels.to_numpy(),
            "sample_id": data.index.to_numpy(),
        }
    )

    fig, ax = plt.subplots(figsize=(8, 6))
    sns.scatterplot(data=plot_df, x="PC1", y="PC2", hue="group", s=80, ax=ax)
    for _, row in plot_df.iterrows():
        ax.text(row["PC1"], row["PC2"], row["sample_id"], fontsize=8, alpha=0.7)
    ax.set_title(f"{omic_name}: PCA")
    fig.tight_layout()
    return _save(fig, save_path)


def run_eda_workflow(dataset: MultiOmicsDataset, config: PipelineConfig) -> list[Path]:
    output_dir = ensure_dir(resolve_path(config, config.paths["plot_eda_dir"]))
    outputs: list[Path] = []
    labels = dataset.sample_metadata.set_index("sample_id")["group"]

    for view_name, matrix in dataset.views.items():
        outputs.append(plot_distributions(matrix, view_name.title(), output_dir / f"{view_name}_distribution.png"))
        outputs.append(
            plot_top_variable_features(
                matrix,
                view_name.title(),
                output_dir / f"{view_name}_top_features.png",
                n_top=int(config.analysis["top_variable_features"]),
            )
        )
        outputs.append(
            plot_pca(
                matrix,
                labels.loc[matrix.columns],
                view_name.title(),
                output_dir / f"pca_{view_name}.png",
                n_components=int(config.analysis["pca_components"]),
            )
        )
    return outputs

