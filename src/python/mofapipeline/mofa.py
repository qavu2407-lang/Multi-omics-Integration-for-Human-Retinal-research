from __future__ import annotations

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from .config import PipelineConfig, resolve_path
from .models import MultiOmicsDataset
from .utils import ensure_dir


LOGGER = logging.getLogger(__name__)


def _import_mofa() -> tuple[object, object]:
    try:
        from mofapy2.run.entry_point import entry_point
        import mofax as mfx
    except Exception as exc:  # pragma: no cover
        raise ImportError(
            "MOFA+ dependencies are unavailable. Install with `pip install .[mofa]`."
        ) from exc
    return entry_point, mfx


def prepare_nested_data(dataset: MultiOmicsDataset, view_order: list[str]) -> tuple[list[list[object]], list[list[str]], list[str]]:
    samples = dataset.sample_metadata["sample_id"].tolist()
    features = [dataset.views[view_name].index.tolist() for view_name in view_order]
    data_nested = [[dataset.views[view_name].T.to_numpy()] for view_name in view_order]
    return data_nested, [samples], features


def train_mofa_run(
    dataset: MultiOmicsDataset,
    config: PipelineConfig,
    run_id: str,
    seed: int | None = None,
    views_override: dict[str, pd.DataFrame] | None = None,
) -> Path:
    entry_point, _ = _import_mofa()
    active_views = views_override or dataset.views
    view_order = list(config.mofa["views"])
    temp_dataset = MultiOmicsDataset(
        views={view_name: active_views[view_name] for view_name in view_order},
        sample_metadata=dataset.sample_metadata.copy(),
        config_path=dataset.config_path,
    )
    data_nested, samples_names, features_names = prepare_nested_data(temp_dataset, view_order)

    out_dir = ensure_dir(resolve_path(config, config.paths["mofa_dir"]))
    output_path = out_dir / f"{run_id}.hdf5"
    training_seed = int(config.analysis["random_seed"] if seed is None else seed)

    ent = entry_point()
    ent.set_data_options(
        scale_views=bool(config.mofa["scale_views"]),
        scale_groups=bool(config.mofa["scale_groups"]),
    )
    ent.set_data_matrix(
        data=data_nested,
        views_names=view_order,
        groups_names=["group1"],
        samples_names=samples_names,
        features_names=features_names,
    )
    ent.set_model_options(
        factors=int(config.mofa["factors"]),
        spikeslab_weights=bool(config.mofa["spikeslab_weights"]),
        ard_factors=bool(config.mofa["ard_factors"]),
        ard_weights=bool(config.mofa["ard_weights"]),
    )
    ent.set_train_options(
        iter=int(config.mofa["maxiter"]),
        convergence_mode=str(config.mofa["convergence_mode"]),
        dropR2=float(config.mofa["drop_r2"]),
        seed=training_seed,
        verbose=False,
    )
    ent.build()
    ent.run()
    ent.save(str(output_path))
    LOGGER.info("Saved MOFA run '%s' to %s.", run_id, output_path)
    return output_path


def load_mofa_model(model_path: str | Path):
    _, mfx = _import_mofa()
    return mfx.mofa_model(str(model_path))


def summarise_variance_explained(model, model_id: str) -> pd.DataFrame:
    variance = model.get_r2()
    summary = variance.groupby("View", as_index=False)["R2"].sum()
    summary.columns = ["view", "r2"]
    summary["model"] = model_id
    return summary


def extract_ranked_weights(model, factor: str = "Factor1", top_n: int = 15) -> pd.DataFrame:
    weights = model.get_weights(views=None, factors=factor, df=True)
    value_columns = [column for column in weights.columns if column not in {"view", "feature", "factor"}]
    if not value_columns:
        raise ValueError("Unable to identify the MOFA weight column.")
    ranked = (
        weights.rename(columns={value_columns[0]: "weight"})
        .assign(abs_weight=lambda df: df["weight"].abs())
        .sort_values("abs_weight", ascending=False)
        .head(top_n)
        .reset_index(drop=True)
    )
    ranked["rank"] = range(1, len(ranked) + 1)
    return ranked


def _save(fig: plt.Figure, path: Path) -> Path:
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_variance_explained_comparison(model_full, model_filtered, path: Path) -> Path:
    bind_df = pd.concat(
        [
            summarise_variance_explained(model_full, "Full data"),
            summarise_variance_explained(model_filtered, "Filtered proteomics"),
        ],
        ignore_index=True,
    )
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.barplot(data=bind_df, x="view", y="r2", hue="model", ax=ax)
    ax.set_title("Variance explained by omics view")
    ax.set_xlabel("Omics layer")
    ax.set_ylabel("Total variance explained (R2)")
    fig.tight_layout()
    return _save(fig, path)


def plot_ranked_weights(model, run_label: str, path: Path, factor: str = "Factor1") -> Path:
    ranked = extract_ranked_weights(model, factor=factor, top_n=15)
    fig, ax = plt.subplots(figsize=(9, 5))
    sns.lineplot(data=ranked, x="rank", y="weight", hue="view", marker="o", ax=ax)
    for _, row in ranked.iterrows():
        ax.text(row["rank"], row["weight"], row["feature"], fontsize=8)
    ax.set_title(f"{run_label} {factor} weights")
    fig.tight_layout()
    return _save(fig, path)


def plot_factor_scatter(
    model,
    sample_metadata: pd.DataFrame,
    run_label: str,
    path: Path,
    factor_x: str = "Factor1",
    factor_y: str = "Factor2",
) -> Path:
    factors_df = model.get_factors(factors=[factor_x, factor_y], df=True).rename(columns={"sample": "sample_id"})
    factors_df = factors_df.rename(columns={factor_x: "factor_x", factor_y: "factor_y"})
    plot_df = factors_df.merge(sample_metadata, on="sample_id", how="left")

    fig, ax = plt.subplots(figsize=(8, 6))
    sns.scatterplot(data=plot_df, x="factor_x", y="factor_y", hue="group", s=80, ax=ax)
    for _, row in plot_df.iterrows():
        ax.text(row["factor_x"], row["factor_y"], row["sample_id"], fontsize=8, alpha=0.7)
    ax.set_title(f"{run_label} latent factor separation")
    ax.set_xlabel(factor_x)
    ax.set_ylabel(factor_y)
    fig.tight_layout()
    return _save(fig, path)


def run_downstream_workflow(
    dataset: MultiOmicsDataset,
    config: PipelineConfig,
    run_full_path: str | Path,
    run_filtered_path: str | Path,
) -> list[Path]:
    plot_dir = ensure_dir(resolve_path(config, config.paths["plot_downstream_dir"]))
    mofa_dir = ensure_dir(resolve_path(config, config.paths["mofa_dir"]))

    full_model = load_mofa_model(run_full_path)
    filtered_model = load_mofa_model(run_filtered_path)

    comparison_df = pd.concat(
        [
            summarise_variance_explained(full_model, "run_01_full"),
            summarise_variance_explained(filtered_model, "run_02_filtered_proteins"),
        ],
        ignore_index=True,
    )
    comparison_path = mofa_dir / "variance_explained_by_view.csv"
    comparison_df.to_csv(comparison_path, index=False)

    return [
        comparison_path,
        plot_variance_explained_comparison(full_model, filtered_model, plot_dir / "variance_explained_comparison.png"),
        plot_ranked_weights(full_model, "MOFA run 01", plot_dir / "run_01_factor1_weights.png"),
        plot_ranked_weights(filtered_model, "MOFA run 02", plot_dir / "run_02_factor1_weights.png"),
        plot_factor_scatter(full_model, dataset.sample_metadata, "MOFA run 01", plot_dir / "run_01_factors_F1_F2.png"),
        plot_factor_scatter(filtered_model, dataset.sample_metadata, "MOFA run 02", plot_dir / "run_02_factors_F1_F2.png"),
    ]

