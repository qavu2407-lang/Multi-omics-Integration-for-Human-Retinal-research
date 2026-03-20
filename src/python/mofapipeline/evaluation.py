from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .config import PipelineConfig, resolve_path
from .models import MultiOmicsDataset
from .mofa import load_mofa_model, summarise_variance_explained, train_mofa_run
from .utils import ensure_dir, safe_corrcoef


def subset_views_by_samples(
    views: dict[str, pd.DataFrame],
    sample_ids: list[str],
) -> dict[str, pd.DataFrame]:
    return {view_name: matrix.loc[:, sample_ids].copy() for view_name, matrix in views.items()}


def seed_stability_analysis(dataset: MultiOmicsDataset, config: PipelineConfig) -> dict[str, object]:
    n_runs = int(config.evaluation["repeated_seed_runs"])
    threshold = float(config.evaluation["stability_correlation_threshold"])
    base_seed = int(config.analysis["random_seed"])

    run_paths = [
        train_mofa_run(dataset, config, run_id=f"stability_seed_{i:02d}", seed=base_seed + i - 1)
        for i in range(1, n_runs + 1)
    ]
    models = [load_mofa_model(path) for path in run_paths]
    variance_df = pd.concat(
        [summarise_variance_explained(model, f"seed_{idx:02d}") for idx, model in enumerate(models, start=1)],
        ignore_index=True,
    )

    factor_tables = [model.get_factors(df=True) for model in models]
    factor_cols = [col for col in factor_tables[0].columns if str(col).lower().startswith("factor")]
    reference = factor_tables[0]
    replication_rows: list[dict[str, object]] = []
    for factor_name in factor_cols:
        ref_vec = reference[factor_name].to_numpy()
        corrs = []
        for table in factor_tables[1:]:
            current_corrs = [abs(safe_corrcoef(ref_vec, table[col].to_numpy())) for col in factor_cols]
            current_corrs = [corr for corr in current_corrs if not np.isnan(corr)]
            corrs.append(max(current_corrs) if current_corrs else np.nan)
        replication_rows.append(
            {
                "factor": factor_name,
                "mean_abs_correlation": float(np.nanmean(corrs)),
                "replication_rate": float(np.nanmean(np.array(corrs) >= threshold)),
            }
        )

    return {
        "run_paths": run_paths,
        "variance": variance_df,
        "replication": pd.DataFrame.from_records(replication_rows),
    }


def sample_size_sensitivity_analysis(dataset: MultiOmicsDataset, config: PipelineConfig) -> pd.DataFrame:
    fractions = list(config.evaluation["sample_size_grid"])
    n_reps = int(config.evaluation["replicates_per_size"])
    base_seed = int(config.analysis["random_seed"])
    sample_ids = dataset.sample_metadata["sample_id"].tolist()

    results: list[pd.DataFrame] = []
    for frac in fractions:
        target_n = max(3, int(np.floor(len(sample_ids) * float(frac))))
        for rep_id in range(1, n_reps + 1):
            rng = np.random.default_rng(base_seed + rep_id + target_n)
            chosen_ids = sorted(rng.choice(sample_ids, size=target_n, replace=False).tolist())
            subset_views = subset_views_by_samples(dataset.views, chosen_ids)
            subset_meta = (
                dataset.sample_metadata.loc[dataset.sample_metadata["sample_id"].isin(chosen_ids)]
                .set_index("sample_id")
                .loc[chosen_ids]
                .reset_index()
            )
            subset_dataset = MultiOmicsDataset(
                views=subset_views,
                sample_metadata=subset_meta,
                config_path=dataset.config_path,
            )
            run_path = train_mofa_run(
                subset_dataset,
                config,
                run_id=f"subset_{target_n:02d}_frac_{str(frac).replace('.', '_')}_rep_{rep_id:02d}",
                seed=base_seed + rep_id,
            )
            result = summarise_variance_explained(load_mofa_model(run_path), f"n_{target_n:02d}_rep_{rep_id:02d}")
            result["sample_fraction"] = float(frac)
            result["n_samples"] = target_n
            result["replicate"] = rep_id
            results.append(result)
    return pd.concat(results, ignore_index=True)


def _save(fig: plt.Figure, path: Path) -> Path:
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_seed_stability(seed_results: dict[str, object], path: Path) -> Path:
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.lineplot(data=seed_results["variance"], x="model", y="r2", hue="view", marker="o", ax=ax)
    ax.set_title("Repeated-seed stability of variance explained")
    ax.tick_params(axis="x", rotation=90)
    fig.tight_layout()
    return _save(fig, path)


def plot_factor_replication(seed_results: dict[str, object], path: Path) -> Path:
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.barplot(data=seed_results["replication"], x="factor", y="replication_rate", color="#54A24B", ax=ax)
    ax.axhline(0.7, linestyle="--", color="#E45756")
    ax.set_title("Factor replication across repeated seeds")
    fig.tight_layout()
    return _save(fig, path)


def plot_sample_sensitivity(sensitivity_df: pd.DataFrame, path: Path) -> Path:
    summary = (
        sensitivity_df.groupby(["view", "sample_fraction", "n_samples"], as_index=False)["r2"]
        .agg(mean_r2="mean", sd_r2="std")
        .fillna({"sd_r2": 0.0})
    )
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.lineplot(data=summary, x="n_samples", y="mean_r2", hue="view", marker="o", ax=ax)
    for view_name, group in summary.groupby("view"):
        ax.fill_between(group["n_samples"], group["mean_r2"] - group["sd_r2"], group["mean_r2"] + group["sd_r2"], alpha=0.15)
    ax.set_title("Sample-size sensitivity analysis")
    fig.tight_layout()
    return _save(fig, path)


def run_evaluation_workflow(dataset: MultiOmicsDataset, config: PipelineConfig) -> list[Path]:
    eval_dir = ensure_dir(resolve_path(config, config.paths["evaluation_dir"]))
    plot_dir = ensure_dir(resolve_path(config, config.paths["plot_evaluation_dir"]))

    seed_results = seed_stability_analysis(dataset, config)
    sensitivity_df = sample_size_sensitivity_analysis(dataset, config)

    seed_variance_path = eval_dir / "repeated_seed_variance.csv"
    seed_replication_path = eval_dir / "repeated_seed_replication.csv"
    sensitivity_path = eval_dir / "sample_size_sensitivity.csv"

    seed_results["variance"].to_csv(seed_variance_path, index=False)
    seed_results["replication"].to_csv(seed_replication_path, index=False)
    sensitivity_df.to_csv(sensitivity_path, index=False)

    return [
        seed_variance_path,
        seed_replication_path,
        sensitivity_path,
        plot_seed_stability(seed_results, plot_dir / "repeated_seed_variance.png"),
        plot_factor_replication(seed_results, plot_dir / "factor_replication.png"),
        plot_sample_sensitivity(sensitivity_df, plot_dir / "sample_size_sensitivity.png"),
    ]

