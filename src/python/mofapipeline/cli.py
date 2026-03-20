from __future__ import annotations

import argparse
import logging

LOGGER = logging.getLogger(__name__)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Reusable MOFA+ pipeline for multi-omics datasets.")
    parser.add_argument("--config", default="config/analysis_python.yaml", help="Path to pipeline YAML config.")
    parser.add_argument(
        "--stage",
        default="all",
        choices=["qc", "eda", "train", "evaluate", "all"],
        help="Pipeline stage to run.",
    )
    parser.add_argument("--log-file", default=None, help="Optional log file path.")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    from .config import read_config, resolve_path
    from .eda import run_eda_workflow
    from .evaluation import run_evaluation_workflow
    from .io import load_multiomics_data
    from .mofa import run_downstream_workflow, train_mofa_run
    from .preprocessing import build_analysis_datasets, run_qc_workflow
    from .utils import configure_logging

    config = read_config(args.config)
    log_path = resolve_path(config, args.log_file) if args.log_file else None
    configure_logging(log_path)

    dataset = load_multiomics_data(config)
    run_qc_workflow(dataset, config)
    if args.stage == "qc":
        LOGGER.info("QC stage complete.")
        return

    base_dataset, filtered_dataset = build_analysis_datasets(dataset, config)

    if args.stage in {"eda", "all"}:
        LOGGER.info("Running EDA workflow.")
        run_eda_workflow(base_dataset, config)

    if args.stage in {"train", "all"}:
        LOGGER.info("Running MOFA training workflow.")
        full_run = train_mofa_run(base_dataset, config, run_id="run_01_full")
        filtered_run = train_mofa_run(filtered_dataset, config, run_id="run_02_filtered_proteins")
        run_downstream_workflow(base_dataset, config, full_run, filtered_run)

    if args.stage in {"evaluate", "all"}:
        LOGGER.info("Running evaluation workflow.")
        run_evaluation_workflow(filtered_dataset, config)


if __name__ == "__main__":
    main()
