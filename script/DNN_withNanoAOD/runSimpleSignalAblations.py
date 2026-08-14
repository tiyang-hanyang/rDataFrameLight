import argparse
import subprocess
import sys
from pathlib import Path


ABLATION_PRESETS = {
    "g1267": ["1", "2", "6", "7"],
    "g24567": ["2", "4", "5", "6", "7"],
    "g13567": ["1", "3", "5", "6", "7"],
    "g3457": ["3", "4", "5", "7"],
    "g34567": ["3", "4", "5", "6", "7"],
    "g12567": ["1", "2", "5", "6", "7"],
    "g123456": ["1", "2", "3", "4", "5", "6"],
    "g126": ["1", "2", "6"],
    "g123457": ["1", "2", "3", "4", "5", "7"],
    "g127": ["1", "2", "7"],
    "g567": ["5", "6", "7"],
    "g56": ["5", "6"],
    "g57": ["5", "7"],
}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Run a batch of simpleSignalModel.py ablation trainings and, optionally, "
            "test/analyze/compare them in one pass."
        )
    )
    parser.add_argument("--train", required=True, help="OS_train.npy")
    parser.add_argument("--val", required=True, help="OS_val.npy")
    parser.add_argument("--test", default=None, help="Optional test dataset for batch evaluation and analysis.")
    parser.add_argument("--dataset-report", default=None, help="mix_report.json")
    parser.add_argument("--output-root", required=True, help="Parent output directory for all ablation jobs.")
    parser.add_argument(
        "--label-subset",
        nargs="+",
        default=None,
        help="Optional sub-task labels, e.g. 0 1. Omit to run the full task, such as four-class training.",
    )
    parser.add_argument(
        "--presets",
        nargs="+",
        default=list(ABLATION_PRESETS.keys()),
        help=f"Ablation preset names. Available: {', '.join(ABLATION_PRESETS)}",
    )
    parser.add_argument("--epochs", type=int, default=100)
    parser.add_argument("--batch-size", type=int, default=1024)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--device", default=None)
    parser.add_argument(
        "--training-weight-mode",
        choices=["balanced", "input", "abs"],
        default="balanced",
    )
    parser.add_argument(
        "--weight-target",
        action="append",
        default=[],
        help="Unified weight interface forwarded to simpleSignalModel.py train.",
    )
    parser.add_argument(
        "--analysis-weight-mode",
        choices=["balanced", "abs", "raw", "unit"],
        default="balanced",
    )
    parser.add_argument("--perm-repeats", type=int, default=5)
    parser.add_argument("--perm-max-events", type=int, default=20000)
    parser.add_argument("--perm-min-events", type=int, default=200)
    parser.add_argument("--gradient-views", choices=["basic", "full"], default="basic")
    parser.add_argument("--score-bins", nargs="+", type=float, default=[0.0, 0.25, 0.5, 0.75, 1.0])
    parser.add_argument(
        "--score-normalization",
        choices=["raw", "global", "shape", "shape_abs"],
        default="shape",
        help="Score-plot normalization passed to simpleSignalModel.py test.",
    )
    parser.add_argument("--top-n", type=int, default=15, help="Top features/groups shown in comparison plots.")
    parser.add_argument("--compare-name", default="comparison", help="Output subdirectory for aggregate comparison.")
    parser.add_argument("--skip-train", action="store_true", help="Reuse existing run directories and skip training.")
    parser.add_argument("--skip-test", action="store_true", help="Skip test-time prediction and score plot generation.")
    parser.add_argument("--skip-analyze", action="store_true", help="Skip feature analysis.")
    parser.add_argument("--skip-compare", action="store_true", help="Skip aggregate comparison products.")
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help=(
            "Skip individual train/test/analyze steps when their expected outputs already exist. "
            "Useful for appending new presets without recomputing finished runs."
        ),
    )
    parser.add_argument("--python-bin", default=sys.executable if sys.executable else "python")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def script_dir():
    return Path(__file__).resolve().parent


def build_run_dir(output_root, preset_name, label_subset):
    label_suffix = "full" if not label_subset else "_".join(label_subset)
    return Path(output_root) / f"{preset_name}_{label_suffix}"


def build_train_command(args, run_dir, feature_groups):
    command = [
        args.python_bin,
        str(script_dir() / "simpleSignalModel.py"),
        "train",
        "--train",
        args.train,
        "--val",
        args.val,
        "--output-dir",
        str(run_dir),
        "--epochs",
        str(args.epochs),
        "--batch-size",
        str(args.batch_size),
        "--lr",
        str(args.lr),
        "--training-weight-mode",
        args.training_weight_mode,
        "--feature-layout",
        "masked",
        "--feature-groups",
        *feature_groups,
    ]
    if args.label_subset:
        command.extend(["--label-subset", *args.label_subset])
    if args.dataset_report is not None:
        command.extend(["--dataset-report", args.dataset_report])
    if args.device is not None:
        command.extend(["--device", args.device])
    for entry in args.weight_target:
        command.extend(["--weight-target", entry])
    return command


def build_test_command(args, run_dir, feature_groups):
    command = [
        args.python_bin,
        str(script_dir() / "simpleSignalModel.py"),
        "test",
        "--test",
        args.test,
        "--model",
        str(run_dir / "model.pt"),
        "--output-dir",
        str(run_dir),
        "--batch-size",
        str(args.batch_size),
        "--feature-layout",
        "masked",
        "--feature-groups",
        *feature_groups,
        "--save-basic-plots",
        "--score-normalization",
        args.score_normalization,
        "--score-bins",
        *[str(edge) for edge in args.score_bins],
    ]
    if args.label_subset:
        command.extend(["--label-subset", *args.label_subset])
    if args.dataset_report is not None:
        command.extend(["--dataset-report", args.dataset_report])
    if args.device is not None:
        command.extend(["--device", args.device])
    return command


def build_analysis_command(args, run_dir):
    command = [
        args.python_bin,
        str(script_dir() / "analyzeSimpleSignalFeatures.py"),
        "--dataset",
        args.test,
        "--model",
        str(run_dir / "model.pt"),
        "--output-dir",
        str(run_dir / "analysis"),
        "--batch-size",
        str(args.batch_size),
        "--weight-mode",
        args.analysis_weight_mode,
        "--perm-repeats",
        str(args.perm_repeats),
        "--perm-max-events",
        str(args.perm_max_events),
        "--perm-min-events",
        str(args.perm_min_events),
        "--gradient-views",
        args.gradient_views,
    ]
    if args.label_subset:
        command.extend(["--label-subset", *args.label_subset])
    if args.dataset_report is not None:
        command.extend(["--dataset-report", args.dataset_report])
    if args.device is not None:
        command.extend(["--device", args.device])
    return command


def build_compare_command(args, run_dirs, preset_names):
    command = [
        args.python_bin,
        str(script_dir() / "compareSimpleSignalRuns.py"),
        "--runs",
        *[str(run_dir) for run_dir in run_dirs],
        "--labels",
        *preset_names,
        "--output-dir",
        str(Path(args.output_root) / args.compare_name),
        "--top-n",
        str(args.top_n),
    ]
    return command


def run_command(command, dry_run):
    print("cmd:", " ".join(command))
    if not dry_run:
        subprocess.run(command, check=True)


def require_existing_file(path, context):
    if not Path(path).is_file():
        raise RuntimeError(f"Missing {context}: {path}")


def train_output_exists(run_dir):
    return (Path(run_dir) / "model.pt").is_file()


def test_output_exists(run_dir):
    return (Path(run_dir) / "predictions.npy").is_file()


def analysis_output_exists(run_dir):
    analysis_dir = Path(run_dir) / "analysis"
    required = [
        analysis_dir / "gradients" / "gradient_importance_summary.csv",
        analysis_dir / "permutation" / "permutation_importance_summary.csv",
        analysis_dir / "group_permutation" / "group_permutation_importance_summary.csv",
    ]
    return all(path.is_file() for path in required)


def main():
    args = parse_args()
    selected_presets = []
    for preset_name in args.presets:
        if preset_name not in ABLATION_PRESETS:
            raise RuntimeError(
                f"Unknown preset '{preset_name}'. Available: {', '.join(sorted(ABLATION_PRESETS))}"
            )
        selected_presets.append((preset_name, ABLATION_PRESETS[preset_name]))

    if args.test is None and (not args.skip_test or not args.skip_analyze or not args.skip_compare):
        raise RuntimeError(
            "--test is required unless you explicitly set --skip-test --skip-analyze --skip-compare."
        )

    run_jobs = []
    for preset_name, feature_groups in selected_presets:
        run_dir = build_run_dir(args.output_root, preset_name, args.label_subset)
        run_jobs.append((preset_name, feature_groups, run_dir))
        print(f"[{preset_name}] output: {run_dir}")

    run_dirs = [run_dir for _, _, run_dir in run_jobs]
    run_labels = [preset_name for preset_name, _, _ in run_jobs]

    for preset_name, feature_groups, run_dir in run_jobs:
        if args.skip_train:
            require_existing_file(run_dir / "model.pt", f"model checkpoint for preset '{preset_name}'")

        if not args.skip_train:
            if args.skip_existing and train_output_exists(run_dir):
                print(f"[{preset_name}] skip train: model.pt already exists")
            else:
                run_command(build_train_command(args, run_dir, feature_groups), args.dry_run)

    for preset_name, feature_groups, run_dir in run_jobs:
        if args.skip_train and not train_output_exists(run_dir) and not args.dry_run:
            raise RuntimeError(f"[{preset_name}] cannot test without existing model.pt in {run_dir}")

        if not args.skip_test:
            if args.skip_existing and test_output_exists(run_dir):
                print(f"[{preset_name}] skip test: predictions.npy already exists")
            else:
                if args.skip_existing and not train_output_exists(run_dir):
                    print(f"[{preset_name}] test requires training output; running test after training/preexisting model")
                run_command(build_test_command(args, run_dir, feature_groups), args.dry_run)

    for preset_name, _, run_dir in run_jobs:
        if args.skip_train and not train_output_exists(run_dir) and not args.dry_run:
            raise RuntimeError(f"[{preset_name}] cannot analyze without existing model.pt in {run_dir}")

        if not args.skip_analyze:
            if args.skip_existing and analysis_output_exists(run_dir):
                print(f"[{preset_name}] skip analyze: analysis summaries already exist")
            else:
                run_command(build_analysis_command(args, run_dir), args.dry_run)

    if not args.skip_compare:
        run_command(build_compare_command(args, run_dirs, run_labels), args.dry_run)


if __name__ == "__main__":
    main()
