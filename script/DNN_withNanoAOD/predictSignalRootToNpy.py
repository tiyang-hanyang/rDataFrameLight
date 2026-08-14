import argparse
from pathlib import Path

import numpy as np
import torch
import uproot

from feature_grouping import select_active_features
from signalExtractionInputFromNanoAOD import extract_event_row, load_arrays
from signal_class_config import DEFAULT_SCHEMA_NAME, get_signal_schema, schema_help_text
from simpleSignalModel import FEATURE_DIM, SimpleMLP, load_model_checkpoint


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run simple signal-classification inference on a ROOT file."
    )
    parser.add_argument("--input", required=True, help="Input ROOT file.")
    parser.add_argument("--model", required=True, help="Trained .pt model.")
    parser.add_argument("--output", required=True, help="Output structured .npy file.")
    parser.add_argument("--tree-name", default="Events", help="Tree name.")
    parser.add_argument("--batch-size", type=int, default=2048, help="Inference batch size.")
    parser.add_argument("--device", default=None, help="cpu or cuda")
    parser.add_argument("--jet-pt-shift-branch", default=None, help="JES shift branch, e.g. CMS_scale_j_FlavorQCD.")
    parser.add_argument("--jet-pt-shift-direction", choices=["up", "down"], default="up")
    parser.add_argument("--jer-direction", choices=["up", "down"], default=None, help="Use Jet_JER_corr_{up/down} to vary Jet_pt.")
    parser.add_argument(
        "--muon-pt-branch",
        default=None,
        help="Optional override for the muon pt branch used by the DNN feature builder.",
    )
    parser.add_argument("--jet-pt-threshold", type=float, default=30.0)
    parser.add_argument("--btag-threshold", type=float, default=None)
    parser.add_argument("--scale-jet-mass-with-pt", action="store_true")
    parser.add_argument(
        "--class-schema",
        default=None,
        help="Optional legacy class schema override. If omitted, class names are read from the checkpoint.",
    )
    return parser.parse_args()


def resolve_device(device_name):
    if device_name is not None:
        return torch.device(device_name)
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")


def build_model(model_path, device):
    state_dict, contract = load_model_checkpoint(model_path, device)
    class_names = list(contract["task_label_names"])
    input_dim = len(contract["active_feature_names"])
    if contract["feature_layout"] == "masked":
        input_dim = FEATURE_DIM
    model = SimpleMLP(input_dim=input_dim, n_classes=len(class_names)).to(device)
    model.load_state_dict(state_dict)
    model.eval()
    return model, contract


def resolve_tree_name(input_path, requested_tree_name):
    with uproot.open(input_path) as root_file:
        if requested_tree_name in root_file:
            return requested_tree_name

        tree_candidates = []
        for key in root_file.keys(cycle=False):
            try:
                obj = root_file[key]
            except Exception:
                continue
            if hasattr(obj, "num_entries") and hasattr(obj, "arrays"):
                tree_candidates.append(key)

    if len(tree_candidates) == 1:
        return tree_candidates[0]
    if not tree_candidates:
        raise RuntimeError(
            f"No TTree-like object found in {input_path}. Available keys do not include '{requested_tree_name}'."
        )
    raise RuntimeError(
        f"Requested tree '{requested_tree_name}' not found in {input_path}. "
        f"Available TTrees: {', '.join(tree_candidates)}"
    )


def extract_features_with_ids(
    input_path,
    tree_name,
    jet_pt_shift_branch=None,
    jet_pt_shift_direction="up",
    jer_direction=None,
    muon_pt_branch=None,
    jet_pt_threshold=30.0,
    btag_threshold=None,
    scale_jet_mass_with_pt=False,
):
    from signalExtractionInputFromNanoAOD import apply_jet_pt_variation

    arrays = load_arrays(
        input_path,
        tree_name=tree_name,
        jet_pt_shift_branch=jet_pt_shift_branch,
        jer_direction=jer_direction,
        muon_pt_branch=muon_pt_branch,
    )
    arrays = apply_jet_pt_variation(
        arrays,
        jet_pt_shift_branch=jet_pt_shift_branch,
        jet_pt_shift_direction=jet_pt_shift_direction,
        jer_direction=jer_direction,
        jet_pt_threshold=jet_pt_threshold,
        btag_threshold=btag_threshold,
        scale_jet_mass_with_pt=scale_jet_mass_with_pt,
    )
    ids = []
    features = []
    for event_index in range(len(arrays["Jet_pt"])):
        row = extract_event_row(arrays, event_index, event_weight=np.float32(1.0))
        if row is None:
            continue
        ids.append(
            (
                int(arrays["run"][event_index]),
                int(arrays["luminosityBlock"][event_index]),
                int(arrays["event"][event_index]),
            )
        )
        features.append(row[:FEATURE_DIM])

    if not features:
        return (
            np.empty((0, 3), dtype=np.int64),
            np.empty((0, FEATURE_DIM), dtype=np.float32),
        )

    return np.asarray(ids, dtype=np.int64), np.asarray(features, dtype=np.float32)


def run_inference(model, features, batch_size, device, n_classes):
    scores = []
    with torch.no_grad():
        for start in range(0, features.shape[0], batch_size):
            stop = start + batch_size
            batch = torch.from_numpy(features[start:stop]).to(device)
            logits = model(batch)
            scores.append(torch.softmax(logits, dim=1).cpu().numpy())

    if not scores:
        return np.empty((0, n_classes), dtype=np.float32)
    return np.concatenate(scores, axis=0).astype(np.float32)


def save_output(output_path, ids, scores, class_names):
    prediction = np.argmax(scores, axis=1).astype(np.int32) if scores.size > 0 else np.empty((0,), dtype=np.int32)
    best_score = np.max(scores, axis=1).astype(np.float32) if scores.size > 0 else np.empty((0,), dtype=np.float32)
    predicted_class_name = np.asarray(
        [class_names[index] for index in prediction],
        dtype=f"<U{max(len(name) for name in class_names) if class_names else 1}",
    )
    dtype = [
        ("run", np.int64),
        ("luminosityBlock", np.int64),
        ("event", np.int64),
        ("prediction", np.int32),
        ("best_score", np.float32),
        ("predicted_class_name", predicted_class_name.dtype),
    ]
    dtype.extend((class_name, np.float32) for class_name in class_names)
    output = np.empty(ids.shape[0], dtype=dtype)
    output["run"] = ids[:, 0]
    output["luminosityBlock"] = ids[:, 1]
    output["event"] = ids[:, 2]
    output["prediction"] = prediction
    output["best_score"] = best_score
    output["predicted_class_name"] = predicted_class_name
    for class_index, class_name in enumerate(class_names):
        output[class_name] = scores[:, class_index]
    np.save(output_path, output)


def main():
    args = parse_args()
    device = resolve_device(args.device)
    resolved_tree_name = resolve_tree_name(args.input, args.tree_name)
    ids, features = extract_features_with_ids(
        args.input,
        resolved_tree_name,
        jet_pt_shift_branch=args.jet_pt_shift_branch,
        jet_pt_shift_direction=args.jet_pt_shift_direction,
        jer_direction=args.jer_direction,
        muon_pt_branch=args.muon_pt_branch,
        jet_pt_threshold=args.jet_pt_threshold,
        btag_threshold=args.btag_threshold,
        scale_jet_mass_with_pt=args.scale_jet_mass_with_pt,
    )
    model, contract = build_model(args.model, device)
    class_names = list(contract["task_label_names"])
    if args.class_schema is not None:
        requested_class_names = list(get_signal_schema(args.class_schema).class_names)
        if requested_class_names != class_names:
            raise RuntimeError(
                f"Requested schema '{args.class_schema}' resolves to {requested_class_names}, "
                f"but checkpoint task labels are {class_names}."
            )
    if contract["feature_layout"] == "compact":
        features = select_active_features(features, contract["active_feature_groups"])
    scores = run_inference(model, features, args.batch_size, device, len(class_names))
    save_output(args.output, ids, scores, class_names)
    print(f"tree name: {resolved_tree_name}")
    print(f"processed events: {features.shape[0]}")
    print(f"saved output: {Path(args.output)}")


if __name__ == "__main__":
    main()
