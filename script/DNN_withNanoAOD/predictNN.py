import argparse
from pathlib import Path

import numpy as np
import torch

from RNNBJAExtraction_staged import extract_inference_features_from_file
from ttbb_BJA_model import mixingDNN


MODEL_KWARGS = {
    "d_globalVar": 16,
    "d_global_hiddens": [50, 50, 50],
    "d_jet_features": 7,
    "d_LSTM_hiddens": [128, 64, 64, 32, 32],
    "d_jet_out": 20,
    "dropout": 0.1,
}

CLASS_NAMES = ["b1b2", "b1b3", "b1b4", "b2b3", "b2b4", "b3b4"]


def parseConfig():
    parser = argparse.ArgumentParser(
        description="Extract ttbb BJA inference features from NanoAOD-like ROOT files and run the trained model."
    )
    parser.add_argument("--input", required=True, help="Input ROOT file path.")
    parser.add_argument(
        "--model",
        required=True,
        help="Path to the trained ttbb BJA PyTorch weights (.pt).",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output .npy path. Default: <input>_ttbb_BJA_scores.npy",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=1024,
        help="Inference batch size.",
    )
    parser.add_argument(
        "--tree-name",
        default="Events",
        help="ROOT tree name.",
    )
    parser.add_argument(
        "--device",
        default=None,
        choices=["cpu", "cuda", None],
        help="Force inference device. Default chooses cuda when available.",
    )
    return parser.parse_args()


def BJA_var(filePath, tree_name="Events"):
    return extract_inference_features_from_file(
        file_name=filePath,
        tree_name=tree_name,
        include_ids=True,
    )


def split_features(feature_matrix):
    jet_part = feature_matrix[:, :28]
    x_node = np.reshape(jet_part, (jet_part.shape[0], 4, 7)).astype(np.float32)
    x_global = feature_matrix[:, 28:].astype(np.float32)
    return x_node, x_global


def build_model(model_path, device):
    model = mixingDNN(**MODEL_KWARGS).to(device)
    state_dict = torch.load(model_path, map_location=device)
    model.load_state_dict(state_dict)
    model.eval()
    return model


def apply_model(model, dataSet, batch_size=1024, device="cpu"):
    x_node, x_global = split_features(dataSet)
    scores = []

    with torch.no_grad():
        for start in range(0, x_node.shape[0], batch_size):
            stop = start + batch_size
            batch_node = torch.from_numpy(x_node[start:stop]).to(device)
            batch_global = torch.from_numpy(x_global[start:stop]).to(device)
            logits = model(batch_node, batch_global)
            batch_scores = torch.softmax(logits, dim=1).cpu().numpy()
            scores.append(batch_scores)

    if not scores:
        return np.empty((0, len(CLASS_NAMES)), dtype=np.float32), np.empty((0,), dtype=np.int64)

    score_array = np.concatenate(scores, axis=0).astype(np.float32)
    prediction = np.argmax(score_array, axis=1).astype(np.int64)
    return score_array, prediction


def saveScoreWithID(output_path, ids, scores, prediction):
    n_events = scores.shape[0]
    dtype = [
        ("run", np.int64),
        ("luminosityBlock", np.int64),
        ("event", np.int64),
        ("prediction", np.int64),
    ]
    dtype.extend((class_name, np.float32) for class_name in CLASS_NAMES)

    output = np.empty(n_events, dtype=dtype)
    output["run"] = ids["run"]
    output["luminosityBlock"] = ids["luminosityBlock"]
    output["event"] = ids["event"]
    output["prediction"] = prediction

    for class_index, class_name in enumerate(CLASS_NAMES):
        output[class_name] = scores[:, class_index]

    np.save(output_path, output)


def main():
    args = parseConfig()

    input_path = Path(args.input)
    model_path = Path(args.model)
    output_path = (
        Path(args.output)
        if args.output is not None
        else input_path.with_name(f"{input_path.stem}_ttbb_BJA_scores.npy")
    )

    device_name = args.device
    if device_name is None:
        device_name = "cuda" if torch.cuda.is_available() else "cpu"
    device = torch.device(device_name)

    extracted = BJA_var(str(input_path), tree_name=args.tree_name)
    features = extracted["features"]
    ids = {
        "run": extracted["run"].astype(np.int64),
        "luminosityBlock": extracted["luminosityBlock"].astype(np.int64),
        "event": extracted["event"].astype(np.int64),
    }

    model = build_model(str(model_path), device)
    scores, prediction = apply_model(
        model=model,
        dataSet=features,
        batch_size=args.batch_size,
        device=device,
    )

    saveScoreWithID(str(output_path), ids, scores, prediction)
    print(f"processed events: {features.shape[0]}")
    print(f"saved output: {output_path}")


if __name__ == "__main__":
    main()
