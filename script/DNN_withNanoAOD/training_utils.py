import numpy as np
from matplotlib import pyplot as plt


def balanceTrainType(y_label, doTrain=1):
    y_array = np.asarray(y_label)
    if y_array.ndim == 2:
        labels = y_array[:, 0].astype(int)
    else:
        labels = y_array.astype(int)

    total_num = labels.shape[0]
    print("total number:", total_num)
    n_classes = int(labels.max()) + 1
    event_count = np.zeros(n_classes, dtype=np.float64)
    for index in range(n_classes):
        event_count[index] = np.sum(labels == index)
        print("type ", index, "events", event_count[index])

    if doTrain:
        with np.errstate(divide="ignore", invalid="ignore"):
            weights = np.divide(
                total_num,
                event_count * n_classes,
                out=np.zeros_like(event_count),
                where=event_count > 0,
            )
    else:
        weights = event_count / max(total_num, 1)
    return weights


def plot_confusion_matrix_fractional(
    y_true,
    y_pred,
    suffix,
    trueLabels,
    predLabels,
    output_name_prefix="confusion_matrix",
):
    y_true = np.asarray(y_true, dtype=int)
    y_pred = np.asarray(y_pred, dtype=int)
    n_classes_true = len(trueLabels)
    n_classes_pred = len(predLabels)

    matrix = np.zeros((n_classes_pred, n_classes_true), dtype=np.float32)
    for truth_label, pred_label in zip(y_true, y_pred):
        if 0 <= truth_label < n_classes_true and 0 <= pred_label < n_classes_pred:
            matrix[pred_label, truth_label] += 1

    true_type_total = matrix.sum(axis=0, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        matrix_norm = np.divide(matrix, true_type_total, where=true_type_total != 0)

    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(matrix_norm, cmap="Blues")
    ax.set_xticks(np.arange(n_classes_true))
    ax.set_yticks(np.arange(n_classes_pred))
    ax.set_xticklabels(trueLabels)
    ax.set_yticklabels(predLabels)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

    for pred_index in range(n_classes_pred):
        for truth_index in range(n_classes_true):
            value = matrix_norm[pred_index, truth_index]
            ax.text(
                truth_index,
                pred_index,
                f"{value:.2f}",
                ha="center",
                va="center",
                color="black" if value < 0.6 else "white",
            )

    ax.set_xlabel("True Label")
    ax.set_ylabel("Predicted Label")
    fig.tight_layout()
    plt.colorbar(im, ax=ax)
    plt.savefig(f"{output_name_prefix}{suffix}.png")
    plt.close(fig)


def plot_loss_curves(train_losses, val_losses, output_name):
    fig, ax = plt.subplots()
    ax.plot(train_losses, label="Train Loss")
    ax.plot(val_losses, label="Val Loss")
    ax.set_xlabel("Epoch")
    ax.set_ylabel("Loss")
    ax.legend()
    ax.grid()
    fig.tight_layout()
    plt.savefig(output_name)
    plt.close(fig)
