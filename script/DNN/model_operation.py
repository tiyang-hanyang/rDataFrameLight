import torch
import torch.nn as nn
from torch.utils.data import DataLoader
from dataclasses import dataclass, replace, fields
from typing import Any, Dict, Union
from pathlib import Path
import json
import numpy as np
from matplotlib import pyplot as plt

from model_helper import count_parameters
from model_helper import plot_loss

from dataset_helper import permuteDS, formDS

##################################################
# train model
##################################################
@dataclass
class TrainCfg:
    device: str = "cpu"
    num_epochs: int = 100
    num_workers: int = 4
    lossFunc: str = "CrossEntropy"
    learning_rate: float = 1e-3

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "TrainCfg":
        allowed = {f.name: f for f in fields(cls)}
        kwargs: Dict[str, Any] = {}
        for k, v in d.items():
            if k not in allowed:
                continue 
            kwargs[k] = v
        return replace(cls(), **kwargs)
    
    @classmethod
    def from_json(cls, src: Union[str, Path]) -> "TrainCfg":
        src = Path(src)
        if src.exists():
            with src.open("r", encoding="utf-8") as f:
                return cls.from_dict(json.load(f))
        return cls.from_dict(json.loads(str(src)))

def categorized_train(model, trainDS, valDS, unpackingFunc, setup, batch_size, modelName, class_weights = None):
    # setup device and model
    if setup.device == "cuda":
        if torch.cuda.is_available():
            device = torch.device("cuda")
            needPin = True
        else:
            print("[WARNING] cuda not available, using CPU for training")
            device = torch.device("cpu")
            needPin = False
    else:
        device = torch.device("cpu")
        needPin = False
    model = model.to(device)
    count_parameters(model)

    # setup dataset loader
    train_loader = DataLoader(
        trainDS, 
        batch_size=batch_size, 
        num_workers=setup.num_workers, 
        pin_memory=needPin, 
        shuffle=True
    )
    val_loader = DataLoader(
        valDS, 
        batch_size=batch_size, 
        num_workers=setup.num_workers, 
        pin_memory=needPin, 
        shuffle=False
    )
    if class_weights is not None:
        class_weights_tensor = torch.tensor(class_weights, dtype=torch.float32).to(device)

    # prepare for training
    train_losses = []
    val_losses = []
    if setup.lossFunc == "CrossEntropy":
        if class_weights is not None:
            criterion = nn.CrossEntropyLoss(weight=class_weights_tensor)
        else:
            criterion = nn.CrossEntropyLoss()
    else:
        print("[ERROR] Other loss functions not set yet, please check")
    optimizer = torch.optim.Adam(model.parameters(), lr=setup.learning_rate)

    # train
    for epoch in range(setup.num_epochs):
        model.train()
        total_loss = 0
        for batch in train_loader:
            optimizer.zero_grad(set_to_none=True)
            batchUPTensor = unpackingFunc(batch, device)
            out = model(*(batchUPTensor[:-1]))
            loss = criterion(out, batchUPTensor[-1])
            loss.backward()
            optimizer.step()
            total_loss += loss.item()
        train_losses.append(total_loss / len(train_loader))

        model.eval()
        val_loss = 0
        correct = 0
        with torch.no_grad():
            for batch in val_loader:
                batchUPTensor = unpackingFunc(batch, device)
                out = model(*(batchUPTensor[:-1]))
                loss = criterion(out, batchUPTensor[-1])
                val_loss += loss.item()
                pred = out.argmax(dim=1)
                correct += (pred == batchUPTensor[-1]).sum().item()
        val_losses.append(val_loss / len(val_loader))
        val_acc = correct / len(valDS)

        print(f"Epoch {epoch+1}/{setup.num_epochs} | Train Loss: {train_losses[-1]:.4f} | "
            f"Val Loss: {val_losses[-1]:.4f} | Val Acc: {val_acc:.4f}")
        
    torch.save(model.state_dict(), modelName+".pt")
    plot_loss(train_losses, val_losses, modelName)

##################################################
# test the trained model
##################################################
def categorized_test(
        model, 
        x_node_test, x_edge_test, x_global_test, y_test, 
        node_var_list, node_merge_var_list, edge_var_list, global_var_list,
        batch_size, suffix, 
        needGlobalVar=True, needEdgeVar=True, doPerm=True, doGrad=True, doConf=True):

    # however will need, put it here
    nEvents = x_node_test.shape[0]
    testDS, testUnpack = formDS(needGlobalVar, needEdgeVar, x_node_test, x_edge_test, x_global_test, y_test)
    prediction, truthLabel = categorized_evaluation(model, testDS, testUnpack, batch_size)
    metric_performance = np.sum(prediction==truthLabel) / nEvents

    # permutation drop
    if doPerm:
        perm_loss_evaluation(
            x_node_test, x_edge_test, x_global_test, y_test, 
            node_var_list, node_merge_var_list, edge_var_list, global_var_list,
            model, batch_size, metric_performance=metric_performance,
            outputName="permutation_acc_drop"+suffix,  needGlobalVar=needGlobalVar, needEdgeVar=needEdgeVar)

    if doGrad:
        gradient_evaluation(
            model, testDS, testUnpack, 
            batch_size, node_var_list, edge_var_list, global_var_list,
            outputName="gradient_imp"+suffix, needEdgeVar=needEdgeVar, needGlobalVar=needGlobalVar)

    if doConf:
        plot_confusion_matrix_fractional(truthLabel, prediction, "confusion_matrix"+suffix)

##################################################
# function to evaluate the model performance
def categorized_evaluation(model, testDS, unpackingFunc, batch_size):
    test_loader = DataLoader(testDS, batch_size=batch_size, shuffle=False)
    prediction = []
    truthLabel = []
    for batch in test_loader:
        batchUPTensor = unpackingFunc(batch) 
        out = model(*(batchUPTensor[:-1]))
        pred = out.argmax(dim=1)
        prediction.extend(pred.cpu().numpy().tolist())
        truthLabel.extend(batchUPTensor[-1].cpu().numpy().tolist())
    prediction = np.array(prediction, dtype=int)
    truthLabel = np.array(truthLabel, dtype=int)
    return prediction, truthLabel

##################################################
# estimating the permutation loss
def perm_loss_evaluation(
        x_node_test, x_edge_test, x_global_test, y_test, 
        node_var_list, node_merge_var_list, edge_var_list, global_var_list,
        trainedModel, 
        batch_size, outputName, metric_performance,
        needGlobalVar=1, needEdgeVar=1):
    # all this structure would need x_node input
    nEvents = x_node_test.shape[0]
    x_node_loss = []
    nNodes = x_node_test.shape[1]
    nVars = x_node_test.shape[2]
    for n_idx in range(nNodes):
        for v_idx in range(nVars):
            x_node_perm = permuteDS(x_node_test, randomizeEntry=[[n_idx, v_idx]])
            testDS_perm, testUnpack = formDS(needGlobalVar, needEdgeVar, x_node_perm, x_edge_test, x_global_test, y_test)
            prediction_perm, truthLabel_perm = categorized_evaluation(trainedModel, testDS_perm, testUnpack, batch_size)
            permed_performance = np.sum(prediction_perm==truthLabel_perm) / nEvents
            x_node_loss.append((metric_performance-permed_performance)/metric_performance)
    x_node_loss = np.array(x_node_loss, dtype=float)
    # prepare plotting
    permPerformanceDrop = [x_node_loss]
    varList = node_var_list

    x_collective_node_loss = []
    for v_idx in range(nVars):
        collectiveEntries =[ [n_idx, v_idx] for n_idx in range(nNodes)]
        x_node_perm = permuteDS(x_node_test, randomizeEntry=collectiveEntries)
        testDS_perm, testUnpack = formDS(needGlobalVar, needEdgeVar, x_node_perm, x_edge_test, x_global_test, y_test)
        prediction_perm, truthLabel_perm = categorized_evaluation(trainedModel, testDS_perm, testUnpack, batch_size)
        permed_performance = np.sum(prediction_perm==truthLabel_perm) / nEvents
        x_collective_node_loss.append((metric_performance-permed_performance)/metric_performance)
    x_collective_node_loss = np.array(x_collective_node_loss, dtype=float)
    # prepare plotting
    collPermPerformanceDrop = [x_collective_node_loss]
    collVarList = node_merge_var_list

    # when need edge input
    if needEdgeVar:
        x_edge_loss = []
        nEdges = x_edge_test.shape[1]
        nVars = x_edge_test.shape[2]
        for n_idx in range(nEdges):
            for v_idx in range(nVars):
                x_edge_perm = permuteDS(x_edge_test, randomizeEntry=[[n_idx, v_idx]])
                testDS_perm, testUnpack = formDS(needGlobalVar, needEdgeVar, x_node_test, x_edge_perm, x_global_test, y_test) # order determines here
                prediction_perm, truthLabel_perm = categorized_evaluation(trainedModel, testDS_perm, testUnpack, batch_size)
                permed_performance = np.sum(prediction_perm==truthLabel_perm) / nEvents
                x_edge_loss.append((metric_performance-permed_performance)/metric_performance)
        permPerformanceDrop.append(x_edge_loss)
        collPermPerformanceDrop.append(x_edge_loss)
        varList.extend(edge_var_list)
        collVarList.extend(edge_var_list)

    # when need global input
    if needGlobalVar:
        x_global_loss = []
        nVars = x_global_test.shape[1]
        for v_idx in range(nVars):
            x_global_perm = permuteDS(x_global_test, randomizeEntry=[v_idx])
            testDS_perm, testUnpack = formDS(needGlobalVar, needEdgeVar, x_node_test, x_edge_test, x_global_perm, y_test) # order determines here
            prediction_perm, truthLabel_perm = categorized_evaluation(trainedModel, testDS_perm, testUnpack, batch_size)
            permed_performance = np.sum(prediction_perm==truthLabel_perm) / nEvents
            x_global_loss.append((metric_performance-permed_performance)/metric_performance)
        permPerformanceDrop.append(x_global_loss)
        collPermPerformanceDrop.append(x_global_loss)
        varList.extend(global_var_list)
        collVarList.extend(global_var_list)

    # merge information
    permPerformanceDrop = np.concatenate(permPerformanceDrop)
    collPermPerformanceDrop = np.concatenate(collPermPerformanceDrop)

    # plot 
    fig1, ax1 = plt.subplots(figsize=(10, 6))
    im1=ax1.bar(np.arange(len(varList)), permPerformanceDrop)
    ax1.set_xticks(np.arange(len(varList)))
    ax1.set_xticklabels(varList)
    # Rotate x labels
    plt.setp(ax1.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    ax1.set_xlabel("Var perm")
    ax1.set_ylabel("Acc drop / %")
    ax1.set_title("")
    fig1.tight_layout()
    plt.savefig(outputName+"_separate.png")
    plt.close()

    fig2, ax2 = plt.subplots(figsize=(10, 6))
    im1=ax2.bar(np.arange(len(collVarList)), collPermPerformanceDrop)
    ax2.set_xticks(np.arange(len(collVarList)))
    ax2.set_xticklabels(collVarList)
    # Rotate x labels
    plt.setp(ax2.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    ax2.set_xlabel("Var perm")
    ax2.set_ylabel("Acc drop / %")
    ax2.set_title("")
    fig2.tight_layout()
    plt.savefig(outputName+"_collective.png")
    plt.close()


##################################################
# estimating the gradient loss
def gradient_evaluation(model, testDS, testUnpack, batch_size, node_var_list, edge_var_list, global_var_list, outputName, needEdgeVar=True, needGlobalVar=True):
    test_loader = DataLoader(testDS, batch_size=batch_size, shuffle=False)
    allBatchImportance = []
    # getting the gradient for each part with the original input structure
    for dataBatch in test_loader:
        gradBatch = testUnpack([comp.detach().clone().requires_grad_(True) for comp in dataBatch[:-1]]+[dataBatch[-1]])
        # get gradients based on the maximum score node score in a batch
        out = model(*(gradBatch[:-1]))
        idx = out.argmax(dim=1)
        batchOut = out[torch.arange(len(gradBatch[0])), idx].sum()
        grads_all = torch.autograd.grad(batchOut, tuple(gradBatch[:-1]), retain_graph=False, create_graph=False)
        # batch average importance from gradients
        meanGrads = []
        for i in range(len(gradBatch)-1):
            meanGrads.append( (grads_all[i] * gradBatch[i]).view(gradBatch[i].shape).abs().mean(dim=0) )
        allBatchImportance.append(meanGrads)

    # cross-batch average & merging
    pos = 0
    totalImp = []
    nodeImp = torch.stack([bGrad[pos] for bGrad in allBatchImportance ], dim=0).mean(dim=0).detach().numpy()
    totalImp.append(nodeImp.reshape([np.prod(nodeImp.shape, dtype=int)]))
    varList = node_var_list
    # for edge processing
    if needEdgeVar:
        pos+=1
        nEdges = testDS[pos].shape[1]
        premeanGradEdge = torch.stack([bGrad[pos] for bGrad in allBatchImportance ], dim=0).mean(dim=0).detach().numpy()
        # batch merging and directional average
        pairs = np.arange(nEdges).reshape([2, int(nEdges/2)]).transpose()
        pairs = torch.tensor(pairs)
        edgeImp = premeanGradEdge[ pairs, : ].mean(axis=1)
        totalImp.append(edgeImp.transpose().reshape([np.prod(edgeImp.shape, dtype=int)]))
        varList.extend(edge_var_list)
    if needGlobalVar:
        pos+=1
        globalImp = torch.stack([bGrad[pos] for bGrad in allBatchImportance ], dim=0).mean(dim=0).detach().numpy()
        totalImp.append(globalImp.reshape([np.prod(globalImp.shape, dtype=int)]))
        varList.extend(global_var_list)
    totalImp = np.concatenate(totalImp)

    fig, ax = plt.subplots(figsize=(10, 6))
    img=ax.bar(np.arange(len(totalImp)), totalImp)
    ax.set_xticks(np.arange(len(totalImp)))
    ax.set_xticklabels(varList)
    # Rotate x labels
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    ax.set_xlabel("Var")
    ax.set_ylabel("gradient importance")
    ax.set_title("")
    fig.tight_layout()
    plt.savefig(outputName+".png")
    plt.close()


##################################################
# function to plot the performance of the dataset
def plot_confusion_matrix_fractional(y_true, y_pred, suffix):
    # collect data
    label_map = {-1: 0, 0:1, 1:2, 2:3, 3:4, 4:5, 5:6}
    print(len(y_true))
    y_true_mapped = np.vectorize(label_map.get)(y_true)
    print(y_true_mapped.shape)
    n_classes_in = 7
    n_classes_out = 6 
    matrix = np.zeros((n_classes_out, n_classes_in), dtype=np.float32)
    for i in range(len(y_true)):
        true_cls = y_true_mapped[i]
        pred_cls = y_pred[i]
        if 0 <= pred_cls < n_classes_out:
            matrix[pred_cls, true_cls] += 1
    print(matrix)

    # Normalize each column
    col_sums = matrix.sum(axis=0, keepdims=True)
    with np.errstate(divide='ignore', invalid='ignore'):
        matrix_norm = np.divide(matrix, col_sums, where=col_sums!=0)

    # Plot
    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(matrix_norm, cmap="Blues")
    ax.set_xticks(np.arange(n_classes_in))
    ax.set_yticks(np.arange(n_classes_out))
    x_labels = ["no matching"] + [f"cate {i}" for i in range(6)]
    y_labels = [f"cate {i}" for i in range(6)]
    ax.set_xticklabels(x_labels)
    ax.set_yticklabels(y_labels)
    # Rotate x labels
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    # Annotate each cell
    for i in range(n_classes_out):
        for j in range(n_classes_in):
            value = matrix_norm[i, j]
            ax.text(j, i, f"{value:.2f}", ha="center", va="center",
                    color="black" if value < 0.6 else "white")
    ax.set_xlabel("True Label")
    ax.set_ylabel("Predicted Label")
    ax.set_title("Confusion Matrix: TTHH DL MINIAOD")
    fig.tight_layout()
    plt.colorbar(im, ax=ax)
    plt.savefig("confusion_matrix_bjet_assignment_GNN"+ suffix+".png")