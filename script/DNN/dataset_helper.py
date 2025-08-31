import numpy as np
import torch
from torch.utils.data import TensorDataset

import copy

'''
This is a model-dependent file
Including functions for dealing with loading and processing of datasets
'''

##################################################
# to load dataset for attention like-structure
def loadDS(filePath, withInvMass=1, onlyMatch=0, onlySS=0, randomizeEntry=[]):
    dataset = np.load(filePath)

    if onlyMatch:
        dataset = dataset[dataset[:,-3]>=0]
    if onlySS:
        dataset = dataset[dataset[:,-1]==1]

    # for evaluation permutation performance
    if len(randomizeEntry) != 0:
        print("load data permutation on ", randomizeEntry)
        nEvents, shapeDim= dataset.shape
        perm = np.random.permutation(nEvents)
        for rIdx in randomizeEntry:
            # DO NOT randomize labels !!!
            if rIdx >= shapeDim-3:
                print("load error: The permutation entry is larger than the shape dimension")
                exit(1)
            dataset[:, rIdx] = dataset[:, rIdx][perm]

    # 0-19 4bjet with (px, py, pz, m, score) or (pt, eta, phi, m, score)
    bjetPart = dataset[:, :20]
    x_node = np.reshape(bjetPart, [bjetPart.shape[0], 4, 5])

    # 22-27 dr, 30-36 invmass
    if withInvMass:
        drinvm = np.concatenate([dataset[:, 22:28],dataset[:, 30:36]], axis=1)
        drinvm = np.transpose(drinvm.reshape([drinvm.shape[0], 2,6]), axes=(0,2,1))
        drinvm = np.concatenate([drinvm, drinvm], axis=1)
        x_edge = drinvm
    else:
        dr = dataset[:, 22:28]
        dr = np.concatenate([dr, dr], axis=1)
        dr = dr.reshape((dr.shape[0], 12,1))
        x_edge = dr
    
    # 20-21 njets, 28-29 HT
    x_global = np.concatenate([dataset[:, 20:22], dataset[:, 28:30]], axis=1)

    # last three entries label
    # [bjet_category, n_Muon, pass_same_sign]
    all_labels = dataset[:,-3:]
    all_labels = all_labels.astype(int)

    return (x_node, x_edge, x_global, all_labels)


##################################################
# semi-general function
# to permute single or several entires inside one data randomly, return a new dataset
# @ x_structure can be anything, in b-jet assignment code could be x_node, x_edge or x_global
# @ randomizeEntry each entry should pointing out the structure index and var index
#   e.g. [ 1 ] when input x_global would be the permutation of n_jets
#   e.g. [ (0,1), (1,1), (2,1), (3,1)] when input x_node would be the permutation of all bjet eta together
def permuteDS(x_structure, randomizeEntry=[]):
    permed_structure = np.array(copy.deepcopy(x_structure))
    if len(randomizeEntry) != 0:
        print("load data permutation on ", randomizeEntry)
        # struc_shape = [sub_struc.shape for sub_struc in x_structure]
        # nEvents = struc_shape[0].shape[0]
        nEvents = permed_structure.shape[0]
        perm = np.random.permutation(nEvents)
        for index in randomizeEntry:
            if isinstance(index, int):
                if index >= permed_structure.shape[1]:
                    print("[ERROR] index overflow when permuting, please check!")
                    exit(1)
                permed_structure[:, index] = permed_structure[:, index][perm]
            elif isinstance(index, list) or isinstance(index, tuple):
                if len(index) != 2:
                    print("[ERROR] currently not back up more than 2-dimensional structure yet.")
                    exit(1)
                if len(permed_structure.shape) != 3:
                    print("[ERROR] to perm structure & indexing dimension not corresponding")
                    exit(1)
                if index[0] >= permed_structure.shape[1] or index[1] >= permed_structure.shape[2]:
                    print("[ERROR] index overflow when permuting, please check!")
                    exit(1)
                permed_structure[:, index[0], index[1]] = permed_structure[:, index[0], index[1]][perm]
        return permed_structure
    else:
        return permed_structure

##################################################
# to form tensor dataset based on the structure of the input for model
def formDS(needGlobalVar, needEdgeVar, x_node, x_edge=[], x_global=[], y=[]):
    if needGlobalVar:
        if needEdgeVar:
            ptDS = TensorDataset(
                torch.from_numpy(x_node).float(), 
                torch.from_numpy(x_edge).float(), 
                torch.from_numpy(x_global).float(), 
                torch.from_numpy(y[:, 0]).long())
            def unpack(batch, device=torch.device("cpu")):
                x, x_edge, x_global, y = batch
                return x.to(device), x_edge.to(device), x_global.to(device), y.to(device)
        else:
            ptDS = TensorDataset(
                torch.from_numpy(x_node).float(), 
                torch.from_numpy(x_global).float(), 
                torch.from_numpy(y[:, 0]).long())
            def unpack(batch, device=torch.device("cpu")):
                x, x_global, y = batch
                return x.to(device), x_global.to(device), y.to(device)
    else:
        ptDS = TensorDataset(
            torch.from_numpy(x_node).float(), 
            torch.from_numpy(y[:, 0]).long())
        def unpack(batch, device=torch.device("cpu")):
            x, y = batch
            return x.to(device), y.to(device)
    return ptDS, unpack

##################################################
# to balance the input weight due to the input dimensions
def balanceTrainType(y_label, doTrain=1):
    # in total 6 types for train dataset
    totalNum = y_label.shape[0]
    print("total number:", totalNum)
    eventCount = np.zeros(6)
    for index in range(6):
        eventCount[index] = y_label[y_label[:, 0] == index].shape[0]
        print("type ", index, "events", eventCount[index])

    # for training, keep input dataset balanced
    if doTrain:
        weights = totalNum / (eventCount*6)
    # for testing, providing the contribution for evaluating total performance
    else:
        weights = eventCount / totalNum
    return weights
