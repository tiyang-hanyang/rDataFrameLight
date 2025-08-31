from dataclasses import dataclass
import torch
from torch.utils.data import DataLoader

from dataset_helper import loadDS, formDS, balanceTrainType

##################################################
# define this job, whether do Train, do which test, what types of input used
@dataclass
class JobDef:
    modelType: str = "GNNTransformer"
    doTrain: bool = True
    # what to do for the test work
    doPerm: bool = True
    doGrad: bool = True
    doConf: bool = True
    # about the dataset format
    doCartesian: bool = True
    needEdgeVar: bool = True
    withInvMass: bool = True
    needGlobalVar: bool = True
    # dataset information
    channel: str = "TTHH_DL"
    era: str = "Run2023C"

##################################################
# main function
if __name__ == '__main__':
    thisJob = JobDef(modelType="GNNTransformer", 
                     doTrain=0, doPerm=1, doGrad=1, doConf=1,
                     doCartesian=0, needEdgeVar=1, withInvMass=1, needGlobalVar=1, channel="TTHH_DL", era="Run2023C")

    # create the model
    if thisJob.modelType == "GNNTransformer":
        from custom_model import Bjet_assignment_GNN_Transformer
        model = Bjet_assignment_GNN_Transformer(feature_dim=5, edge_dim=2, d_model=32, nhead=8, num_layers=3, n_classes=6, global_dim=4, enc_dropout=0.3, pool_dropout=0.3, use_norm_first=False)
    elif thisJob.modelType == "Transformer":
        from custom_model import Bjet_assignment_Transformer
        model = Bjet_assignment_Transformer(feature_dim=5, d_model=64, nhead=8, num_layers=3, n_classes=6, global_dim=4, enc_dropout=0.3, pool_dropout=0.3, use_norm_first=False)
    elif thisJob.modelType == "GNN":
        from custom_model import Bjet_assignment_GNN
        model = Bjet_assignment_GNN(input_dim=5, edge_dim=2, hidden_dim=32, num_classes=6)
    else:
        print("The model type is not registered, please check!")
        exit(1)

    if thisJob.doCartesian:
        suffix = "_cartesian"
    else:
        suffix = "_polar"
    if thisJob.withInvMass:
        suffix += "_withInvMass"

    # now put it here
    batch_size=256

    if thisJob.doTrain:
        from model_operation import TrainCfg, categorized_train
        trainSetup=TrainCfg.from_json("trainConfig/bjetAssignment.json")
        trainSetup.batch_size = batch_size
        x_node, x_edge, x_global, y_train = loadDS(
            "datasets/"+thisJob.channel+"/"+thisJob.era+"_train"+suffix+".npy", 
            thisJob.withInvMass, 
            onlyMatch=1, onlySS=0)
        x_node_val, x_edge_val, x_global_val, y_val = loadDS(
            "datasets/"+thisJob.channel+"/"+thisJob.era+"_val"+suffix+".npy", 
            thisJob.withInvMass, 
            onlyMatch=1, onlySS=0)

        trainDS, trainUnpack = formDS(thisJob.needGlobalVar, thisJob.needEdgeVar, x_node, x_edge, x_global, y_train)
        valDS, valUnpack = formDS(thisJob.needGlobalVar, thisJob.needEdgeVar, x_node_val, x_edge_val, x_global_val, y_val)
        typeWeight = balanceTrainType(y_train)
        categorized_train(model, trainDS, valDS, trainUnpack, trainSetup, batch_size, thisJob.modelType+suffix, class_weights=typeWeight)
    else:
        from model_operation import categorized_test
        trainedPara = torch.load(thisJob.modelType+suffix+".pt")
        model.load_state_dict(trainedPara)
        x_node_test, x_edge_test, x_global_test, y_test = loadDS(
            "datasets/"+thisJob.channel+"/"+thisJob.era+"_test"+suffix+".npy", 
            thisJob.withInvMass, 
            onlyMatch=0, onlySS=1)

        if thisJob.doCartesian:
            node_var_list = [
                r'$\text{b-jet}_{1}$ $p_{x}$', r'$\text{b-jet}_{1}$ $p_{y}$', r'$\text{b-jet}_{1}$ $p_{z}$', r'$\text{b-jet}_{1}$ m', r'$\text{b-jet}_{1}$ score',
                r'$\text{b-jet}_{2}$ $p_{x}$', r'$\text{b-jet}_{2}$ $p_{y}$', r'$\text{b-jet}_{2}$ $p_{z}$', r'$\text{b-jet}_{2}$ m', r'$\text{b-jet}_{2}$ score',
                r'$\text{b-jet}_{3}$ $p_{x}$', r'$\text{b-jet}_{3}$ $p_{y}$', r'$\text{b-jet}_{3}$ $p_{z}$', r'$\text{b-jet}_{3}$ m', r'$\text{b-jet}_{3}$ score',
                r'$\text{b-jet}_{4}$ $p_{x}$', r'$\text{b-jet}_{4}$ $p_{y}$', r'$\text{b-jet}_{4}$ $p_{z}$', r'$\text{b-jet}_{4}$ m', r'$\text{b-jet}_{4}$ score']
            node_merge_var_list = [r'b-jet $p_{x}$', r'b-jet $p_{y}$', r'b-jet $p_{z}$', r'b-jet m', r'b-jet score']
        else:
            node_var_list = [
                r'$\text{b-jet}_{1}$ $p_{T}$', r'$\text{b-jet}_{1}$ $\eta$', r'$\text{b-jet}_{1}$ $\phi$', r'$\text{b-jet}_{1}$ m', r'$\text{b-jet}_{1}$ score',
                r'$\text{b-jet}_{2}$ $p_{T}$', r'$\text{b-jet}_{2}$ $\eta$', r'$\text{b-jet}_{2}$ $\phi$', r'$\text{b-jet}_{2}$ m', r'$\text{b-jet}_{2}$ score',
                r'$\text{b-jet}_{3}$ $p_{T}$', r'$\text{b-jet}_{3}$ $\eta$', r'$\text{b-jet}_{3}$ $\phi$', r'$\text{b-jet}_{3}$ m', r'$\text{b-jet}_{3}$ score',
                r'$\text{b-jet}_{4}$ $p_{T}$', r'$\text{b-jet}_{4}$ $\eta$', r'$\text{b-jet}_{4}$ $\phi$', r'$\text{b-jet}_{4}$ m', r'$\text{b-jet}_{4}$ score']
            node_merge_var_list = [r'b-jet $p_{T}$', r'b-jet $\eta$', r'b-jet $\phi$', r'b-jet m', r'b-jet score']

        edge_var_list = [
            r'$b_{1,2}$ $\Delta r$',r'$b_{1,3}$ $\Delta r$',r'$b_{1,4}$ $\Delta r$',r'$b_{2,3}$ $\Delta r$',r'$b_{2,4}$ $\Delta r$',r'$b_{3,4}$ $\Delta r$',
            r'$b_{1,2}$ m',r'$b_{1,3}$ m',r'$b_{1,4}$ m',r'$b_{2,3}$ m',r'$b_{2,4}$ m',r'$b_{3,4}$ m'
        ]

        global_var_list = [
            r'$n_{\text{jets}}$', r'$n_{\text{b-jets}}$', r'$H_{T,\text{jets}}$', r'$H_{T,\text{b-jets}}$'
        ]

        categorized_test(
            model, 
            x_node_test, x_edge_test, x_global_test, y_test, 
            node_var_list, node_merge_var_list, edge_var_list, global_var_list,
            batch_size, "_"+thisJob.modelType+suffix, 
            needGlobalVar=thisJob.needGlobalVar, needEdgeVar=thisJob.needEdgeVar, 
            doPerm=thisJob.doPerm, doGrad=thisJob.doGrad, doConf=thisJob.doConf)