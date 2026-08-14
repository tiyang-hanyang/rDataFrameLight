import torch
import torch.nn as nn
from torch.utils.data import TensorDataset, DataLoader

from matplotlib import pyplot as plt
import numpy as np
import copy

from training_utils import balanceTrainType, plot_confusion_matrix_fractional

##################################################
def fixDRComputation(x_node, x_edge):
    # eta computation would be a bit easier: no need of considering the sign of phi
    #   no need of scaling back
    # eta = arcsinh( pz / pt )
    eta = np.asinh(x_node[:, :, 2] / np.sqrt(x_node[:, :, 0]**2 + x_node[:, :, 1]**2))

    # extract information
    order = np.array([[0,1], [0,2], [0,3], [1,2], [1,3], [2,3]])
    deta = eta[:, order[:, 0]] - eta[:, order[:, 1]]
    x_edge_copy = copy.deepcopy(x_edge)
    dr = x_edge_copy[:,:6,0]

    # correct dr
    orig_dPhi2 = dr**2 - deta**2
    orig_dPhi2[orig_dPhi2<0] = 0.0
    orig_dPhi = np.sqrt(orig_dPhi2)
    orig_dPhi[orig_dPhi>np.pi] = 2*np.pi - orig_dPhi[orig_dPhi>np.pi]
    corrected_dr = np.sqrt(deta**2 + orig_dPhi**2)
    corrected_dr = np.tile(corrected_dr,2)

    x_edge_copy[:,:,0] = corrected_dr

    return x_edge_copy

def hardCode_de_std(x_node):
    # x_node in form of [N, 4, 5], for the [0-3] is the three-momentum, just times STD
    # 4 and 5 are mass and score, need to times STD then plus Mean
    # those values are hard-coded here
    x_node_copy = copy.deepcopy(x_node)

    p_std = 146.2224842392
    m_mean = 14.4182656164
    m_std = 10.6231408469
    PNet_mean = 0.7925129098
    PNet_std = 0.2690405989

    x_node_copy[:, :, 0:3] = x_node_copy[:,:,0:3] * p_std
    x_node_copy[:, :, 3] = x_node_copy[:,:,3]*m_std + m_mean
    x_node_copy[:, :, 4] = x_node_copy[:,:,4]*PNet_std + PNet_mean

    return x_node_copy

##################################################
def fixDRComputation_polar(x_node, x_edge):
    eta =  x_node[:, :, 1]
    phi =  x_node[:, :, 2]

    # re-compute dr
    order = np.array([[0,1], [0,2], [0,3], [1,2], [1,3], [2,3]])
    deta = eta[:, order[:, 0]] - eta[:, order[:, 1]]
    dphi = np.abs(phi[:, order[:, 0]] - phi[:, order[:, 1]])
    dphi[dphi>np.pi] = 2*np.pi - dphi[dphi>np.pi]
    corrected_dr = np.sqrt(deta**2 + dphi**2)
    corrected_dr = np.tile(corrected_dr,2)

    # packed back
    x_edge_copy = copy.deepcopy(x_edge)
    x_edge_copy[:,:,0] = corrected_dr

    return x_edge_copy

def hardCode_de_std_polar(x_node):
    # x_node in form of [N, 4, 5], for the [0-3] is the three-momentum, just times STD
    # 4 and 5 are mass and score, need to times STD then plus Mean
    # those values are hard-coded here
    x_node_copy = copy.deepcopy(x_node)

    pt_mean = 97.1714688314
    pt_std = 79.3298290615
    m_mean = 14.4182656164
    m_std = 10.6231408469
    PNet_mean = 0.7925129098
    PNet_std = 0.2690405989

    x_node_copy[:, :, 0] = x_node_copy[:,:,0] * pt_std + pt_mean
    x_node_copy[:, :, 3] = x_node_copy[:,:,3] * m_std + m_mean
    x_node_copy[:, :, 4] = x_node_copy[:,:,4] * PNet_std + PNet_mean

    return x_node_copy


##################################################
# util, for printing the model structure
def count_parameters(model):
    total = 0
    print(f"{'Layer':<40} {'Param #':>10}")
    print("="*52)
    for name, param in model.named_parameters():
        if param.requires_grad:
            num_params = param.numel()
            total += num_params
            print(f"{name:<40} {num_params:>10}")
    print("="*52)
    print(f"{'Total Trainable Params':<40} {total:>10}")

# model structure
# the concept of the additional bjets
# additional jets: jets not originate from the top
# generator-level: from the decay chain
# details in AN2021_040_v7: GenHFHadronMatcher tool,  Jets with ghost-matched b or c hadrons that do not have any parton-level top quarks in their history are referred to as "additional" jets
# This information does not directly exist in the NANOAOD, so we are using the trace-back with the NANOAOD level mother particle tracing, then this can be done only in the 
# detector-level: the aim of this ttbb study (AN2019_158_v7)
# the match between the reco jet and the particle-level jets is less than 0.4

def init_lstm_(lstm: nn.LSTM):
    for name, param in lstm.named_parameters():
        if "weight_ih" in name:
            nn.init.xavier_uniform_(param)
        elif "weight_hh" in name:
            nn.init.orthogonal_(param)
        elif "bias" in name:
            nn.init.constant_(param, 0.0)

class JetAggregate(nn.Module):
    # by default, input features include jet four-momentum and the b-tagging score
    def __init__(self, d_in_feature=5, d_hiddens=[64,32,32], d_out = 20, dropout=0.0):
        super().__init__()
        self.inDims = [d_in_feature]+d_hiddens[:-1]
        self.outDims = d_hiddens
        self.dropout = nn.AlphaDropout(p=dropout)
        self.convLayers = nn.ModuleList([
            nn.Conv1d(d_in, d_out, kernel_size=1, stride=1, bias=True) for (d_in, d_out) in zip(self.inDims, self.outDims)
        ])
        # lecun normalization
        for d_in, layer in zip(self.inDims, self.convLayers):
            lecun_std = 1.0 / (d_in ** 0.5)
            nn.init.normal_(layer.weight, mean=0.0, std=lecun_std)
            nn.init.constant_(layer.bias, 0.0)
        self.mixJet = nn.LSTM(d_hiddens[-1], d_out, batch_first=True)
        init_lstm_(self.mixJet)

    def forward(self, x):       # x: [B, N, d]
        x = x.transpose(1, 2)   # x: [B, d, N]
        for embL in self.convLayers: 
            x = embL(x)         # x: [B, F, N]
            x = torch.selu(x)
            x = self.dropout(x)
        x = x.transpose(1, 2)   # x: [B, N, F]
        x = torch.flip(x, dims=[1])
        x_out, (hn, cn) = self.mixJet(x) 
        h_out = self.dropout(hn[-1])    # h_out: [B, D]
        return h_out


class mixingDNN(nn.Module):
    def __init__(self, d_globalVar=4, d_global_hiddens=[], d_jet_features=5, d_LSTM_hiddens=[64,32,32], d_jet_out = 20, dropout=0.0):
        super().__init__()
        # premixing layers
        self.dropout = nn.Dropout(p=dropout)
        if len(d_global_hiddens)>0:
            self.inDims = [d_globalVar]+d_global_hiddens[:-1]
            self.outDims = d_global_hiddens
            self.premixed_dim = d_global_hiddens[-1]
        else:
            self.inDims = []
            self.outDims= []
            self.premixed_dim = d_globalVar
        self.premixing = nn.ModuleList([
            nn.Linear(d_in, d_out) for (d_in, d_out) in zip(self.inDims, self.outDims)
        ])
        self.jet_processing = JetAggregate(d_in_feature=d_jet_features, d_hiddens=d_LSTM_hiddens, d_out = d_jet_out, dropout=dropout)
        # 4 layers add the 
        self.fusion_head = nn.Linear(self.premixed_dim + d_jet_out, 6)
    
    def forward(self, jet_var, global_var): # jet_var = [B, N, d], global_var = [B, D]
        jetMixing = self.jet_processing(jet_var)
        for globLayer in self.premixing:
            global_var = globLayer(global_var)
            global_var = torch.relu(global_var)
            global_var = self.dropout(global_var)
        # ordinary RNN+DNN
        x_out = torch.concat((jetMixing, global_var), 1)
        x_out = self.fusion_head(x_out)
        return x_out


##################################################
# to load dataset for attention like-structure
def loadDS(filePath, onlyMatch=0, onlySS=0, randomizeEntry=[]):
    dataset = np.load(filePath)
    if onlyMatch:
        dataset = dataset[dataset[:,-3]>=0]
    if onlySS:
        dataset = dataset[dataset[:,-1]==1]

        #np.save(shufflePath, permutation)

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

    # 0-27 4 jets with (pt, eta, b-tagscore, invmass_mu1, dr_mu1, invmass_mu2, dr_mu2)
    jetPart = dataset[:, :28]
    x_node = np.reshape(jetPart, [jetPart.shape[0], 4, 7])

    # 28 29 multiplicity, 30 31 total pt, 32-37 deltaR, 38-43 invmass
    x_global = dataset[:, 28:-3] 

    # last three entries label
    # [bjet_category, n_Muon, pass_same_sign]
    all_labels = dataset[:,-3:]
    all_labels = all_labels.astype(int)

    return (x_node, x_global, all_labels)


def train():
    model = mixingDNN(
        d_globalVar=16, 
        d_global_hiddens=[50,50,50],
        d_jet_features=7, 
        d_LSTM_hiddens=[128, 64, 64, 32, 32], 
        d_jet_out = 20, 
        dropout=0.1
    )
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = model.to(device)

    count_parameters(model)

    # new 2024 training
    train_2024_Path = "../ApplyTo2024/RNN_model/TTHH_DL_2B2W_batch1_RNN_features_train.npy"
    val_2024_Path = "../ApplyTo2024/RNN_model/TTHH_DL_2B2W_batch1_RNN_features_val.npy"

    x_node, x_global, y_train = loadDS(train_2024_Path, onlyMatch=1, onlySS=0)
    x_node_val, x_global_val, y_val = loadDS(val_2024_Path, onlyMatch=1, onlySS=0)

    # >>> x_node.shape
    # (62877, 4, 7)
    # >>> x_global.shape
    # (62877, 19)
    # >>> y_train.shape
    # (62877, 3)

    # reweighting 
    typeWeight = balanceTrainType(y_train, doTrain=1)
    class_weights_tensor = torch.tensor(typeWeight, dtype=torch.float32).to(device)

    # preparing the dataset for the usage
    trainDS = TensorDataset(
        torch.from_numpy(x_node).float(), 
        torch.from_numpy(x_global).float(), 
        torch.from_numpy(y_train[:, 0]).long())

    valDS = TensorDataset(
        torch.from_numpy(x_node_val).float(), 
        torch.from_numpy(x_global_val).float(), 
        torch.from_numpy(y_val[:, 0]).long())

    # define training options
    num_epochs = 200
    batch_size = 1024
    train_losses = []
    val_losses = []
    criterion = nn.CrossEntropyLoss(weight=class_weights_tensor)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)

    train_loader = DataLoader(
        trainDS,
        batch_size=batch_size,
        shuffle=True,
        drop_last=False
    )
    val_loader = DataLoader(
        valDS,
        batch_size=batch_size,
        shuffle=False,
        drop_last=False
    )
    for epoch in range(num_epochs):
        model.train()
        total_loss = 0
        for batch_node, batch_glob, batch_y  in train_loader:
            batch_node = batch_node.to(device)
            batch_glob = batch_glob.to(device)
            batch_y = batch_y.to(device)

            optimizer.zero_grad()
            out = model(batch_node, batch_glob)
            loss = criterion(out, batch_y)
            loss.backward()
            optimizer.step()
            total_loss += loss.item()

        train_losses.append(total_loss / len(train_loader))

        model.eval()
        val_loss = 0
        correct = 0
        with torch.no_grad():
            for val_batch_node, val_batch_glob, val_batch_y in val_loader:
                val_batch_node = val_batch_node.to(device)
                val_batch_glob = val_batch_glob.to(device)
                val_batch_y = val_batch_y.to(device)
                out = model(val_batch_node, val_batch_glob)
                loss = criterion(out, val_batch_y)
                val_loss += loss.item()

                pred = out.argmax(dim=1)
                correct += (pred == val_batch_y).sum().item()

        val_losses.append(val_loss / len(val_loader))
        val_acc = correct / len(valDS)

        print(f"Epoch {epoch+1}/{num_epochs} | Train Loss: {train_losses[-1]:.4f} | "
            f"Val Loss: {val_losses[-1]:.4f} | Val Acc: {val_acc:.4f}")
        
    torch.save(model.state_dict(), "ttbb_model_2024_dimuon_RNN_BJA.pt")

    plt.plot(train_losses, label='Train Loss')
    plt.plot(val_losses, label='Val Loss')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.legend()
    plt.grid()
    plt.title('Training and Validation Loss')
    plt.savefig("ttbb_model_2024_dimuon_RNN_BJA_train_loss.png")


def test():
    model = mixingDNN(
        d_globalVar=16, 
        d_global_hiddens=[50,50,50],
        d_jet_features=7, 
        d_LSTM_hiddens=[128, 64, 64, 32, 32], 
        d_jet_out = 20, 
        dropout=0.1
    )
    model_folder = "./"
    # model_folder = "noRNN" 
    trainedPara = torch.load(model_folder+ "/ttbb_model_2024_dimuon_RNN_BJA.pt")
    model.load_state_dict(trainedPara)

    test_2024_Path = "../ApplyTo2024/RNN_model/TTHH_DL_2B2W_batch2_RNN_features.npy"

    x_node_test, x_global_test, y_test = loadDS(test_2024_Path, onlyMatch=1, onlySS=1)
    nEvents = x_node_test.shape[0]

    testDS = TensorDataset(
        torch.from_numpy(x_node_test).float(), 
        torch.from_numpy(x_global_test).float(), 
        torch.from_numpy(y_test[:, 0]).long())

    batch_size=256
    test_loader = DataLoader(testDS, batch_size=batch_size, shuffle=False)
    prediction = []
    truthLabel = []
    for batch_node, batch_glob, batch_y  in test_loader:
        out = model(batch_node, batch_glob)
        pred = out.argmax(dim=1)
        prediction.extend(pred.cpu().numpy().tolist())
        truthLabel.extend(batch_y.cpu().numpy().tolist())
    prediction = np.array(prediction, dtype=int)
    truthLabel = np.array(truthLabel, dtype=int)

    metricPerformance = np.sum(prediction==truthLabel) / nEvents
    print("acc:", metricPerformance)
    
    # evaluate succeed rate
    plot_confusion_matrix_fractional(truthLabel, prediction, "_RNN_dimuon", ["unmatch", "b1b2", "b1b3", "b1b4", "b2b3", "b2b4", "b3b4"], ["b1b2", "b1b3", "b1b4", "b2b3", "b2b4", "b3b4"], 1)

if __name__ == "__main__":
    # train()
    test()
