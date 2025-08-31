import torch
import torch.nn as nn
from torch_geometric.nn import NNConv
import torch.nn.functional as F

from math import comb

##################################################
# structures needed by model
##################################################

# attention pooling for [batch_size, token_size, dimension] structure output
class AttnPool(nn.Module):
    def __init__(self, d_model, nhead=4, dropout=0.0):
        super().__init__()
        # randomize a virtual query for the pooler, shared by the whole batch and has one word
        self.query = nn.Parameter(torch.randn(1, 1, d_model)) 
        self.mha = nn.MultiheadAttention(embed_dim=d_model, num_heads=nhead, dropout=dropout, batch_first=True)

    def forward(self, x):  # x: [B,S,D]
        B = x.size(0)
        q = self.query.expand(B, -1, -1)         # [B,1,D]
        y, _ = self.mha(q, x, x)                 # [B,1,D]
        return y.squeeze(1)                      # [B,D]
    
# Positional encoding after the initial encoding of nodes
class PositionalEncoding(nn.Module):
    def __init__(self, d_model, max_len=5000):
        super().__init__()
        pe = torch.zeros(max_len, d_model)
        position = torch.arange(0, max_len, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, d_model, 2).float() * (-torch.log(torch.tensor(10000.0)) / d_model))
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        pe = pe.unsqueeze(0)  # shape: (1, max_len, d_model)
        self.register_buffer('pe', pe)

    def forward(self, x):
        x = x + self.pe[:, :x.size(1)]
        return x

# MLP head to merge global variables
class LateFusionHead(nn.Module):
    def __init__(self, d_model, global_dim, n_classes, hidden=128, dropout=0.1):
        super().__init__()
        assert global_dim > 0, "global_dim must be > 0 for this head"
        self.proj_global = nn.Sequential(
            nn.LayerNorm(global_dim),
            nn.Linear(global_dim, d_model), nn.GELU(), nn.Dropout(dropout)
        )
        fuse_in = d_model + d_model           # h_pool(D) + g_emb(D)
        self.cls = nn.Sequential(
            nn.LayerNorm(fuse_in),
            nn.Linear(fuse_in, hidden), nn.GELU(), nn.Dropout(dropout),
            nn.Linear(hidden, n_classes)
        )
    def forward(self, h_pool, x_global):      # h_pool:[B,D], x_global:[B,G]
        g_emb = self.proj_global(x_global)    # [B,D]
        z = torch.cat([h_pool, g_emb], dim=-1)
        return self.cls(z)                    # [B,C]

# assist classes for mixing the layers output
class LayerMix(nn.Module):
    def __init__(self, n_layers: int):
        super().__init__()
        self.w = nn.Parameter(torch.zeros(n_layers))
    def forward(self, layer_list):                   # list of [B,S,D]
        a = torch.softmax(self.w, dim=0)             # [L]
        out = sum(w * x for w, x in zip(a, layer_list))
        return out                                   # [B,S,D]

def copy_layer(layer: nn.Module) -> nn.Module:
    import copy
    return copy.deepcopy(layer)

# transformer encoder with layer mixing
class CollectingTransformerEncoder(nn.Module):
    def __init__(self, encoder_layer: nn.TransformerEncoderLayer, num_layers: int, norm: nn.Module | None = None):
        super().__init__()
        self.layers = nn.ModuleList([copy_layer(encoder_layer) for _ in range(num_layers)])
        self.norm = norm

    def forward(self, src, src_mask=None, src_key_padding_mask=None):
        out_list = []
        x = src
        for layer in self.layers:
            x = layer(x, src_mask=src_mask, src_key_padding_mask=src_key_padding_mask)
            out_list.append(x)
        if self.norm is not None:
            x = self.norm(x)
            out_list[-1] = x 
        return out_list

# basic GNN model, provide batch encoding processing
class FullyConnectGNN(nn.Module):
    def __init__(self, num_nodes=4, input_dim=5, edge_dim=1, hidden_dim=16, dropout=0.3, initialEncoding=False):
        super().__init__()
        self.hidden_dim = hidden_dim
        self.num_nodes  = num_nodes
        self.initialEncoding = initialEncoding

        # fully connect GNN edge preparation, registered as buffer to the model, do not write it again!
        self.num_edges_per_graph = comb(num_nodes, 2) * 2
        edge_index_base = torch.cat([(ij:=torch.triu_indices(num_nodes, num_nodes, 1, dtype=torch.long)), ij.flip(0)], dim=1)
        self.register_buffer('edge_index_base', edge_index_base, persistent=False)

        # node encoding (only need for the first GNN)
        self.nodeEncode = nn.Sequential(
            nn.LayerNorm(input_dim),
            nn.Linear(input_dim, hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout)
        )

        # edge encoding
        self.edgeEncode = nn.Sequential(
            nn.Linear(edge_dim, 64),
            nn.GELU(),
            nn.Linear(64, hidden_dim * hidden_dim)
        )

        # NNConv
        self.conv = NNConv(
            in_channels=hidden_dim,
            out_channels=hidden_dim,
            nn=self.edgeEncode,
            aggr='mean'
        )

        # FFN
        self.ffn = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim * 2),
            nn.GELU(),
            nn.Linear(hidden_dim * 2, hidden_dim),
            nn.Dropout(dropout)
        )

        self.ln_out = nn.LayerNorm(hidden_dim)

    # In GNN, nodes and edges from different events inside one batch are represented by offset
    @torch.no_grad()
    def _make_batched_edge_index(self, batch_size: int, device: torch.device):
        N  = self.num_nodes
        E  = self.num_edges_per_graph
        ei = self.edge_index_base  # [2,E] on correct device because registered buffer
        # offsets: [B,1] * [2,E] → [B,2,E] → reshape to [2,B*E]
        offsets = (torch.arange(batch_size, device=device) * N).view(-1,1,1)  # [B,1,1], [[[0]], [[4]], ...]
        ei_b = ei.unsqueeze(0) + offsets                                      # [B,2,E], [ei, ei+4, ...]
        ei_b = ei_b.permute(1,0,2).reshape(2, batch_size*E)                   # [2, B*E], all edges index 0, 1
        return ei_b

    def forward(self, x: torch.Tensor, x_edge: torch.Tensor):
        n_batches, num_nodes, Fin = x.shape
        # safety checking
        assert num_nodes == self.num_nodes, f"expected {self.num_nodes} nodes, got {num_nodes}"
        assert x_edge.size(1) == self.num_edges_per_graph, \
            f"x_edge second dim must be {self.num_edges_per_graph} (directed edges)"
        dev = x.device
       
        # initial encoding
        if self.initialEncoding:
            x = self.nodeEncode(x)                                              # [B, n_node, H]

        # for batching graph computation
        x_flat      = x.reshape(n_batches*num_nodes, self.hidden_dim)           # [B*n_node, H]
        edge_attr   = x_edge.reshape(n_batches*self.num_edges_per_graph, -1)    # [B*n_edge, edge_dim]
        edge_index  = self._make_batched_edge_index(n_batches, dev)             # [2, B*n_edge]

        x_msg = self.conv(x_flat, edge_index, edge_attr)                        # [B*4, H]
        x_msg = F.gelu(x_msg)
        x_out = x_flat + x_msg
        x_out = x_out + self.ffn(x_out)
        x_out = self.ln_out(x_out)

        return x_out.view(n_batches, num_nodes, self.hidden_dim)

##################################################
# utilize models
##################################################
# pure GNN
class Bjet_assignment_GNN(nn.Module):
    def __init__(self, num_nodes=4, input_dim=5, edge_dim=1, hidden_dim=16, global_dim=4, num_classes=6, dropout=0.3):
        super().__init__()
        self.hidden_dim=hidden_dim
        self.global_dim=global_dim
        self.num_classes=num_classes

        # first GNN need initial encoding for rising dimension, later not
        self.GNN1 = FullyConnectGNN(num_nodes, input_dim, edge_dim, hidden_dim, dropout, 1)
        self.GNN2 = FullyConnectGNN(num_nodes, hidden_dim, edge_dim, hidden_dim, dropout, 0)

        # then go to final pooling from DNN
        self.pool = AttnPool(hidden_dim, nhead=8, dropout=dropout)
        # concatenate global to categorization MLP and output
        self.head = LateFusionHead(hidden_dim, global_dim=global_dim, n_classes=num_classes, hidden=128, dropout=dropout)

    def forward(self, x_node, x_edge, x_global):
        x_node = self.GNN1(x_node, x_edge)
        x_node = self.GNN2(x_node, x_edge)
        x_node = self.pool(x_node)
        out = self.head(x_node, x_global)
        return out

##################################################
# Pure Transformer
class Bjet_assignment_Transformer(nn.Module):
    def __init__(self, feature_dim, d_model, nhead, num_layers,
                 n_classes, global_dim, enc_dropout=0.3, pool_dropout=0.3, use_norm_first=False):
        super().__init__()
        # simple initial encoder
        self.input_proj = nn.Linear(feature_dim, d_model)
        self.layer_norm = nn.LayerNorm(d_model)
        self.pos_encoder = PositionalEncoding(d_model)

        encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model, nhead=nhead, dim_feedforward=d_model*4,
            activation = "gelu", dropout=enc_dropout, 
            batch_first=True, norm_first=use_norm_first
        )

        # for applying LayerMix
        self.encoder = CollectingTransformerEncoder(encoder_layer, num_layers=num_layers, norm=None)
        self.mix = LayerMix(num_layers)
        self.pool = AttnPool(d_model, nhead=nhead, dropout=pool_dropout)

        # global variables fusion
        self.head = LateFusionHead(d_model, global_dim, n_classes, hidden=128, dropout=0.1)

    def forward(self, x_node, x_global): 
        x = self.input_proj(x_node)
        x = self.pos_encoder(x)
        x = self.layer_norm(x)

        layers_out = self.encoder(x)
        out = self.mix(layers_out)
        out = self.pool(out)

        out = self.head(out, x_global)
        return out

##################################################
# utilize model: GNN+Transformer
class Bjet_assignment_GNN_Transformer(nn.Module):
    def __init__(self, feature_dim, edge_dim, d_model, nhead, num_layers,
                 n_classes, global_dim, enc_dropout=0.3, pool_dropout=0.3, use_norm_first=False):
        super().__init__()
        # encoding the bjet information together with di-bjet variables
        self.input_proj = FullyConnectGNN(input_dim=feature_dim, edge_dim=edge_dim, hidden_dim=d_model, dropout=enc_dropout, initialEncoding=True)
        self.pos_encoder = PositionalEncoding(d_model)
        self.layer_norm = nn.LayerNorm(d_model)

        encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model, nhead=nhead, dim_feedforward=d_model*4,
            activation = "gelu", dropout=enc_dropout, 
            batch_first=True, norm_first=use_norm_first
        )

        # for applying LayerMix
        self.encoder = CollectingTransformerEncoder(encoder_layer, num_layers=num_layers, norm=None)
        self.mix = LayerMix(num_layers)
        self.pool = AttnPool(d_model, nhead=nhead, dropout=pool_dropout)

        # global variables fusion
        self.head = LateFusionHead(d_model, global_dim, n_classes, hidden=128, dropout=0.1)

    def forward(self, x_node, x_edge, x_global): 
        x = self.input_proj(x_node, x_edge)
        x = self.pos_encoder(x)
        x = self.layer_norm(x)

        layers_out = self.encoder(x)
        out = self.mix(layers_out)
        out = self.pool(out)

        out = self.head(out, x_global)
        return out