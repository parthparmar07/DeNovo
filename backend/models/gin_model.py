#!/usr/bin/env python3
"""
GIN Model Loader for Pre-trained PyTorch Geometric Models
==========================================================
Loads pre-trained GIN models that were trained with PyTorch Geometric
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import global_mean_pool, global_add_pool
from torch_geometric.data import Data
import numpy as np

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("⚠️ RDKit not available")


class GINConv(nn.Module):
    """
    GIN Convolution Layer matching PyTorch Geometric's structure
    """
    def __init__(self, input_dim, output_dim, edge_dim1=5, edge_dim2=3):
        super(GINConv, self).__init__()
        
        # MLP for node features (2 layers with hidden_dim = 2*output_dim)
        hidden_dim = 2 * output_dim
        self.mlp = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, output_dim)
        )
        
        # Two edge embeddings with different dimensions
        # edge_embedding1: 5 types, edge_embedding2: 3 types
        self.edge_embedding1 = nn.Embedding(edge_dim1, output_dim)
        self.edge_embedding2 = nn.Embedding(edge_dim2, output_dim)
    
    def forward(self, x, edge_index, edge_attr=None):
        """
        Forward pass with neighbor aggregation
        """
        # Manual message passing for GIN
        row, col = edge_index
        
        # Aggregate messages from neighbors
        out = torch.zeros_like(x)
        for i in range(x.size(0)):
            neighbors = col[row == i]
            if len(neighbors) > 0:
                neighbor_features = x[neighbors]
                if edge_attr is not None and self.edge_embedding1 is not None:
                    edge_features = self.edge_embedding1(edge_attr[row == i])
                    neighbor_features = neighbor_features + edge_features
                out[i] = neighbor_features.sum(dim=0)
        
        # Add self-loops (GIN's key feature: (1 + epsilon) * x)
        out = out + x
        
        # Apply MLP
        out = self.mlp(out)
        
        return out


class GINModel(nn.Module):
    """
    GIN Model matching the pre-trained checkpoint architecture
    Architecture from checkpoints:
    - x_embedding1, x_embedding2: Node feature embeddings
    - gnns: 5 GIN layers with edge embeddings
    - batch_norms: 5 batch norm layers
    - feat_lin: Feature projection layer
    - pred_head: Prediction head (3 layers)
    """
    
    def __init__(self, num_tasks=12, num_node_features=119, emb_dim=300, num_layers=5):
        super(GINModel, self).__init__()
        
        self.num_layers = num_layers
        self.emb_dim = emb_dim
        self.num_tasks = num_tasks
        
        # Node feature embeddings (matching x_embedding1, x_embedding2)
        # x_embedding2 has 3 features in the checkpoint
        self.x_embedding1 = nn.Embedding(num_node_features, emb_dim)
        self.x_embedding2 = nn.Embedding(3, emb_dim)
        
        # GIN layers (matching gnns.0 to gnns.4)
        # edge_embedding1: 5 types, edge_embedding2: 3 types
        self.gnns = nn.ModuleList()
        for layer in range(num_layers):
            self.gnns.append(GINConv(emb_dim, emb_dim, edge_dim1=5, edge_dim2=3))
        
        # Batch normalization layers (matching batch_norms.0 to batch_norms.4)
        self.batch_norms = nn.ModuleList()
        for layer in range(num_layers):
            self.batch_norms.append(nn.BatchNorm1d(emb_dim))
        
        # Feature linear projection (matching feat_lin) - projects to 512
        self.feat_lin = nn.Linear(emb_dim, 512)
        
        # Prediction head (matching pred_head.0, pred_head.2, pred_head.4)
        # 512 → 256 → 256 → num_tasks
        self.pred_head = nn.Sequential(
            nn.Linear(512, 256),
            nn.ReLU(),
            nn.Linear(256, 256),
            nn.ReLU(),
            nn.Linear(256, num_tasks)
        )
    
    def forward(self, data):
        """
        Forward pass through GIN model
        
        Args:
            data: PyTorch Geometric Data object with x, edge_index, edge_attr, batch
        
        Returns:
            predictions: Tensor of shape (batch_size, num_tasks)
        """
        x, edge_index, edge_attr, batch = data.x, data.edge_index, data.edge_attr, data.batch
        
        # Embed node features (assuming x contains integer indices)
        if x.dtype == torch.long:
            h = self.x_embedding1(x[:, 0]) + self.x_embedding2(x[:, 1] if x.size(1) > 1 else x[:, 0])
        else:
            # If x is continuous, use first embedding layer's weight as projection
            h = torch.matmul(x, self.x_embedding1.weight.T)
        
        # Apply GIN layers
        for i in range(self.num_layers):
            h = self.gnns[i](h, edge_index, edge_attr)
            h = self.batch_norms[i](h)
            h = F.relu(h)
        
        # Global pooling (sum over all nodes in each graph)
        h = global_add_pool(h, batch)
        
        # Feature projection
        h = self.feat_lin(h)
        h = F.relu(h)
        
        # Prediction head
        out = self.pred_head(h)
        
        return out


def smiles_to_graph_pyg(smiles):
    """
    Convert SMILES to PyTorch Geometric Data object
    
    Args:
        smiles: SMILES string
    
    Returns:
        Data object with x, edge_index, edge_attr
    """
    if not RDKIT_AVAILABLE:
        raise ImportError("RDKit is required for SMILES processing")
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    
    # Add hydrogens
    mol = Chem.AddHs(mol)
    
    # Node features (atoms)
    num_atoms = mol.GetNumAtoms()
    x = []
    
    for atom in mol.GetAtoms():
        # Simple atomic number encoding
        atomic_num = atom.GetAtomicNum()
        x.append([min(atomic_num, 118), 0])  # Use atomic number as feature
    
    x = torch.tensor(x, dtype=torch.long)
    
    # Edge features (bonds)
    edge_index = []
    edge_attr = []
    
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        
        # Bond type encoding
        bond_type = bond.GetBondType()
        if bond_type == Chem.rdchem.BondType.SINGLE:
            bond_feature = 0
        elif bond_type == Chem.rdchem.BondType.DOUBLE:
            bond_feature = 1
        elif bond_type == Chem.rdchem.BondType.TRIPLE:
            bond_feature = 2
        else:
            bond_feature = 0
        
        # Add both directions (undirected graph)
        edge_index.extend([[i, j], [j, i]])
        edge_attr.extend([bond_feature, bond_feature])
    
    if len(edge_index) == 0:
        # Single atom molecule - add self-loop
        edge_index = [[0, 0]]
        edge_attr = [0]
    
    edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
    edge_attr = torch.tensor(edge_attr, dtype=torch.long)
    
    # Create Data object
    data = Data(x=x, edge_index=edge_index, edge_attr=edge_attr)
    
    return data


def load_gin_model(checkpoint_path, num_tasks, device='cpu'):
    """
    Load pre-trained GIN model from checkpoint
    
    Args:
        checkpoint_path: Path to .pth file
        num_tasks: Number of prediction tasks
        device: 'cpu' or 'cuda'
    
    Returns:
        Loaded GIN model
    """
    # Create model
    model = GINModel(num_tasks=num_tasks)
    
    # Load checkpoint
    checkpoint = torch.load(checkpoint_path, map_location=device)
    
    # Load state dict
    model.load_state_dict(checkpoint, strict=True)
    
    model.to(device)
    model.eval()
    
    return model


if __name__ == "__main__":
    # Test model loading
    print("Testing GIN Model Loader...")
    
    # Test SMILES to graph conversion
    if RDKIT_AVAILABLE:
        smiles = "CCO"  # Ethanol
        data = smiles_to_graph_pyg(smiles)
        print(f"✅ Converted SMILES to graph: {data}")
        print(f"   Nodes: {data.x.shape}, Edges: {data.edge_index.shape}")
    else:
        print("⚠️ RDKit not available - skipping SMILES test")
    
    # Test model creation
    model = GINModel(num_tasks=12)
    print(f"✅ Created GIN model with {sum(p.numel() for p in model.parameters())} parameters")
