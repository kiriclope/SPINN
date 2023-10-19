import os
import numpy as np
import torch
from torch import nn
from yaml import safe_load

torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False

class Connections(nn.Module):
    def __init__(self, conf_file, sim_name, repo_root, **kwargs):
        super().__init__()

        # load parameters
        self.loadConfig(conf_file, sim_name, repo_root, **kwargs)

        # set csts
        self.initConst()
        
        # Create recurrent layer with nn.Linear
        self.Wab = [[None]*self.N_POP for _ in range(self.N_POP)]
        
        for i_pop in range(self.N_POP):
            for j_pop in range(self.N_POP):
                Wij = nn.Linear(self.Na[i_pop], self.Na[j_pop], bias=(i_pop==j_pop), dtype=self.FLOAT, device=self.device)
                Wij.weight.data = self.initWeights(i_pop, j_pop)
                self.Wab[i_pop][j_pop] = Wij
                self.Wab[i_pop][j_pop] = self.Wab[i_pop][j_pop]

        block_matrix = self.form_block_matrix(self.Wab)
        sparse_matrix = block_matrix.to_sparse()
        
    def form_block_matrix(self, layers_2d):
        row_blocks = [torch.cat([layer.weight.view(layer.out_features, layer.in_features) for layer in row], dim=1) for row in layers_2d]
        blocks = torch.cat(row_blocks, dim=0)
        return blocks
                                
    def loadConfig(self, conf_file, sim_name, repo_root, **kwargs):
        # Loading configuration file
        conf_path = repo_root + '/conf/'+ conf_file
        print('Loading config from', conf_path)
        param = safe_load(open(conf_path, "r"))
        
        param["FILE_NAME"] = sim_name
        param.update(kwargs)

        for k, v in param.items():
            setattr(self, k, v)

        self.DATA_PATH = repo_root + "/data/simul/"
        self.MAT_PATH = repo_root + "/data/matrix/"

        if not os.path.exists(self.DATA_PATH):
            os.makedirs(self.DATA_PATH)

        if not os.path.exists(self.MAT_PATH):
            os.makedirs(self.MAT_PATH)

        if self.FLOAT_PRECISION == 32:
            self.FLOAT = torch.float
        else:
            self.FLOAT = torch.float64

        self.device = torch.device(self.DEVICE)
    
    def initWeights(self, i_pop, j_pop):

        Na = self.Na[i_pop]
        Nb = self.Na[j_pop]
        Kb = self.Ka[j_pop]

        Pij = 1.0

        if 'cos' in self.STRUCTURE[i_pop][j_pop]:
            theta = torch.arange(0, 2.0 * torch.pi, 2.0 * torch.pi / Na, dtype=self.FLOAT, device=self.device)
            phi = torch.arange(0, 2.0 * torch.pi, 2.0 * torch.pi / Nb, dtype=self.FLOAT, device=self.device)

            i, j = torch.meshgrid(torch.arange(Na, device=self.device), torch.arange(Nb, device=self.device), indexing="ij")
            
            theta_diff = theta[i] - phi[j]

            if 'spec' in self.STRUCTURE[i_pop][j_pop]:
                self.KAPPA[i_pop][j_pop] = self.KAPPA[i_pop][j_pop] / np.sqrt(self.Ka[j_pop])

            Pij = 1.0 + 2.0 * self.KAPPA[i_pop][j_pop] * torch.cos(theta_diff - self.PHASE)

        if 'sparse' in self.CONNECTIVITY:
            if self.VERBOSE:
                print('Sparse random connectivity ')
            Cij = self.Jab[i_pop][j_pop] * (torch.rand(Na, Nb, device=self.device) < Kb / Nb * Pij)

        if 'all2all' in self.CONNECTIVITY:
            if self.VERBOSE:
                print('All to all connectivity ')
            Cij = self.Jab[i_pop][j_pop] * Pij / Nb

        if self.VERBOSE:
            if "cos" in self.STRUCTURE[i_pop][j_pop]:
                if "spec" in self.STRUCTURE[i_pop][j_pop]:
                    print('with weak cosine structure')
                else:
                    print('with strong cosine structure')

        return Cij

    def initConst(self):
        self.Na = []
        self.Ka = []

        for i_pop in range(self.N_POP):
            self.Na.append(int(self.N_NEURON * self.frac[i_pop]))
            # self.Ka.append(self.K * const.frac[i_pop])
            self.Ka.append(self.K)
            
        self.Na = torch.tensor(self.Na, dtype=torch.int)
        self.Ka = torch.tensor(self.Ka, dtype=self.FLOAT)
        self.csumNa = torch.cat((torch.tensor([0]), torch.cumsum(self.Na, dim=0)))
        
        self.STRUCTURE = np.array(self.STRUCTURE).reshape(self.N_POP, self.N_POP)
        self.SIGMA = torch.tensor(self.SIGMA, dtype=self.FLOAT).view(self.N_POP, self.N_POP)
        self.KAPPA = torch.tensor(self.KAPPA, dtype=self.FLOAT).view(self.N_POP, self.N_POP)
        self.PHASE = self.PHASE * torch.pi / 180.0
