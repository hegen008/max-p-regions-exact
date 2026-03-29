from pulp import *
import highspy
import numpy as np
from libpysal import weights
from dataclasses import dataclass, field

@dataclass
class MaxPConfig:
    bound_num_regions: bool = False
    bound_contiguity: bool = False
    merge_leaves: bool = False
    independent_split: bool = False
    reduce_size_of_t: bool = False
    preassign_roots: bool = False
    exclude_roots: bool = False
    max_attr_for_root: bool = False
    min_adj_order: bool = False
    sort_regions: bool = False

class MaxPExact():
    # array initialization
    def __init__(self, adj_mat, sim_mat, spatial_attr, threshold):
        # Initialize data
        self.adj_mat = adj_mat
        self.sim_mat = sim_mat
        self.spatial_attr = spatial_attr
        self.threshold = threshold

        # Initialize derived information
        self.num_areas = len(spatial_attr)
        self.spatial_weights = weights.full2W(adj_mat)
        self.regions = np.empty(self.num_areas, dtype=object)
        self.maxp = None
        self.obj = None
        self.weight_factor = 10**(1 + np.floor(np.log10(np.sum(np.triu(sim_mat, k=1)))))

    # geopandas initialization
    @classmethod
    def from_gpd(self):
        pass

    # construct MIP model
    def construct(self):
        # Initialize model
        self.model = LpProblem("MaxP_Model", LpMaximize)

        # Define iterative ranges
        self.I_set = range(self.num_areas)
        self.K_set = range(self.num_areas)
        self.C_set = range(self.num_areas)

        # Decision variables
        self.x = LpVariable.dicts("var_x", (self.I_set, self.K_set, self.C_set), cat="Binary")
        self.t = LpVariable.dicts("var_t", [(i, j) for i in self.I_set for j in self.I_set if j > i], cat="Binary")

        # Objective function
        self.model += lpSum(self.x[i][j][0] for i in self.I_set for j in self.K_set) * self.weight_factor + lpSum(self.t[i,j] * self.sim_mat[i][j] for i in self.I_set for j in self.I_set if j > i), "Objective"

        # Define base model constraints
        self.model.extend([ # Single Root Constraints
            lpSum(self.x[i][k][0] for i in self.I_set) <= 1
            for k in self.K_set
        ]) 
        self.model.extend([ # Threshold Constraints
            lpSum(self.x[i][k][c] * self.spatial_attr[i] for i in self.I_set for c in self.C_set) >= self.threshold * lpSum(self.x[i][k][0] for i in self.I_set)
            for k in self.K_set
        ])
        self.model.extend([ # Single Assignment Constraints
            lpSum(self.x[i][k][c] for k in self.K_set for c in self.C_set) == 1
            for i in self.I_set
        ])
        self.model.extend([ # Adjacency Constraints
            self.x[i][k][c] <= lpSum(self.x[j][k][c-1] for j in self.spatial_weights.neighbors[i])
            for i in self.I_set 
            for k in self.K_set 
            for c in self.C_set if c > 0
        ])
        self.model.extend([ # x-t Matching Constraints
            self.t[i,j] <= lpSum(self.x[i][k][c] - self.x[j][k][c] for c in self.C_set) + 1
            for i in self.I_set 
            for j in self.I_set if j > i 
            for k in self.K_set
        ])

    # solve MIP model
    def solve(self, time):
        self.model.solve(HiGHS(timeLimit=time, keepFiles=False))
        self.maxp = int(value(lpSum(self.x[i][k][0] for i in self.I_set for k in self.K_set)))
        self.obj = value(self.model.objective)

        assigned = {(i,k) for i in self.I_set for k in self.K_set for c in self.C_set if value(self.x[i][k][c]) > 0.9}
        for i,k in assigned:
            self.regions[i] = k

    