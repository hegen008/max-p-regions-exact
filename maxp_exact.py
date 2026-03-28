from pulp import *
import highspy
import numpy as np
from libpysal import weights

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
        self.weight_factor = 1 + np.floor(np.log10(np.sum(np.triu(sim_mat, k=1))))

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
        self.x = LpVariable.dicts("x", (self.I_set, self.K_set, self.C_set), cat="Binary")
        self.t = LpVariable.dicts("t", (self.I_set, self.I_set), cat="Binary")

        # Objective function
        self.model += lpSum(self.x[i][j][0] for i in self.I_set for j in self.K_set), "Objective"

        # Define base model constraints
        for k in self.K_set:
            self.model += lpSum(self.x[i][k][0] for i in self.I_set) <= 1, f"Single_Root_{k}"
            self.model += lpSum(self.x[i][k][c] * self.spatial_attr[i] for i in self.I_set for c in self.C_set) >= self.threshold * lpSum(self.x[i][k][0] for i in self.I_set), f"Threshold_{k}"
        for i in self.I_set:
            self.model += lpSum(self.x[i][k][c] for k in self.K_set for c in self.C_set) == 1, f"Single_Assignment_{i}"
            for k in self.K_set:
                for c in self.C_set:
                    if c > 0:
                        self.model += self.x[i][k][c] <= lpSum(self.x[j][k][c-1] for j in self.spatial_weights.neighbors[i]), f"Adjacency_{i}_{k}_{c}"
                for j in self.I_set:
                    if j < i:
                        self.t[i][j] <= lpSum(self.x[i][k][c] - self.x[j][k][c] for c in self.C_set) + 1, f"x_t_Matching_{i}_{j}_{k}"

    # solve MIP model
    def solve(self, time):
        self.solver = getSolver("PULP_CBC_CMD", timeLimit=time)
        self.model.solve(self.solver)
        self.maxp = int(value(lpSum(self.x[i][k][0] for i in self.I_set for k in self.K_set)))
        assigned = {(i,k) for i in self.I_set for k in self.K_set for c in self.C_set if value(self.x[i][k][c]) > 0.9}
        for i,k in assigned:
            self.regions[i] = k

    