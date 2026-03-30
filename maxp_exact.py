from pulp import *
import highspy
import numpy as np
from libpysal import weights
from dataclasses import dataclass, field

def _bound_num_regions(spatial_attr, threshold):
    over_thres = np.sum(spatial_attr > threshold)
    under_thres = np.sum(spatial_attr[spatial_attr < threshold])
    under_thres //= threshold
    return over_thres + under_thres

def _can_split(spatial_attr, threshold, path):
    attr_list = spatial_attr[path]
    from_head = np.cumsum(attr_list) >= threshold
    from_tail = np.cumsum(attr_list[::-1])[::-1] >= threshold
    split_points = from_head[:-1] & from_tail[1:]
    return np.any(split_points)

def _recursive_step(weights, spatial_attr, threshold, path, excluded):
    if _can_split(spatial_attr, threshold, path):
        return len(path) - 1
    if not (set(weights.neighbors[path[0]]) - excluded):
        return len(path)
    max_q = 0
    for next in weights.neighbors[path[0]]:
        if next not in excluded:
            depth = _recursive_step(weights, spatial_attr, threshold, [next] + path, excluded | set(weights.neighbors[path[0]]))
            if depth > max_q:
                max_q = depth
    return max_q

def _bound_contiguity(weights, spatial_attr, threshold):
    max_q = 0
    for i in range(len(spatial_attr)):
        depth = _recursive_step(weights, spatial_attr, threshold, [i], set())
        if depth > max_q:
            max_q = depth
    return max_q

def _standardize_solution(solution):
    id_map = {}
    counter = 0
    for item in solution:
        if item not in id_map: # new region ID
            id_map[item] = counter
            counter += 1
    return [id_map[s] for s in solution]

def _find_excluded_roots(weights, spatial_attr, threshold):
    excluded = set()
    for i in sorted(range(len(spatial_attr)), key=lambda i: spatial_attr[i]):
        if spatial_attr[i] >= threshold:
            break
        contig_excl = {i}
        while True:
            size = len(contig_excl)
            contig_excl |= {neigh for ex in contig_excl for neigh in weights.neighbors[ex] if neigh in excluded}
            if len(contig_excl) == size:
                break
        if sum(spatial_attr[ex] for ex in contig_excl) < threshold:
            excluded.add(i)
    return excluded

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
    def construct(self, config):
        # Initialize model
        self.model = LpProblem("MaxP_Model", LpMaximize)

        # Define iterative ranges
        self._I_set = range(self.num_areas)
        if config.bound_num_regions:
            self._K_set = range(_bound_num_regions(self.spatial_attr, self.threshold))
        else:
            self._K_set = range(self.num_areas)
        if config.bound_contiguity:
            self._C_set = range(_bound_contiguity(self.spatial_weights, self.spatial_attr, self.threshold))
        else:
            self._C_set = range(self.num_areas)

        # Decision variables
        self.x = LpVariable.dicts("var_x", (self._I_set, self._K_set, self._C_set), cat="Binary")
        self.t = LpVariable.dicts("var_t", [(i, j) for i in self._I_set for j in self._I_set if j > i], cat="Binary")

        # Objective function
        self.model += lpSum(self.x[i][j][0] for i in self._I_set for j in self._K_set) * self.weight_factor + lpSum(self.t[i,j] * self.sim_mat[i][j] for i in self._I_set for j in self._I_set if j > i), "Objective"

        # Define base model constraints
        self.model.extend([ # Single Root Constraints
            lpSum(self.x[i][k][0] for i in self._I_set) <= 1
            for k in self._K_set
        ]) 
        self.model.extend([ # Threshold Constraints
            lpSum(self.x[i][k][c] * self.spatial_attr[i] for i in self._I_set for c in self._C_set) >= self.threshold * lpSum(self.x[i][k][0] for i in self._I_set)
            for k in self._K_set
        ])
        self.model.extend([ # Single Assignment Constraints
            lpSum(self.x[i][k][c] for k in self._K_set for c in self._C_set) == 1
            for i in self._I_set
        ])
        self.model.extend([ # Adjacency Constraints
            self.x[i][k][c] <= lpSum(self.x[j][k][c-1] for j in self.spatial_weights.neighbors[i])
            for i in self._I_set 
            for k in self._K_set 
            for c in self._C_set if c > 0
        ])
        self.model.extend([ # x-t Matching Constraints
            self.t[i,j] <= lpSum(self.x[i][k][c] - self.x[j][k][c] for c in self._C_set) + 1
            for i in self._I_set 
            for j in self._I_set if j > i 
            for k in self._K_set
        ])
        excluded_roots = set()
        if config.exclude_roots: # Exclude Roots
            excluded_roots = _find_excluded_roots(self.spatial_weights, self.spatial_attr, self.threshold)
            self.model.extend([
                lpSum(self.x[i][k][0] for k in self._K_set) == 0
                for i in excluded_roots
            ])
        if config.preassign_roots: # Preassign Roots
            if max(self.spatial_attr) < self.threshold:
                temp_attr = np.array(self.spatial_attr, copy=True)[list(excluded_roots)] = -np.inf
                self.model += self.x[np.argmax(temp_attr)][0][0] == 1
            else:
                over_thresh = [i for i in self._I_set if self.spatial_attr[i] >= self.threshold]
                self.model.extend([
                    self.x[i][idx][0] == 1
                    for idx,i in enumerate(over_thresh)
                ])


    # solve MIP model
    def solve(self, time):
        self.model.solve(HiGHS(timeLimit=time, keepFiles=False))
        self.maxp = int(value(lpSum(self.x[i][k][0] for i in self._I_set for k in self._K_set)))
        self.obj = value(self.model.objective)

        assigned = {(i,k) for i in self._I_set for k in self._K_set for c in self._C_set if value(self.x[i][k][c]) > 0.9}
        for i,k in assigned:
            self.regions[i] = k
        self.regions = _standardize_solution(self.regions)

    