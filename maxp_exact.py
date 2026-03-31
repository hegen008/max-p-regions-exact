from pulp import *
import highspy
import numpy as np
from libpysal import weights
from dataclasses import dataclass, field
from copy import deepcopy

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
    for next_ind in weights.neighbors[path[0]]:
        if next_ind not in excluded:
            depth = _recursive_step(weights, spatial_attr, threshold, [next_ind] + path, excluded | set(weights.neighbors[path[0]]))
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
    return np.array([id_map[s] for s in solution])

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
    sort_region_roots: bool = False

class MaxPExact():
    # array initialization
    def __init__(self, adj_mat, sim_mat, spatial_attr, threshold):
        # Initialize data - sort on intake
        _idx = np.argsort(spatial_attr)[::-1]
        self._orig_map = np.argsort(_idx)
        self._adj_mat = adj_mat[_idx, :][:, _idx]
        self._sim_mat = sim_mat[_idx, :][:, _idx]
        self._spatial_attr = spatial_attr[_idx]
        self.threshold = threshold

        # Initialize derived information
        self.num_areas = len(self._spatial_attr)
        self._spatial_weights = weights.full2W(self._adj_mat)
        self._regions = np.full(self.num_areas, None, dtype=object)
        self.maxp = None
        self.obj = None
        self.weight_factor = 10**(1 + np.floor(np.log10(np.sum(np.triu(self._sim_mat, k=1)))))

    @property
    def spatial_attr(self):
        return self._spatial_attr[self._orig_map]
    
    @property
    def regions(self):
        return self._regions[self._orig_map]
    
    @property
    def sim_mat(self):
        return self._sim_mat[self._orig_map, :][:, self._orig_map]
    
    @property
    def adj_mat(self):
        return self._adj_mat[self._orig_map, :][:, self._orig_map]

    # geopandas initialization
    @classmethod
    def from_gpd(self):
        pass

    # construct MIP model
    def construct(self, config):
        # Define iterative ranges
        self._I_set = range(self.num_areas)
        if config.bound_num_regions:
            self._K_set = range(_bound_num_regions(self._spatial_attr, self.threshold))
        else:
            self._K_set = range(self.num_areas)
        if config.bound_contiguity:
            self._C_set = range(_bound_contiguity(self._spatial_weights, self._spatial_attr, self.threshold))
        else:
            self._C_set = range(self.num_areas)

        # Initialize model
        self.model = LpProblem("MaxP_Model", LpMaximize)

        # Decision variables
        self.x = LpVariable.dicts("var_x", (self._I_set, self._K_set, self._C_set), cat="Binary")
        self.t = LpVariable.dicts("var_t", [(i, j) for i in self._I_set for j in self._I_set if j > i], cat="Binary")

        # Objective function
        self.model += lpSum(self.x[i][j][0] for i in self._I_set for j in self._K_set) * self.weight_factor + lpSum(self.t[i,j] * self._sim_mat[i][j] for i in self._I_set for j in self._I_set if j > i), "Objective"

        # Define base model constraints
        self.model.extend([ # Single Root Constraints
            lpSum(self.x[i][k][0] for i in self._I_set) <= 1
            for k in self._K_set
        ]) 
        self.model.extend([ # Threshold Constraints
            lpSum(self.x[i][k][c] * self._spatial_attr[i] for i in self._I_set for c in self._C_set) >= self.threshold * lpSum(self.x[i][k][0] for i in self._I_set)
            for k in self._K_set
        ])
        self.model.extend([ # Single Assignment Constraints
            lpSum(self.x[i][k][c] for k in self._K_set for c in self._C_set) == 1
            for i in self._I_set
        ])
        self.model.extend([ # Adjacency Constraints
            self.x[i][k][c] <= lpSum(self.x[j][k][c-1] for j in self._spatial_weights.neighbors[i])
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
            excluded_roots = _find_excluded_roots(self._spatial_weights, self._spatial_attr, self.threshold)
            self.model.extend([
                lpSum(self.x[i][k][0] for k in self._K_set) == 0
                for i in excluded_roots
            ])
        if config.preassign_roots: # Preassign Roots
            if max(self._spatial_attr) < self.threshold:
                temp_attr = np.array(self._spatial_attr, copy=True)[list(excluded_roots)] = -np.inf
                self.model += self.x[np.argmax(temp_attr)][0][0] == 1
            else:
                over_thresh = [i for i in self._I_set if self._spatial_attr[i] >= self.threshold]
                self.model.extend([
                    self.x[i][ind][0] == 1
                    for ind,i in enumerate(over_thresh)
                ])
        if config.max_attr_for_root:
            self.model.extend([
                lpSum(self._spatial_attr[j] * self.x[j][k][0] for j in self._I_set) >= self._spatial_attr[i] * self.x[i][k][c]
                for c in self._C_set if c > 0
                for i in self._I_set if i not in excluded_roots
                for k in self._K_set
            ])
        if config.min_adj_order:
            self.model.extend([
                self.num_areas * (1 - self.x[i][k][c]) >= lpSum(self.x[j][k][d] for j in self._spatial_weights.neighbors[i] for d in range(0,c-1))
                for c in self._C_set if c > 1
                for i in self._I_set
                for k in self._K_set
            ])
        if config.sort_region_roots:
            self.model.extend([
                lpSum(self._spatial_attr[i] * self.x[i][k-1][0] for i in self._I_set) >= lpSum(self._spatial_attr[i] * self.x[i][k][0] for i in self._I_set)
                for k in self._K_set if k > 0
            ])

    # solve MIP model
    def solve(self, time):
        self.model.solve(HiGHS(timeLimit=time, msg=True, keepFiles=False, options=['mip_abs_gap=1e-4', 'mip_rel_gap=1e-10']))
        self.maxp = int(value(lpSum(self.x[i][k][0] for i in self._I_set for k in self._K_set)))
        self.obj = value(self.model.objective)

        assigned = {(i,k) for i in self._I_set for k in self._K_set for c in self._C_set if value(self.x[i][k][c]) > 0.9}
        for i,k in assigned:
            self._regions[i] = k
        self._regions = _standardize_solution(self._regions)

    