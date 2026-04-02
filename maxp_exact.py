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

def _recursive_step(spatial_weights, spatial_attr, threshold, path, excluded):
    if _can_split(spatial_attr, threshold, path):
        return len(path) - 1
    if not (set(spatial_weights.neighbors[path[0]]) - excluded):
        return len(path)
    max_q = 0
    for next_ind in spatial_weights.neighbors[path[0]]:
        if next_ind not in excluded:
            depth = _recursive_step(spatial_weights, spatial_attr, threshold, [next_ind] + path, excluded | set(spatial_weights.neighbors[path[0]]))
            if depth > max_q:
                max_q = depth
    return max_q

def _bound_contiguity(spatial_weights, spatial_attr, threshold):
    max_q = 0
    for i in range(len(spatial_attr)):
        depth = _recursive_step(spatial_weights, spatial_attr, threshold, [i], set())
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

def _find_excluded_roots(spatial_weights, spatial_attr, threshold):
    excluded = set()
    for i in sorted(range(len(spatial_attr)), key=lambda i: spatial_attr[i]):
        if spatial_attr[i] >= threshold:
            break
        contig_excl = {i}
        while True:
            size = len(contig_excl)
            contig_excl |= {neigh for ex in contig_excl for neigh in spatial_weights.neighbors[ex] if neigh in excluded}
            if len(contig_excl) == size:
                break
        if sum(spatial_attr[ex] for ex in contig_excl) < threshold:
            excluded.add(i)
    return excluded

def _merge_leaf_nodes(spatial_weights, spatial_attr, sim_mat, index_mapping, threshold):
    adjustment = 0
    while True:
        merge_candidates = [ind for ind,i in enumerate(spatial_attr) if i < threshold and len(spatial_weights.neighbors[ind]) == 1]
        if not merge_candidates:
            return spatial_weights, spatial_attr, sim_mat, index_mapping, adjustment
        to_merge = merge_candidates[0]
        merge_neigh = spatial_weights.neighbors[to_merge][0]
        adjustment += sim_mat[to_merge, merge_neigh]
        mask = np.ones(spatial_weights.sparse.shape[0], dtype=bool)
        mask[to_merge] = False
        spatial_weights = weights.W.from_sparse(spatial_weights.sparse[mask, :][:, mask])
        spatial_attr[merge_neigh] += spatial_attr[to_merge]
        spatial_attr = spatial_attr[mask]
        sim_mat = sim_mat[mask, :][:, mask]
        index_mapping[merge_neigh] |= index_mapping[to_merge]
        index_mapping = {(k-1 if k > to_merge else k):v for k,v in index_mapping.items() if k != threshold}

        
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
    def __init__(self, adj_mat, sim_mat, spatial_attr, threshold, dissimilarity=False):
        # Initialize data
        self.adj_mat = deepcopy(adj_mat)
        self.sim_mat = deepcopy(sim_mat)
        self.spatial_attr = deepcopy(spatial_attr)
        self.threshold = threshold
        self.dissimilarity = dissimilarity

        # Initialize derived information
        self.num_areas = len(self.spatial_attr)
        self.spatial_weights = weights.full2W(self.adj_mat)
        self.regions = np.full(self.num_areas, None, dtype=object)
        self.maxp = None
        self.obj = None
        self.weight_factor = 10**(1 + np.floor(np.log10(np.sum(np.triu(self.sim_mat, k=1)))))
        self._obj_adj = 0

        if self.dissimilarity:
            self._obj_adj += np.sum(np.triu(self.sim_mat, k=1))
            self.sim_mat *= -1

    # geopandas initialization
    @classmethod
    def from_gpd(self):
        pass

    # construct MIP model
    def construct(self, config):
        # Copy input for modification
        copy_spatial_attr = deepcopy(self.spatial_attr)
        copy_sim_mat = deepcopy(self.sim_mat)
        copy_spatial_weights = deepcopy(self.spatial_weights)

        self.index_mapping = {i:{i} for i in range(self.num_areas)}

        if config.merge_leaves:
            copy_spatial_weights, copy_spatial_attr, copy_sim_mat, self.index_mapping, merged_obj = _merge_leaf_nodes(copy_spatial_weights, copy_spatial_attr, copy_sim_mat, self.index_mapping, self.threshold)
            self._obj_adj += merged_obj

        if config.sort_region_roots and config.preassign_roots:
            sort_idx = np.argsort(copy_spatial_attr)[::-1]
            copy_spatial_attr = copy_spatial_attr[sort_idx]
            copy_sim_mat = copy_sim_mat[sort_idx,:][:,sort_idx]
            copy_spatial_weights = weights.W.from_sparse(copy_spatial_weights.sparse[sort_idx,:][:,sort_idx])
            self.index_mapping = {ind:self.index_mapping[i] for ind,i in enumerate(sort_idx)}

        # Define iterative ranges
        self._I_set = range(len(copy_spatial_attr))
        if config.bound_num_regions:
            self._K_set = range(_bound_num_regions(copy_spatial_attr, self.threshold))
        else:
            self._K_set = range(len(copy_spatial_attr))
        if config.bound_contiguity:
            self._C_set = range(_bound_contiguity(copy_spatial_weights, copy_spatial_attr, self.threshold))
        else:
            self._C_set = range(len(copy_spatial_attr))

        # Initialize model
        self.model = LpProblem("MaxP_Model", LpMaximize)

        # Decision variables
        self.x = LpVariable.dicts("var_x", (self._I_set, self._K_set, self._C_set), cat="Binary")
        self.t = LpVariable.dicts("var_t", [(i, j) for i in self._I_set for j in self._I_set if j > i], cat="Binary")

        # Objective function
        self.model += lpSum(self.x[i][j][0] for i in self._I_set for j in self._K_set) * self.weight_factor + lpSum(self.t[i,j] * copy_sim_mat[i][j] for i in self._I_set for j in self._I_set if j > i), "Objective"

        # Define base model constraints
        self.model.extend([ # Single Root Constraints
            lpSum(self.x[i][k][0] for i in self._I_set) <= 1
            for k in self._K_set
        ]) 
        self.model.extend([ # Threshold Constraints
            lpSum(self.x[i][k][c] * copy_spatial_attr[i] for i in self._I_set for c in self._C_set) >= self.threshold * lpSum(self.x[i][k][0] for i in self._I_set)
            for k in self._K_set
        ])
        self.model.extend([ # Single Assignment Constraints
            lpSum(self.x[i][k][c] for k in self._K_set for c in self._C_set) == 1
            for i in self._I_set
        ])
        self.model.extend([ # Adjacency Constraints
            self.x[i][k][c] <= lpSum(self.x[j][k][c-1] for j in copy_spatial_weights.neighbors[i])
            for i in self._I_set 
            for k in self._K_set 
            for c in self._C_set if c > 0
        ])
        if self.dissimilarity:
            self.model.extend([ # x-t Matching Constraints
                self.t[i,j] >= lpSum(self.x[i][k][c] + self.x[j][k][c] for c in self._C_set) - 1
                for i in self._I_set 
                for j in self._I_set if j > i 
                for k in self._K_set
            ])
        else:
            self.model.extend([ # x-t Matching Constraints
                self.t[i,j] <= lpSum(self.x[i][k][c] - self.x[j][k][c] for c in self._C_set) + 1
                for i in self._I_set 
                for j in self._I_set if j > i 
                for k in self._K_set
            ])
        excluded_roots = set()
        if config.exclude_roots: # Exclude Roots
            excluded_roots = _find_excluded_roots(copy_spatial_weights, copy_spatial_attr, self.threshold)
            self.model.extend([
                lpSum(self.x[i][k][0] for k in self._K_set) == 0
                for i in excluded_roots
            ])
        if config.preassign_roots: # Preassign Roots
            if max(copy_spatial_attr) < self.threshold:
                temp_attr = np.array(copy_spatial_attr, copy=True)[list(excluded_roots)] = -np.inf
                self.model += self.x[np.argmax(temp_attr)][0][0] == 1
            else:
                over_thresh = [i for i in self._I_set if copy_spatial_attr[i] >= self.threshold]
                self.model.extend([
                    self.x[i][ind][0] == 1
                    for ind,i in enumerate(over_thresh)
                ])
        if config.max_attr_for_root:
            self.model.extend([
                lpSum(copy_spatial_attr[j] * self.x[j][k][0] for j in self._I_set) >= copy_spatial_attr[i] * self.x[i][k][c]
                for c in self._C_set if c > 0
                for i in self._I_set if i not in excluded_roots
                for k in self._K_set
            ])
        if config.min_adj_order:
            self.model.extend([
                len(copy_spatial_attr) * (1 - self.x[i][k][c]) >= lpSum(self.x[j][k][d] for j in copy_spatial_weights.neighbors[i] for d in range(0,c-1))
                for c in self._C_set if c > 1
                for i in self._I_set
                for k in self._K_set
            ])
        if config.sort_region_roots:
            self.model.extend([
                lpSum(copy_spatial_attr[i] * self.x[i][k-1][0] for i in self._I_set) >= lpSum(copy_spatial_attr[i] * self.x[i][k][0] for i in self._I_set)
                for k in self._K_set if k > 0
            ])

    # solve MIP model
    def solve(self, time):
        self.model.solve(HiGHS(timeLimit=time, msg=True, keepFiles=False, options=['mip_abs_gap=1e-4', 'mip_rel_gap=1e-10']))
        self.maxp = int(value(lpSum(self.x[i][k][0] for i in self._I_set for k in self._K_set)))
        self.obj = value(self.model.objective) + self._obj_adj 

        assigned = {(i,k) for i in self._I_set for k in self._K_set for c in self._C_set if value(self.x[i][k][c]) > 0.9}
        for i,k in assigned:
            for ind in self.index_mapping[i]:
                self.regions[ind] = k
        self.regions = _standardize_solution(self.regions)

    