from pulp import *
import highspy
import numpy as np
from libpysal import weights
from dataclasses import dataclass
from scipy.spatial.distance import pdist, squareform
from copy import deepcopy
import geopandas as gpd

# This class constructs and solves the max-p-regions problem using an exact MILP solvers
# The strengthening methods and algorithms are based on my honors thesis at the Univeristy of Minnesota: 
# Strengthening the Max-P-Regions Problem for the Confidentiality of Census Microdata (2026)
# Author: Arlan Hegenbarth

def _bound_num_regions(spatial_attr, threshold):
    """
    Calculates an upper bound for the number of regions in the optimal solution

    Parameters
    ----------

    spatial_attr : ndarray of shape (n,), required
        spatial extensive attribute values for the input areas used for 
        thresholding the regions

    threshold : {int, float}, required
        minimum spatially extensive attribute for each region

    Returns
    -------

    region_bound : int
        An upper bound on the number of regions in the problem instance

    """

    # Find number of areas over the threshold
    over_thres = np.sum(spatial_attr >= threshold)

    # Find number of region that can be created from areas under threshold
    under_thres = np.sum(spatial_attr[spatial_attr < threshold])
    under_thres //= threshold

    region_bound = over_thres + under_thres
    return int(region_bound)


def _can_split(spatial_attr, threshold, path):
    """
    Determines if a path of input areas can be split into two regions that meet the minimum threshold requirement

    Parameters
    ----------

    spatial_attr : ndarray of shape (n,), required
        spatial extensive attribute values for the input areas used for 
        thresholding the regions

    threshold : {int, float}, required
        minimum spatially extensive attribute for each region

    path : list, required
        list of indices representing the path that is being tested

    Returns
    -------

    can_split : boolean
        ``True`` if the path can be split into two parts that meet the minimum region threshold

    """

    attr_list = spatial_attr[path]

    # Calculate cummulative sum from head and tail
    from_head = np.cumsum(attr_list) >= threshold
    from_tail = np.cumsum(attr_list[::-1])[::-1] >= threshold

    # Determine if any split points are above threshold on both sides
    split_points = from_head[:-1] & from_tail[1:]
    can_split = np.any(split_points)
    return can_split


def _recursive_step(spatial_weights, spatial_attr, threshold, path, excluded):
    """
    Take a recursive step on the depth first search to bound the maximum contiguity order

    Parameters
    ----------

    spatial_weights : libpysal.weights.W
        libpysal spatial weights object for input areas

    spatial_attr : ndarray of shape (n,), required
        spatial extensive attribute values for the input areas used for 
        thresholding the regions

    threshold : {int, float}, required
        minimum spatially extensive attribute for each region

    path : list, required
        list of indices representing the path that is being tested.

    excluded : set, required
        list of indices that cannot be added to this path, because they are adjacent
        to an area that is not the head 

    Returns
    -------

    max_q : int
        Maximum path length found by recursively adding areas to this path

    """

    # Check exit conditions
    if _can_split(spatial_attr, threshold, path):
        return len(path) - 1
    if not (set(spatial_weights.neighbors[path[0]]) - excluded): # Dead-end
        return len(path)
    
    # Iterate though head's neighbors and expand the path
    max_q = 0
    for next_ind in spatial_weights.neighbors[path[0]]:
        if next_ind not in excluded:
            depth = _recursive_step(spatial_weights, spatial_attr, threshold, [next_ind] + path, excluded | set(spatial_weights.neighbors[path[0]]))
            if depth > max_q:
                max_q = depth
    return max_q


def _bound_contiguity(spatial_weights, spatial_attr, threshold):
    """
    Find an upper bound for the maximum contiguity order

    Parameters
    ----------

    spatial_weights : libpysal.weights.W
        libpysal spatial weights object for input areas

    spatial_attr : ndarray of shape (n,), required
        spatial extensive attribute values for the input areas used for 
        thresholding the regions

    threshold : {int, float}, required
        minimum spatially extensive attribute for each region

    Returns
    -------

    max_q : int
        Maximum path length found by recursively adding areas in a depth first search
    """
    max_q = 0
    # Iterate through each area as an initial head for a path
    for i in range(len(spatial_attr)):
        depth = _recursive_step(spatial_weights, spatial_attr, threshold, [i], set())
        if depth > max_q:
            max_q = depth
    return max_q


def _standardize_solution(solution):
    """
    Standardizes the solution to use the same index to represent each region, 
    regardless of indices used by MILP solution.

    Parameters
    ----------

    solution : ndarray (n,)
        A region assignment for the max-p-regions problem

    Returns
    -------

    std_solution : ndarray (n,)
        An index standardized max-p-regions problem solution
    """
    id_map = {}
    counter = 0
    for item in solution:
        if item not in id_map: # new region ID
            id_map[item] = counter
            counter += 1
    std_solution = np.array([id_map[s] for s in solution]) # Assign standardized region numbers
    return std_solution


def _find_excluded_roots(spatial_weights, spatial_attr, threshold):
    """
    Finds a set of areas that can excluded from being roots of an optimal solution

    Parameters
    ----------

    spatial_weights : libpysal.weights.W
        libpysal spatial weights object for input areas

    spatial_attr : ndarray of shape (n,), required
        spatial extensive attribute values for the input areas used for 
        thresholding the regions

    threshold : {int, float}, required
        minimum spatially extensive attribute for each region

    Returns
    -------

    excluded : set
        A set of indices representing areas that can be excluded from being roots
    """
    excluded = set()

    # Iterate in ascending order of spatially extensive attribute
    for i in sorted(range(len(spatial_attr)), key=lambda i: spatial_attr[i]):
        if spatial_attr[i] >= threshold:
            break

        contig_excl = {i}

        # Determine attribute size of contiguous excluded area created
        while True:
            size = len(contig_excl)
            contig_excl |= {neigh for ex in contig_excl for neigh in spatial_weights.neighbors[ex] if neigh in excluded}
            if len(contig_excl) == size:
                break

        # Add area only if contiguous excluded area does not exceed the threshold
        if sum(spatial_attr[ex] for ex in contig_excl) < threshold:
            excluded.add(i)

    return excluded


def _merge_leaf_nodes(spatial_weights, spatial_attr, sim_mat, index_mapping, threshold):
    """
    Merges input areas that are below the region threshold and only have one neighbor

    Parameters
    ----------

    spatial_weights : libpysal.weights.W
        Original libpysal spatial weights object for input areas

    spatial_attr : ndarray of shape (n,), required
        Original spatial extensive attribute values for the input areas used for 
        thresholding the regions

    sim_mat : ndarray of shape (n,n), required
        Original non-negative symmetric adjacency matrix representing the similarities 
        (or dissimilarities) between input areas.

    index_mapping : dict, required
        Original dictionary mapping MILP indices to input area indices

    threshold : {int, float}, required
        minimum spatially extensive attribute for each region

    Returns
    -------

    spatial_weights : libpysal.weights.W
        Updated libpysal spatial weights object for input areas

    spatial_attr : ndarray of shape (n,)
        Updated spatial extensive attribute values for the input areas used for 
        thresholding the regions

    sim_mat : ndarray of shape (n,n)
        Updated non-negative symmetric adjacency matrix representing the similarities 
        (or dissimilarities) between input areas.

    index_mapping : dict
        Updated dictionary mapping MILP indices to input area indices

    adjustment : float
        similarity attribute removed, needs to be added to adjust objective function to original
    """
    adjustment = 0
    while True:
        # Determine merge candidates, and exit if none
        merge_candidates = [ind for ind,i in enumerate(spatial_attr) if i < threshold and len(spatial_weights.neighbors[ind]) == 1]
        if not merge_candidates:
            return spatial_weights, spatial_attr, sim_mat, index_mapping, adjustment
        
        # Identify a merge candidate and its neighbor
        to_merge = merge_candidates[0]
        merge_neigh = spatial_weights.neighbors[to_merge][0]
        adjustment += sim_mat[to_merge, merge_neigh] # Adjustment to objective value

        # Create mask to remove merged area
        mask = np.ones(spatial_weights.sparse.shape[0], dtype=bool)
        mask[to_merge] = False

        # Update spatial weights object
        spatial_weights = weights.W.from_sparse(spatial_weights.sparse[mask, :][:, mask])

        # Update spatial attribute array
        spatial_attr[merge_neigh] += spatial_attr[to_merge]
        spatial_attr = spatial_attr[mask]

        # Update similarity matrix
        sim_mat[merge_neigh, :] += sim_mat[to_merge, :]
        sim_mat[:, merge_neigh] += sim_mat[:, to_merge]
        sim_mat = sim_mat[mask, :][:, mask]
        sim_mat[merge_neigh, merge_neigh] = 0

        # Update index mapping
        index_mapping[merge_neigh] |= index_mapping[to_merge]
        index_mapping = {(k-1 if k > to_merge else k):v for k,v in index_mapping.items()}

        
@dataclass
class MaxPConfig:
    """This class defines the max-p regions problem construction configuration.
    Each parameter is a different strategy that can be used to strenghen the problem

    Parameters
    ----------

    bound_num_regions :  boolean
        A tighter upper bound for the number of regions in the optimal solution will be applied

    bound_contiguity : boolean
        A tighter upper bound for the maximum contiguity order will be applied

    merge_leaves : boolean
        Leaf nodes (areas with only one neighbor) will be merged with their neighbor when applicable
    
    preassign_roots : boolean
        Some input areas will be preassigned as roots of certain regions

    exclude_roots : boolean
        Some input areas will excluded from being roots of a region
    
    min_index_for_root : boolean
        The input area with the lowest index must be the region root

    min_adj_order : boolean
        The smallest possible adjacency order for each area must be used

    sort_region_roots : boolean
        The regions must be sorted by ascending index of the roots

    """
    bound_num_regions: bool = False
    bound_contiguity: bool = False
    merge_leaves: bool = False
    preassign_roots: bool = False
    exclude_roots: bool = False
    min_index_for_root: bool = False
    min_adj_order: bool = False
    sort_region_roots: bool = False


class MaxPExact():
    """The max-p-regions problem involves the aggregation of n areas into an unknown
    maximum number of homogeneous regions, while ensuring that each region is contiguous
    and satisfies a minimum threshold value imposed on a predefined spatially extensive
    attribute. This class is designed to solve the max-p-regions problem using a exact
    optimization approach with optimization solvers.

    Parameters
    ----------

    adj_mat : ndarray of shape (n,n), required
        binary symmetric adjacency matrix between input areas.

    sim_mat : ndarray of shape (n,n), required
        non-negative symmetric adjacency matrix representing the similarities 
        (or dissimilarities) between input areas.

    spatial_attr : ndarray of shape (n,), required
        spatial extensive attribute values for the input areas used for 
        thresholding the regions

    threshold : {int, float}, required
        minimum spatially extensive attribute for each region

    dissimilarity : boolean
        Set to ``True`` if sim_mat represents dissimilarities.
        If true, within-region dissimilarity will be minimized.
        If false, within-region similarity will be maxmimized.

    Attributes
    ----------

    num_areas : int
        The number of input areas in the problem

    spatial_weights : libpysal.weights.W
        libpysal spatial weights object for input areas

    regions : ndarray (n,)
        Region assignments for the best solution with standardized indexing

    maxp : int
        The number of regions in the best solution

    obj : float
        The objective value in the best solution

    weight_factor : float
        The weighting factor (10^h) for the objective function

    status : string
        The current solve status of the problem

    optimal : bool
        Optimal solution status, for quick access

    """
    # array initialization
    def __init__(self, adj_mat, sim_mat, spatial_attr, threshold, dissimilarity=False):
        """
        Initializes a max-p-regions problem from numpy arrays
        """
        try:
            # Check input data
            if not isinstance(adj_mat, np.ndarray):
                raise TypeError("Adjacency matrix must be a numpy array")
            if not np.array_equal(adj_mat, adj_mat.T):
                raise ValueError("Adjacency matrix is not symmetric")
            if not np.all((adj_mat == 0) | (adj_mat == 1)):
                raise ValueError("Adjacency matrix is not binary")
            if not isinstance(sim_mat, np.ndarray):
                raise TypeError("Similarity matrix must be a numpy array")
            if not np.array_equal(sim_mat, sim_mat.T):
                raise ValueError("Similarity matrix is not symmetric")
            if np.min(sim_mat) < 0:
                raise ValueError("Similarity matrix must be non-negative")
            if not isinstance(spatial_attr, np.ndarray):
                raise TypeError("Spatial attribute values must be a numpy array")
            if np.min(spatial_attr) < 0:
                raise ValueError("Spatial attribute values must be non-negative")
            if not isinstance(threshold, (int, float, np.number)):
                raise TypeError("Threshold must be a integer or float")
            if not isinstance(dissimilarity, bool):
                raise TypeError("Dissimilarity flag must be a boolean")
            if threshold <= 0:
                raise ValueError("Threshold must be positive")
            if adj_mat.shape[0] != sim_mat.shape[0]:
                raise ValueError("Adjacency matrix and similarity matrix must be the same size")
            if adj_mat.shape[0] != len(spatial_attr):
                raise ValueError("Length of spatial attribute value must match spatial adjacency matrix")
        except (ValueError, TypeError):
            raise

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
        self.optimal = False
        self.weight_factor = 10**(1 + np.floor(np.log10(np.sum(np.triu(self.sim_mat, k=1)))))
        self._obj_adj = 0
        self.status = "Unconstructed"

        if self.dissimilarity: # Use negative values, but adjust the objective result
            self._obj_adj += np.sum(np.triu(self.sim_mat, k=1))
            self.sim_mat *= -1


    # geopandas initialization
    @classmethod
    def from_gdf(cls, gdf, sp_weights, dissim_attr, threshold_attr, threshold):
        """
        Initializes a max-p-regions problem from a geopandas dataframe

        Parameters
        ----------

        gdf : geopandas.GeoDataFrame, required
            Geodataframe containing the original input areas

        sp_weights : libpysal.weights.W, required
            Weights object created from the given geodataframe

        dissim_attr : list, required
            Strings for attribute names (columns of gdf) used for 
            dissimilarity calculations

        threshold_attr : string, required
            The name of the spatial extensive attribute in the gdf

        threshold : {int, float}, required
            minimum spatially extensive attribute for each region
        """

        # Check input data
        if not isinstance(gdf, gpd.GeoDataFrame):
            raise TypeError("gdf must be a geopandas dataframe")
        if not isinstance(sp_weights, weights.W):
            raise TypeError("sp_weights must be a libpysal.weights.W object")
        if len(gdf) != sp_weights.n:
            raise ValueError("Geodataframe and spatial weights must have same number of indices")

        attr = np.atleast_2d(gdf[dissim_attr].values)
        if attr.shape[0] == 1:
            attr = attr.T
        dist_matrix = squareform(pdist(attr, metric="cityblock"))
        threshold_array = gdf[threshold_attr].values
        instance = cls(sp_weights.full()[0], dist_matrix, threshold_array, threshold, dissimilarity=True)
        return instance


    # construct MILP model
    def construct(self, config):
        """
        Constructs a MILP formulation for the max-p-regions problem instance

        Parameters
        ----------

        config: maxp_exact.MaxPConfig object, required
            A MaxPConfig object for the strategies used to construct the MILP Problem
        """
        # Copy input for modification
        copy_spatial_attr = deepcopy(self.spatial_attr)
        copy_sim_mat = deepcopy(self.sim_mat)
        copy_spatial_weights = deepcopy(self.spatial_weights)

        self._index_mapping = {i:{i} for i in range(self.num_areas)}

        # Merge leaf nodes
        if config.merge_leaves:
            copy_spatial_weights, copy_spatial_attr, copy_sim_mat, self._index_mapping, merged_obj = _merge_leaf_nodes(copy_spatial_weights, copy_spatial_attr, copy_sim_mat, self._index_mapping, self.threshold)
            self._obj_adj += merged_obj

        # If both sort region roots and preassign roots are used, all data structures must be sorted
        if config.sort_region_roots and config.preassign_roots:
            sort_idx = np.argsort(copy_spatial_attr)[::-1]
            copy_spatial_attr = copy_spatial_attr[sort_idx]
            copy_sim_mat = copy_sim_mat[sort_idx,:][:,sort_idx]
            copy_spatial_weights = weights.W.from_sparse(copy_spatial_weights.sparse[sort_idx,:][:,sort_idx])
            self._index_mapping = {ind:self._index_mapping[i] for ind,i in enumerate(sort_idx)}

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
            self.model.extend([ # x-t Matching Constraints for dissimilarity
                self.t[i,j] >= lpSum(self.x[i][k][c] + self.x[j][k][c] for c in self._C_set) - 1
                for i in self._I_set 
                for j in self._I_set if j > i 
                for k in self._K_set
            ])
        else:
            self.model.extend([ # x-t Matching Constraints for similarity
                self.t[i,j] <= lpSum(self.x[i][k][c] - self.x[j][k][c] for c in self._C_set) + 1
                for i in self._I_set 
                for j in self._I_set if j > i 
                for k in self._K_set
            ])

        # Define constraints for strengthened formulations
        excluded_roots = set()
        if config.exclude_roots: # Exclude Roots Constraints
            excluded_roots = _find_excluded_roots(copy_spatial_weights, copy_spatial_attr, self.threshold)
            self.model.extend([
                lpSum(self.x[i][k][0] for k in self._K_set) == 0
                for i in excluded_roots
            ])
        if config.preassign_roots: # Preassign Roots Constraints
            if max(copy_spatial_attr) < self.threshold: # No areas above regional attribute threshold
                temp_attr = np.array(copy_spatial_attr, copy=True)
                temp_attr[list(excluded_roots)] = -1
                self.model += self.x[np.argmax(temp_attr)][0][0] == 1
            else: # Some areas above regional attribute threshold
                self.model.extend([
                    self.x[i][ind][0] == 1
                    for i,ind in enumerate(self._I_set) if copy_spatial_attr[i] >= self.threshold
                ])
        if config.min_index_for_root:
            self.model.extend([ # Ensure Root is Minimum Index Constraints
                lpSum((len(copy_spatial_attr) - j) * self.x[j][k][0] for j in self._I_set) >= (len(copy_spatial_attr) - i) * self.x[i][k][c]
                for c in self._C_set if c > 0
                for i in self._I_set if i not in excluded_roots
                for k in self._K_set
            ])
        if config.min_adj_order:
            self.model.extend([ # Ensure Minimum Possible Adjacency Order Contraints
                len(copy_spatial_attr) * (1 - self.x[i][k][c]) >= lpSum(self.x[j][k][d] for j in copy_spatial_weights.neighbors[i] for d in range(0,c-1))
                for c in self._C_set if c > 1
                for i in self._I_set
                for k in self._K_set
            ])
        if config.sort_region_roots:
            self.model.extend([ # Sort Regions by Root Index Constraints
                lpSum((len(copy_spatial_attr) - i) * self.x[i][k-1][0] for i in self._I_set) >= lpSum((len(copy_spatial_attr) - i) * self.x[i][k][0] for i in self._I_set)
                for k in self._K_set if k > 0
            ])

        # Save status as no solution
        self.status = LpSolution[0]

    # solve MIP model
    def solve(self, time, abs_gap=1e-7, rel_gap=1e-4):
        """
        Solve a constructed MILP for a max-p-region problem formulation

        Parameters
        ----------

        time : {int, float}, required
            Length of maximum solve time in seconds

        abs_gap : float
            Absolute gap at which the MILP is considered solved to optimality

        rel_gap : float
            Relative gap at which the MILP is considered solved to optimality
        """
        self.model.solve(HiGHS(timeLimit=time, msg=True, keepFiles=False, options=[f'mip_abs_gap={abs_gap}', f'mip_rel_gap={rel_gap}']))

        # Derive Result
        self.maxp = int(value(lpSum(self.x[i][k][0] for i in self._I_set for k in self._K_set)))
        self.obj = value(self.model.objective) + self._obj_adj 
        self.status = LpSolution[self.model.sol_status]
        if self.status == "Optimal Solution Found":
            self.optimal = True

        # Extract and Standardize Solution
        assigned = {(i,k) for i in self._I_set for k in self._K_set for c in self._C_set if value(self.x[i][k][c]) > 0.9}
        for i,k in assigned:
            for ind in self._index_mapping[i]:
                self.regions[ind] = k
        self.regions = _standardize_solution(self.regions)

    