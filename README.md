# Max-P-Regions Problem

The max-p-regions problem clusters "a set of geographic areas into the maximum number of homogeneous regions such that the value of a spatially extensive regional attribute is above a predefined threshold value" (Duque 2012). This repository implements an exact solver for the max-p-regions problem and strengthens to formulation for larger problem instances to be solved using exact methods.

## Similarity Problem Formulation

### Parameters

$$i, j, I = \text{ indices and set of input areas, } I = \{1,\dots,n\}$$

$$k, K = \text{ index and set of potential regions, } K = \{1,\dots,m\} \text{, with } m=n \text{ by default }$$

$$c, C = \text{ index and set of contiguity orders, } C = \{0,\dots,q\} \text{, with } q = n-1 \text{ by default }$$

$$w_{ij} = \begin{cases} 
1, & \text{if areas } i \text{ and } j \text{ are adjacent, with } i, j \in I \text{ and } i \neq j \\ 
0, & \text{otherwise} 
\end{cases}$$

$$N_i = \{j|w_{ij}=1\} \text{, the set of areas that are adjacent to area } i \text{, this is an alternate representation of } w$$

$$s_{ij} = \text{ similarity relationship between areas } i \text{ and } j \text{, with } i, j \in I \text{ and } i<j$$

$$h  = 1 + \lfloor log(\sum_i \sum_{j \mid j>i} s_{ij})\rfloor \text{, which is the number of integer digits of} \sum_i \sum_{j \mid j>i} s_{ij}\text{, with } i,j \in I $$

$$l_i = \text{ spatially extensive attribute value of area } i \text{, with } i \in I$$

$$\tau  = \text{ threshold, minimum value for attribute } l \text{ at regional scale} $$

### Decision Variables

$$t_{ij}  = \begin{cases}
1, & \text{ if areas } i \text{ and } j \text{ belong to the same region, with } j > i \\
0, & \text{otherwise}
\end{cases}$$

$$x_i^{kc} = \begin{cases}
1, & \text{ if area } i \text{ is assigned to region } k \text{ in contiguity order } c \\
0, & \text{otherwise}
\end{cases}$$

### Objective Function

Maximize

$$Z = \sum_{k \in K}\sum_{i \in I}x_i^{k0} \cdot 10^h + \sum_{i \in I}\sum_{j \in I \mid j>i}s_{ij}t_{ij}$$

### Constraints

Regions cannot have multiple roots

$$\sum_{i \in I}x_i^{k0} \leq 1 \quad \forall k \in K$$

Each are is assigned to exactly one region

$$\sum_{k \in K}\sum_{c \in C}x_i^{kc} = 1 \quad \forall i \in I$$

Contiguity constraints

$$x_i^{kc}\leq\sum_{j\in N_i}x_j^{k(c-1)} \quad \forall i \in I,k \in K, c \in C \mid c  > 0$$

Regions meet minimum attribute threshold

$$\sum_{i \in I}\sum_{c \in C}x_i^{kc}l_i \geq \tau \cdot \sum_{i \in I}x_i^{k0} \quad \forall k\in K$$

Ensure alignment of $x$ and $t$

$$t_{ij} \leq \sum_{c \in C}x_i^{kc} - \sum_{c \in C}x_j^{kc} + 1  \quad \forall i \in I, j \in I, k \in K \mid j > i$$

Binary constraints for $x$ and $t$

$$x_i^{kc} \in \{0,1\} \quad \forall i \in I, k \in K, c \in C$$

$$t_{ij} \in \{0,1\} \quad \forall i \in I, j \in I \mid j > i$$


## Dissimilarity Problem Formulation
If the max-p-regions problem is being formulated to minimize within-region dissimilarity, the following changes are made

### Parameters
$$d_{ij} = \text{ dissimilarity relationship between areas } i \text{ and } j \text{, with } i, j \in I \text{ and } i<j$$

### Objective Function

Maximize

$$Z = \sum_{k \in K}\sum_{i \in I}x_i^{k0} \cdot 10^h - \sum_{i \in I}\sum_{j \in I \mid j>i}d_{ij}t_{ij}$$

### Constraints

Ensure alignment of $x$ and $t$

$$t_{ij} \geq \sum_{c \in C}x_i^{kc} + \sum_{c \in C}x_j^{kc} - 1  \quad \forall i \in I, j \in I, k \in K \mid j > i$$

## Strengthening Methods

### Bound Number of Regions
Controlled by the parameter `bound_num_regions`

This method provides an tighter upper bound for the number of regions $m$, with the equation below. Finding a smaller value of $m$ reduce the dimensionality of the optimization problem.

$$m = \sum_{i \in I} (\mathbb{𝟙}\{l_i \geq \tau\}) + \left\lfloor \frac{\sum_{i \in I} (l_i \cdot \mathbb{𝟙}\{l_i < \tau\})}{\tau} \right\rfloor$$

In this equation, $\mathbb{𝟙}\{A\}$ represents an indicator function, where the value is 1 if condition $A$ is met, otherwise the value is 0. The first term of this equation is the total number of areas in the problem that meet the threshold constraint on their own. The second term is the maximum number of regions that could be formed out of the total remaining areas that are under the threshold.

### Bound Maximum Contiguity Order
Controlled by the parameter `bound_contiguity`

This method finds a tighter bound for the maximium contiguity order $q$, again reducing the dimensionality of the optimization problem. The bound can be determined by finding the longest acyclic path of input areas that cannot be split into multiple regions. In this implementation, this is accomplished through and exhaustive breadth-first search of the adjacency graph. In very large problem instances, this is not computationally feasible. While these cases are also likely too large to solve using exact methods, an alternate aspatial bound is possible. At this time, an aspatial bound is left for future improvement.

### Merge Leaf Nodes
Controlled by the parameter `merge_leaves`

This method performs a preprocessing step where some input areas are merged together before model construction. For ayd input area $i$ such that $l_i < \tau$ and $|N_i| = 1$, the input area can be merged with its only neighbor. In all feasible problem solutions these two areas must be assigned to the same region

### Preassign Root Areas
Controlled by the parameter `preassign_roots`

A region in the optimal solution will never contain muliple areas that have a spatially extensive attribute the meets the region threshold. This means thateach input area $i$ that meets the threshold can be assigned as a root of their own region $k$ with the constraint $x_i^{k0} = 1$. In the case where no areas meet the threshold, one arbitrary area can be assigned as a region root. In this implementation, the area with the largest spatially extensive attribute will be preassigned.

### Exclude Areas as Roots
Controlled by the parameter `exclude_roots`

Similarity to preassigning roots, some areas can be excluded from being roots. We know that for any contiguous set of input areas with a total spatially extensive attribute below the threshold, a version of the optimal solution exists where none of these areas are the root of any region. Multiple contiguuous sets of these input areas can exist in the same problem. In this implementation, the set of excluded roots $X$ is identified through a greedy approach. To exclude set $X$ from being region roots, the following constraints can be added.

$$\sum_{k \in K}x_i^{k0} = 0 \quad \forall i \in X$$

### Minimize Index for Root Area
Controlled by the parameter `min_index_for_root`

### Minimize Contiguity Order
Controlled by the parameter `min_cont_order`

### Sort Region Roots
Controlled by the parameter `sort_region_roots`
