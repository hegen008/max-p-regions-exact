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
If there max-p-regions problem is being formulated to minimize within region similarity, the following change are made

### Parameters
$$d_{ij} = \text{ similarity relationship between areas } i \text{ and } j \text{, with } i, j \in I \text{ and } i<j$$

### Objective Function

Maximize

$$Z = \sum_{k \in K}\sum_{i \in I}x_i^{k0} \cdot 10^h - \sum_{i \in I}\sum_{j \in I \mid j>i}d_{ij}t_{ij}$$

### Constraints

Ensure alignment of $x$ and $t$

$$t_{ij} \geq \sum_{c \in C}x_i^{kc} + \sum_{c \in C}x_j^{kc} - 1  \quad \forall i \in I, j \in I, k \in K \mid j > i$$

## Strengthening Methods

### Bound Number of Regions
Controlled by the parameter `bound_num_regions`

### Bound Maximum Contiguity Order
Controlled by the parameter `bound_contiguity`

### Merge Leaf Nodes
Controlled by the parameter `merge_leaves`

### Preassign Root Areas
Controlled by the parameter `preassign_roots`

### Exclude Areas as Roots
Controlled by the parameter `exclude_roots`

### Max Attribute for Root Area
Controlled by the parameter `max_attr_for_root`

### Minimize Adjacency Order
Controlled by the parameter `min_adj_order`

### Sort Region Roots
Controlled by the parameter `sort_region_roots`