# Max-P-Regions Problem

The max-p-regions problem clusters "a set of geographic areas into the maximum number of homogeneous regions such that the value of a spatially extensive regional attribute is above a predefined threshold value" (Duque 2012). This repository implements an exact solver for the max-p-regions problem and strengthens to formulation for larger problem instances to be solved using exact methods.

## Problem Formulation

### Parameters

$i, j, I = \text{ indices and set of input areas, } I = \{1,\dots,n\}$

$k, K = \text{ index and set of potential regions, } K = \{1,\dots,m\} \text{, with } m=n \text{ by default }$

$c, C = \text{ index and set of contiguity orders, } C = \{0,\dots,q\} \text{, with } q = n-1 \text{ by default }$

$w_{ij} = \begin{cases} 
1, & \text{if areas } i \text{ and } j \text{ are adjacent, with } i, j \in I \text{ and } i \neq j \\ 
0, & \text{otherwise} 
\end{cases}$

$N_i = \{j|w_{ij}=1\} \text{, the set of areas that are adjacent to area } i \text{, this is an alternate representation of } w$

$s_{ij} = \text{ similarity relationship between areas } i \text{ and } j \text{, with } i, j \in I \text{ and } i<j$

$h  = 1 + \lfloor log(\sum_i \sum_{j \mid j>i} s_{ij})\rfloor \text{, which is the number of integer digits of} \sum_i \sum_{j \mid j>i} s_{ij}\text{, with } i,j \in I $

$l_i = \text{ spatially extensive attribute value of area } i \text{, with } i \in I$

$\tau  = \text{ threshold, minimum value for attribute } l \text{ at regional scale} $

### Decision Variables

$t_{ij}  = \begin{cases}
1, & \text{ if areas } i \text{ and } j \text{ belong to the same region, with } j > i \\
0, & \text{otherwise}
\end{cases}$

$x_i^{kc} = \begin{cases}
1, & \text{ if area } i \text{ is assigned to region } k \text{ in contiguity order } c \\
0, & \text{otherwise}
\end{cases}$

### Objective Function

Maximize

$Z = \sum_{k \in K}\sum_{i \in I}x_i^{k0} \cdot 10^h + \sum_{i \in I}\sum_{j \in I \mid j>i}s_{ij}t_{ij}$

### Constraints

$\sum_{i \in I}x_i^{k0} \leq 1 \quad \forall k \in K$

$\sum_{k \in K}\sum_{c \in C}x_i^{kc} = 1 \quad \forall i \in I$

$x_i^{kc}\leq\sum_{j\in N_i}x_j^{k(c-1)} \quad \forall i \in I,k \in K, c \in C \mid c  > 0$

$\sum_{i \in I}\sum_{c \in C}x_i^{kc}l_i \geq \tau \cdot \sum_{i \in I}x_i^{k0} \quad \forall k\in K$

$t_{ij} \leq \sum_{c \in C}x_i^{kc} - \sum_{c \in C}x_j^{kc} + 1  \quad \forall i \in I, j \in I, k \in K \mid j > i$

$x_i^{kc} \in \{0,1\} \quad \forall i \in I, k \in K, c \in C$

$t_{ij} \in \{0,1\} \quad \forall i \in I, j \in I \mid j > i$
