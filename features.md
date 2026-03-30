# List of Desired Features

### Basic Functionality
- [x] Intake data as vectors
- [x] Construct a model
- [x] Solve a model

### Beginning Strenthening Strategies from Thesis
- [x] Upper bound on number of regions
- [x] Upper bound on contiguity order
- [ ] Merge Leaf Nodes
- [x] Preassign Roots
- [ ] Exclude Roots
- [ ] Sort Regions
- [ ] Maximize Root Node
- [ ] Minimize Adjacency Order

### Additional Strengthening Strategies from Thesis
- [ ] Split independent subproblems
- [ ] Reduce size of t
- [ ] Double bound t
- [ ] Lower bound solution
- [ ] Inverse objective

### Broad Application
- [ ] Use similarity or disimilarity
- [ ] Non-adjacent zero similarity tracking
- [ ] Connect non-adjacent parts
- [ ] Gracefully handle memory issues

### User Inputs
- [ ] Input gap tolerance
- [ ] Ability to input geopandas
- [ ] Ability to write to geopandas column
- [ ] Automatically calculate dissimilarity from attribute
- [ ] Automatically calculate border lengths
- [ ] Specifiy type of contiguity

### Accessing the Solution
- [ ] Seperately store optimal and best known solution
- [ ] Seperately store optimal and best known objective value
- [ ] Store problem status
- [ ] Access to MPS file

### Solve Tracking
- [ ] Store solver properties as class attributes
- [ ] Access to converted upper bound
- [ ] Access to other model attributes
- [ ] Temporal tracking
- [ ] Calculation of primal and dual integrals
- [ ] Custom callback functions

### Optimization Connection
- [ ] Replace pulp with pyomo
- [ ] Allow selection between multiple solvers
- [ ] Allow any solver path including commercial
- [ ] Allow disjoint solving where progress is saved

### Extra Reformulations not from Thesis
- [ ] Independent splits with non-adjacent similarities
- [ ] Parallelize independent splits
- [ ] Convert t to continuous <= 1
- [ ] Two phase solve
- [ ] Iterative bounding of similarities
- [ ] Seed with heuristic

### Testing
- [ ] Full solves feasible, infeasible, time limit, memory limit
- [ ] Individually test each algorithm
- [ ] Test all combinations of interacting methods