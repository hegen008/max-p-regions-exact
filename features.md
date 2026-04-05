# List of Desired Features

### Basic Functionality
- [x] Intake data as vectors
- [x] Construct a model
- [x] Solve a model

### Beginning Strenthening Strategies from Thesis
- [x] Upper bound on number of regions
- [x] Upper bound on contiguity order
- [x] Merge Leaf Nodes
- [x] Preassign Roots
- [x] Exclude Roots
- [x] Sort Regions
- [x] Maximize Root Node
- [x] Minimize Adjacency Order

### Additional Strengthening Strategies from Thesis
- [ ] Split independent subproblems
- [ ] Reduce size of t
- [ ] Double bound t
- [ ] Lower bound solution
- [ ] Inverse objective

### Broad Application
- [x] Use similarity or disimilarity - **MVP**
- [ ] Non-adjacent zero similarity tracking
- [ ] Connect non-adjacent parts
- [ ] Gracefully handle memory issues
- [ ] Handle default and trivial solutions
- [ ] Weight similarity to always be right of decimal
- [ ] Index sorts rather than attribute sorts

### User Inputs
- [x] Input gap tolerance - **MVP**
- [x] Ability to input geopandas - **MVP**
- [x] Ability to write to geopandas column - **MVP**
- [ ] Use integer or float as threshold attribute
- [ ] Error checking in the init function
- [x] Automatically calculate dissimilarity from attribute
- [ ] Automatically calculate border lengths
- [x] Specifiy type of contiguity
- [ ] Store Unique IDs

### Accessing the Solution
- [ ] Seperately store optimal and best known solution
- [ ] Seperately store optimal and best known objective value
- [x] Store problem status - **MVP**
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
- [ ] Both types of sorting regions
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
- [ ] Add example data to repo
