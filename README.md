# admm-nonconvex-heuristic
This script implements a heuristic for solving a mixed-integer quadratic programming problem using ADMM. The heuristic is presented in Takapoui et al. (2015). A Simple Effective Heuristic for Embedded Mixed-Integer Quadratic Programming.

More specifically, the sample problem has a chain of n nodes. Each node has a load and a generator. The generators can be turned on or off. The model seeks to meet the load at each node by minimizing the cost of generation. Depending on the cost parameters of the generators, the nodes can import from or export to neighbouring nodes. Mathematically, the heuristic solves the following problem

```
Minimize (1/2)x'Qx + c'x + alpha
s.t. Ax = b
x are continuous or integer
```

# Running the sample model
Run model.m block-by-block in Matlab. You need to have Gurobi installed and added to your Matlab path to compare the results of the heuristic to the optimal solution. The heuristic itself does not require Gurobi.

