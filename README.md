# Computation of Convex Hull Prices using Column-And-Row Generation method.
This github repository you can find three methods to compute the Convex Hull Prices (CHP).
- [Exteded Formulation](https://doi.org/10.1287/ijoc.2017.0802)
- [Row Generation method](https://doi.org/10.1287/ijoc.2017.0802)
- [Column Generation method](https://arxiv.org/abs/2012.13331)
- [Column-And-Row Generation method](https://doi.org/10.1007/s13675-013-0009-9)

The code is in [julia](https://julialang.org/downloads/).
Required packages are : CSV, DataFrames, JuMP, Gurobi, Plots, GLPK, TickTock, PyPlot, LaTeXStrings, JSON.

It can be installed in the julia shell by running these commands

```
using Pkg;
Pkg.add("package_to_add");
```
# Annotations
Row Generation method : RG  
Column Generation method : CG  
Column-And-Row Generation method : CRG