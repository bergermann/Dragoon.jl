# <img src="docs/img/Dragoon.png" alt="" width=300> <!--Dragoon.jl-->

# Dragoon

Tools and algorithms for optimization of a [MADMAX](https://madmax.mpp.mpg.de/)-like
booster setup.

## Installation
As it is not part of the official Julia package registry, first manually install the
[BoostFractor](https://github.com/mppmu/BoostFractor.jl) package:
```julia
using Pkg
Pkg.add(url="https://github.com/mppmu/BoostFractor.jl.git")
```

Then install the main package:
```julia
Pkg.add(url="https://github.com/bergermann/Dragoon.jl.git")
```

## Usage
See [examples](./examples) for comprehensive guides and read docstrings for additional
information.

Implements three different, customizable optimization algorithm classes, namely
[Nelder-Mead](https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method), 
[Linesearch](https://en.wikipedia.org/wiki/Line_search) and 
[Simulated Annealing](https://en.wikipedia.org/wiki/Simulated_annealing).
Each algorithm was modified to suit the specific needs and properties of a physical
MADMAX-like booster. This package supports analytical 1d and 3d boostfactor calculations
from [BoostFractor](https://github.com/mppmu/BoostFractor.jl).

NOTE: The 3d calculations use distributed computing, hence be sure to do
```julia
julia> using Distributed
julia> addprocs(n) # n is the number of cores/workers you want to use
julia> println(nworkers()) # check if workers are correctly available
julia> @everywhere using BoostFractor, Dragoon
```
before using any 3d function (not necessary if you want to work on a single core (default)). 

For use with an actual physical booster and VNA, additional packages are available:
[Motor Control](https://git.rwth-aachen.de/nick1/XIMC-jl), 
[VNA Control](https://git.rwth-aachen.de/nick1/KeyVNA-jl).

### Modular Options
- Objective functions: `ObjAnalytical`, `ObjRef`, `ObjRefLin`, `ObjRefSquare`, `ObjRefExp`,
    `Obj3d`
- Unstuckinators: `UnstuckDont`, `UnstuckRandom`, `UnstuckCoord`

#### Nelder-Mead
- Simplex Initialization: `InitSimplexCoord`, `InitSimplexAffine`, `InitSimplexRegular`
- Simplex Scheduler: `DefaultSimplexSampler`
- Unstuckinator: `UnstuckDont`, `UnstuckNew`, `UnstuckExpand`

#### Linesearch
- Derivative Calculator: `Derivator1`, `Derivator2`
- Solvers: `SolverSteep`, `SolverNewton`, `SolverBFGS`, `SolverHybrid`
- Step Manager: `StepNorm`
- Searches: `SearchStandard`, `SearchExtendedSteps`, `SearchExtendedDist`, `SearchTrueNewton`

#### Simulated Annealing
- Temperature Manager: `TempLinear`, `TempExp`
- Unstuckinator: `UnstuckDont`, `UnstuckTemp`
