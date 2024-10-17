# This file is a part of Dragoon.jl, licensed under the MIT License (MIT).

__precompile__(true)

"""
    Dragoon

MADMAX disc position optimization package. `Dragoon` contains multiple MADMAX specific
implementations of optimization algorithms, namely `nelderMead`, `linesearch` and
`simulatedAnnealing`. Each algorithm features high customizability through modular options.

See [`here`](insertlinkformasterthesis), mainly chapter 5, for a comprehensive guide.
"""
module Dragoon

using BoostFractor
using Distributed
using LinearAlgebra
using Plots
using LaTeXStrings

import Dates: DateTime, now, UTC, canonicalize, format,
    Year, Month, Day, Hour, Minute, Second, Nanosecond
import DateFormats: /ₜ, *ₜ
import ProgressBars: ProgressBar
import Peaks: findmaxima

# General optimizers
include("types.jl")
include("helper_functions.jl")
include("objectivefunctions.jl")
include("unstuckinators.jl")

include("Analytical1D.jl")
include("general_utilities.jl")

# Linesearch options
include("linesearch/optimizer.jl")
include("linesearch/functions.jl")
include("linesearch/solvers.jl")
include("linesearch/derivatives.jl")
include("linesearch/steps.jl")
include("linesearch/searches.jl")

# Nelder-Mead options
include("neldermead/optimizer.jl")
include("neldermead/functions.jl")
include("neldermead/simplexinit.jl")
include("neldermead/simplexobj.jl")
include("neldermead/unstuckinators.jl")

# Simulated Annealing options
include("annealing/optimizer.jl")
include("annealing/functions.jl")
include("annealing/tempmanagers.jl")
include("annealing/unstuckinators.jl")

include("dragoon_.jl")

end
