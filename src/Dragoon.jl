# This file is a part of Dragoon.jl, licensed under the MIT License (MIT).

__precompile__(true)

module Dragoon

using BoostFractor

import Dates: DateTime, now, UTC, canonicalize, Second, Nanosecond
import DateFormats: /ₜ, *ₜ

# General optimizers
include("types.jl")
include("functions.jl")
include("optimizers.jl")
include("objective.jl")
include("unstuckinators.jl")

include("Analytical1D.jl")
include("utilities.jl")

# Linesearch options
include("linesearch/functions.jl")
include("linesearch/solvers.jl")
include("linesearch/derivatives.jl")
include("linesearch/steps.jl")
include("linesearch/searches.jl")

# Nelder-Mead options
include("neldermead/functions.jl")
include("neldermead/simplexinit.jl")
include("neldermead/simplexobj.jl")

# Simulated Annealing options
include("annealing/functions.jl")

end
