

### get optimizers and stuff
wget "https://github.com/bergermann/Dragoon.jl/blob/main/examples/cluster/standard_settings.jl" -O "optimizers/standard_settings.jl"
wget "https://github.com/bergermann/Dragoon.jl/blob/main/examples/cluster/variants.txt" -O "optimizers/variants.txt"

wget "https://github.com/bergermann/Dragoon.jl/blob/main/examples/cluster/nm1.jl" -O "optimizers/nm1.jl"



### get cluster scheduler scripts
wget "https://github.com/bergermann/Dragoon.jl/blob/main/examples/cluster/optimization.sh" -O "optimization.sh"
wget "https://github.com/bergermann/Dragoon.jl/blob/main/examples/cluster/scheduler.sh" -O "scheduler.sh"

### update Dragoon
julia -e "using Pkg; Pkg.add(\"https://github.com/mppmu/BoostFractor.jl.git\"); Pkg.update()"