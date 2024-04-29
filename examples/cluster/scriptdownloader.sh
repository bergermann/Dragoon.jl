

### get optimizers and stuff
wget "https://raw.githubusercontent.com/bergermann/Dragoon.jl/main/examples/cluster/standard_settings.jl" -O "optimizers/standard_settings.jl"
wget "https://raw.githubusercontent.com/bergermann/Dragoon.jl/main/examples/cluster/variants.txt" -O "optimizers/variants.txt"

wget "https://raw.githubusercontent.com/bergermann/Dragoon.jl/main/examples/cluster/nm1.jl" -O "optimizers/nm1.jl"



### get cluster scheduler scripts
wget "https://raw.githubusercontent.com/bergermann/Dragoon.jl/main/examples/cluster/optimization.sh" -O "optimization.sh"
wget "https://raw.githubusercontent.com/bergermann/Dragoon.jl/main/examples/cluster/scheduler.sh" -O "scheduler.sh"

### update Dragoon
julia -e "using Pkg; Pkg.add(url=\"https://github.com/bergermann/Dragoon.jl.git\"); Pkg.update()"