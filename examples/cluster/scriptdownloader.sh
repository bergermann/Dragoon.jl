
DIR = "https://raw.githubusercontent.com/bergermann/Dragoon.jl/main/examples/cluster/"

### get optimizers and stuff
wget "${DIR}standard_settings.jl" -O "optimizers/standard_settings.jl"
wget "${DIR}variants.txt" -O "optimizers/variants.txt"

wget "${DIR}optimizers/nm1.jl" -O "optimizers/nm1.jl"



### get cluster scheduler scripts
wget "${DIR}optimization.sh" -O "optimization.sh"
wget "${DIR}scheduler.sh" -O "scheduler.sh"

### update Dragoon
julia -e "using Pkg; Pkg.add(url=\"https://github.com/bergermann/Dragoon.jl.git\"); Pkg.update()"