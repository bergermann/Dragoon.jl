
### get optimizers and stuff
wget "https://raw.githubusercontent.com/bergermann/Dragoon.jl/main/examples/cluster/standard_settings.jl" -O "optimizers/standard_settings.jl"
wget "https://raw.githubusercontent.com/bergermann/Dragoon.jl/main/examples/cluster/variants.txt" -O "optimizers/variants.txt"

wget "https://raw.githubusercontent.com/bergermann/Dragoon.jl/main/examples/cluster/optimizers/nm1.jl" -O "optimizers/nm1.jl"
wget "https://raw.githubusercontent.com/bergermann/Dragoon.jl/main/examples/cluster/optimizers/nm1ref.jl" -O "optimizers/nm1ref.jl"

wget "https://raw.githubusercontent.com/bergermann/Dragoon.jl/main/examples/cluster/optimizers/sa1.jl" -O "optimizers/sa1.jl"

### get cluster scheduler scripts
wget "https://raw.githubusercontent.com/bergermann/Dragoon.jl/main/examples/cluster/optimization.sh" -O "optimization.sh"
wget "https://raw.githubusercontent.com/bergermann/Dragoon.jl/main/examples/cluster/scheduler.sh" -O "scheduler.sh"

### update Dragoon
/home/jn226467/julia-1.10.2/bin/julia -e "using Pkg; Pkg.add(url=\"https://github.com/bergermann/Dragoon.jl.git\"); Pkg.update()"