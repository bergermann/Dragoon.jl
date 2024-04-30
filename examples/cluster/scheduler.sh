#!/bin/bash

#SBATCH --output=./joboutputs/%x_%j.out
#SBATCH --error=./joboutputs/%x_%j.err

#SBATCH --time=02:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 47
#SBATCH --mem-per-cpu 1G

julia -p 47 optimizers/${1}.jl "${@:2}"