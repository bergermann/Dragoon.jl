#!/bin/bash

#SBATCH --account=

#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2000M

cd /home/jn226467"
/home/jn226467/julia-1.10.2/bin/julia -p 47 /home/jn226467/${SCRIPT}