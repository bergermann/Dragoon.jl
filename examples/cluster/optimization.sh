### Script to start optimization jobs on cluster. First, a file algorithmvariant.jl needs to be prepared
### if not already present.
###
### Usage:
### sh optimization.sh algorithmvariant sigx nsig [fcenter fwidth npoints ndisk eps tand]
###
### Example:
### sh optimization.sh nm1 100e-6 1000 20.025e9 50e6 10 20 24.0 0.
###
### parameters in brackets are optional and standard parameters (same as in example) are chosen if omitted
### parameter order needs to be retained, singe parameters can be omitted with _, e.g.
### sh optimization.sh nm1 100e-6 1000 _ 50e6 100 _


### add memory/time scaling in future?

# cd /home/jn226467/

if ! test -f /optimizers/${1}.jl; then
    echo "Optimizer file ${1}.jl does not exist. Aborting."
    exit 1
fi

if [[ "$#" -lt 3 ]]; then
    echo "Invalid number of arguments. At least 3 required."
    echo "ARGS: algorith sigx nsig [fcenter fwidth npoints ndisk eps tand]"
    exit 1
fi

NAME=opt_${1}

echo "Jobname: $NAME"

sbatch -J $NAME scheduler.sh "$@"