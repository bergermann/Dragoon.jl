### sh optimization.sh algorithmvariant sigx nsig [fcenter fwidth npoints ndisk eps tand]
### sh optimization.sh nm1 100e-6 1000 [20.025e9 50e6 10 20 24.0 0.]

### add memory/time scaling in future?

cd /home/jn226467/

NAME=opt_${1}

echo "Jobname: $NAME"

sbatch -J $NAME scheduler.sh "$@"