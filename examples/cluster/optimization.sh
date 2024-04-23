### sh optimization.sh algorithmvariant sigx nsig [fcenter fwidth npoints ndisk eps tand]
### sh optimization.sh nm1 100e-6 1000 [20.025e9 50e6 10 20 24.0 0.]

NAME=opt_${1}_${2}
OUT=out_${1}_${2}

echo "Jobname: $NAME"
echo "Output:  $OUT"

sh --jobname=$NAME --output=$OUT scheduler.sh "$@"