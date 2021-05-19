#!/bin/bash


PINBALLEX=$@


if [ -z "$PINBALLEX" ]; then
	echo "Please provide the path to the executable to the command line, for example:"
	echo "  ./run.sh ~/bin/pinball.ex"
	echo "  ./run.sh mpirun -np 8 ~/git/pinball/espresso-5.2.0/bin/pw.x"
	exit
fi
# running scf
cd scf-delithiated
$PINBALLEX < scf.in | tee scf.out
cd ..

cp -r scf-delithiated/out pinball-nofit/out


cd pinball-nofit 
$PINBALLEX < pinball.in | tee pinball.out
cd ..


