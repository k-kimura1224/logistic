#!/bin/sh
#chmod u+x run.sh 
run()
{
	name=$1
   echo " -- running ${name} -- "
	../read_sol/bin/scip -s ../read_sol/settings/default.set -f ../data/${name}.logreg > ./LOG/${name}.log
}

run madelon_train
run musk2
#run pima-indians-diabetes
#
#run mammo
#
#run seismic-bumps
#
#run parkin
#
#run statH
#
#run breast
#run biodeg
#run spectf
#run statG
