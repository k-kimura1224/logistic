#!/bin/sh
#chmod u+x run.sh 
run()
{
	name=$1
   echo " -- running ${name} -- "
	../../read_sol/bin/scip -s ../../read_sol/settings/default.set -f ../../data/${name}.logreg > ./${name}.log
}

run seismic-bumps
run breast
#run biodeg
run spectf
#run musk2
#run statG
run madelon_train
