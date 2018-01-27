#!/bin/sh

#cd exp_cplex4
#~/ILOG/CPLEX_Studio1263/cplex/bin/x86-64_linux/cplex < run.txt
#cd ..
#
#cd ug5/exp
#export OMP_NUMTHREADS=1
#../bin/fscip ../settings/large.set ../../data/musk2.logreg -q -sth 16 > ./musk2.log
#export OMP_NUMTHREADS=1
#../bin/fscip ../settings/noheur.set ../../data/musk2.logreg -q -sth 16 > ./musk2.log
#cd ../..

cd ./stepwise/exp
./run ../../data
cd ../..


cd ./ug5/exp_split
./run2 16
cd ../..
