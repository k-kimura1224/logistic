#!/bin/sh
run()
{
   buf=${1##*/}
   filename=${buf%.*}

   ../bin/GENE_POINTS $1 9 ../../sw_solution/${filename}_b.sol ../../sw_solution/${filename}_f.sol > ./${filename}_09.lp
   ../bin/GENE_POINTS $1 17 ../../sw_solution/${filename}_b.sol ../../sw_solution/${filename}_f.sol > ./${filename}_17.lp
}

datadir=$1

rm -f *.lp

for dir in `ls -d ${datadir}/*.logreg`
do
   run $dir
done

