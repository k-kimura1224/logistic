#!/bin/sh

#pred()
#{
#   i=$1
#   dirname=$2
#   subdirname=$3
#   subdir=$4
#   samplefile=${i}_sample_${dirname}.linereg
#   predictfile=${i}_predict_${dirname}.linereg
#   solfile=${i}_sample_${dirname}.sol
#   outputfile=${i}_${dirname}.pred
#   #buf=${1##*/}
#   #filename=${buf%.*}
#
#   echo "../bin/prediction ${subdir}/${samplefile} ${subdir}/${predictfile} ../aaaa/${dirname}/${subdirname}/${solfile} $outputfile"
#   mv $outputfile ${dirname}/${subdirname}
#}
#
#subrun()
#{
#   subdir=$1
#   dirname=$2
#   subdirname=${1##*/}
#
#   mkdir ${dirname}/${subdirname}
#
#   for i in `seq 0 4`
#   do
#      pred $i $dirname $subdirname $subdir
#   done
#}
#
#run()
#{
#   dir=$1
#   dirname=${1##*/}
#
#   rm -rf ${dirname}
#   mkdir ${dirname}
#
#   for subdir in `ls -d ${dir}/*`
#   do
#      subrun ${subdir} ${dirname}
#   done
#}

run()
{
   buf=${1##*/}
   filename=${buf%.*}

   ../bin/GENE_POINTS $1 5 > ./${filename}_5.lp
   ../bin/GENE_POINTS $1 9 > ./${filename}_9.lp
   ../bin/GENE_POINTS $1 17 > ./${filename}_17.lp
}

datadir=$1

for dir in `ls -d ${datadir}/*.logreg`
do
   run $dir
done

