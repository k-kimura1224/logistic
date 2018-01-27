#!/bin/sh

pred()
{
   i=$1
   dirname=$2
   subdirname=$3
   subdir=$4
   soldir=$5
   samplefile=${i}_sample_${dirname}.logreg
   predictfile=${i}_predict_${dirname}.logreg
   solfile=${i}_sample_${dirname}.sol
   outputfile=${i}_${dirname}.pred
   #buf=${1##*/}
   #filename=${buf%.*}

   #echo "../bin/prediction ${subdir}/${samplefile} ${subdir}/${predictfile} ${soldir}/${subdirname}/edit_sol/${solfile} $outputfile"
   ../bin/prediction ${subdir}/${samplefile} ${subdir}/${predictfile} ${soldir}/${subdirname}/edit_sol/${solfile} $outputfile
   mv $outputfile ${dirname}/${subdirname}
}

subrun()
{
   subdir=$1
   dirname=$2
   subdirname=${1##*/}
   soldir=$3

   mkdir ${dirname}/${subdirname}

   for i in `seq 0 4`
   do
      pred $i $dirname $subdirname $subdir $soldir
   done
}

run()
{
   dir=$1
   dirname=${1##*/}
   soldir=$2

   rm -rf ${dirname}
   mkdir ${dirname}

   for subdir in `ls -d ${dir}/*`
   do
      subrun ${subdir} ${dirname} ${soldir}
   done
}

datadir=$1
soldir=$2

run $datadir $soldir

