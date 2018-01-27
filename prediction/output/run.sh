#!/bin/sh

pred()
{
   i=$1
   dirname=$2
   subdirname=$3
   subdir=$4
   samplefile=${i}_sample_${dirname}.logreg
   predictfile=${i}_predict_${dirname}.logreg
   solfile=${i}_sample_${dirname}.sol
   outputfile=${i}_${dirname}.pred
   #buf=${1##*/}
   #filename=${buf%.*}

   #echo "../bin/prediction ${subdir}/${samplefile} ${subdir}/${predictfile} ../../exp_split_ug/${dirname}/${dirname}/${subdirname}/edit_sol/${solfile} $outputfile"
   ../bin/prediction ${subdir}/${samplefile} ${subdir}/${predictfile} ../../exp_split_ug/${dirname}/${subdirname}/edit_sol/${solfile} $outputfile
   mv $outputfile ${dirname}/${subdirname}
}

subrun()
{
   subdir=$1
   dirname=$2
   subdirname=${1##*/}

   mkdir ${dirname}/${subdirname}

   for i in `seq 0 4`
   do
      pred $i $dirname $subdirname $subdir
   done
}

run()
{
   dir=$1
   dirname=${1##*/}

   rm -rf ${dirname}
   mkdir ${dirname}

   for subdir in `ls -d ${dir}/*`
   do
      subrun ${subdir} ${dirname}
   done
}

datadir=$1

run $datadir

