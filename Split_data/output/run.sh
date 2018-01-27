#!/bin/sh
#chmod u+x gen.sh
generate(){
   file=$1
   buf=${1##*/}
   name=${buf%.*}
   k=$2
   r=$3
   samplename=sample_${buf}
   predictname=predict_${buf}
   rm -rf $name
   mkdir $name
	for i in `seq 0 $r`
   do
      mkdir $name/$i
      ../bin/scip $file $k $samplename $predictname
      mv *.logreg $name/$i
   done
}

k=5
r=9

#for file in `ls ../../data/*.logreg`
#do
#   generate $file $k $r
#done
generate ../../data/biodeg.logreg $k $r


