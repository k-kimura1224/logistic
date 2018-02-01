#!/bin/sh

run3()
{
   data=$1
   method=$2
   para=$3
   for file in `ls -d ./${method}/${para}/${data}/*/*.pred`
   do
      sens=$(sed -n 10p $file)
      echo "${method} $sens" >> ${data}_${para}.bp
      spec=$(sed -n 12p $file)
      echo "${method} $spec" >> ${data}_${para}.bp
   done
}

run2()
{
   data=$1
   method=$2
   run3 $data $method para_1
   run3 $data $method para_2
   run3 $data $method para_3
}

run1()
{
   data=$1
   run2 $data step_b
   run2 $data step_f
   run2 $data ug
}

rm -f *.bp

run1 breast
run1 biodeg
#run1 musk2
run1 parkin
run1 seismic-bumps
run1 spectf
#run1 statG
#run1 statH
#run1 madelon_train
