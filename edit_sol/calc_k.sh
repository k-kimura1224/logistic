#!/bin/sh

calc()
{
   solfile=$1

   echo "$solfile ...\c"
   grep "z" $solfile | wc -l
}

soldir=$1

for sol in `ls -d $soldir/*.sol`
do
   calc $sol
done

