#!/bin/sh

edit()
{
   solfile=$1

   cp $solfile cp.sol
   numline=$(cat cp.sol | wc -l)
   for i in `grep -e objective -n cp.sol | sed -e 's/:.*//g'`
   do
      x=`expr $numline - $i + 1`
      tail -n $x cp.sol > buf.sol
   done
   rm -f cp.sol
   mv buf.sol result/${solfile}
}

soldir=$1
cd $soldir
rm -rf result
mkdir result

for sol in `ls -d *.sol`
do
   edit $sol
done

cd ..
