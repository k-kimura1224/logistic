#/bin/sh
#chmod u+x gen.sh
generate(){
   file=$1
   buf=${1##*/}
   name=${buf%.*}
   txtname=$2

   echo "read $file" >> $txtname
   echo "read ${name}_init.sol" >> $txtname
   echo "opt" >> $txtname
   echo "write ${name}.sol" >> $txtname
   echo "" >> $txtname
}

time=$2
txtname=run.txt
dirname=$1

rm -f $txtname

echo "generate $txtname"
echo "set timelimit $time" >> $txtname
echo "set mip tolerances mipgap 1e-10" >> $txtname
echo "set threads 16" >> $txtname
echo "" >> $txtname

for file in `ls ${dirname}/*.lp`
do
   generate $file $txtname
done

echo "q" >> $txtname
