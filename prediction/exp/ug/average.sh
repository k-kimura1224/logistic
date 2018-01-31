#!/bin/sh

run()
{
   dir=$1
   echo "** dir = $1 **********************"

   ct=0
   sum_sens=0
   sum_spec=0
   for file in `ls -d ${dir}/*/*.pred`
   do
      sens=$(sed -n 10p $file)
      sum_sens=$(echo "$sum_sens + $sens" | bc -l)

      spec=$(sed -n 12p $file)
      sum_spec=$(echo "$sum_spec + $spec" | bc -l)

      ct=`expr $ct + 1`
   done

   ave_sens=$(echo "$sum_sens / $ct" | bc -l)
   ave_spec=$(echo "$sum_spec / $ct" | bc -l)

   echo "Sensitivity: $ave_sens"
   echo "Specificity: $ave_spec"
}

print_result()
{
   name=$1

   echo "-- $name ---------"
   run para_1/$name
   run para_2/$name
   run para_3/$name
   echo ""
}

print_result breast
#print_result biodeg
#print_result musk2
#print_result parkin
#print_result seismic-bumps
#print_result spectf
#print_result statG
#print_result statH
print_result madelon_train
