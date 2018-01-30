#!/bin/sh

run()
{
   name=$1
   cd para_1
   ./run.sh ../../../../Split_data/output/${name} ../../../../exp_split_step/${name} 1
   cd ..

   cd para_2
   ./run.sh ../../../../Split_data/output/${name} ../../../../exp_split_step/${name} 2
   cd ..

   cd para_3
   ./run.sh ../../../../Split_data/output/${name} ../../../../exp_split_step/${name} 3
   cd ..
}

#run breast
#run biodeg
run musk2
