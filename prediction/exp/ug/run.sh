#!/bin/sh

run()
{
   name=$1

   cd para_1
   ./run.sh ../../../../Split_data/output/${name} ../../../../exp_split_ug/${name} 1
   cd ..

   cd para_2
   ./run.sh ../../../../Split_data/output/${name} ../../../../exp_split_ug/${name} 2
   cd ..

   cd para_3
   ./run.sh ../../../../Split_data/output/${name} ../../../../exp_split_ug/${name} 3
   cd ..
}

#run biodeg
#run breast
#run musk2
#run madelon_train
#run parkin
#run seismic-bumps
run spectf
