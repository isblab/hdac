#!/bin/bash

py=/bin/python
start=$1
end=$2

#/home/kartik/imp-clean/build/setup_environment.sh 
for i in $(seq $start $end)
do	
    cd Output_$i/ 
    mkdir Output_$i/
    cp -r rmfs/ Output_$i
    echo Done for $i
    cd ../
done
