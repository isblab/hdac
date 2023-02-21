#!/bin/bash

py=/bin/python

start=1
end=30
cores=4
for runid in $(seq $start $end)
do
	mpirun -np $cores $IMP python ./coREST_hdac1_modeling.py prod $runid &
done

wait
sleep 2m
echo May the Force be with you...
