#!/bin/bash

py=/bin/python

start=1
end=15
cores=8
for runid in $(seq $start $end)
do
	mpirun -np $cores $IMP python ./sin3a_hdac1_modeling.py prod $runid &
done

wait
sleep 20s
echo May the Force be with you...
