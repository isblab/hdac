#!/bin/bash

py=/bin/python

cores=8
for runid in 1 2 3 4 5 6 7 8 9 10
do
	mpirun -np $cores $IMP python ./hdac1_modeling.py prod $runid &
done

wait
sleep 2m
echo May the Force be with you...
