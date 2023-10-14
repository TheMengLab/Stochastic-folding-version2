#!/bin/bash


id=1

for repid in `seq 1 600`
do

statenum=`nl ../../prepare/assign/assign_10kb_rep${repid}.dat | tail -n 1 | awk '{print $1}'`
#statenum=`nl ../../prepare/assign/detail/assign_10kb_rep${repid}.dat | tail -n 1 | awk '{print $1}'`

tail -n ${statenum} ../../simulation/replica/${repid}/output_10.dat ${repid}.dat
cd mergepoint
echo -e ../${repid}.dat '\n' ../../../prepare/assign/detail/assign_10kb_rep${repid}.dat '\n' ${repid}.dat | ./mergecoor.o
cd ..
rm ${repid}.dat



done

