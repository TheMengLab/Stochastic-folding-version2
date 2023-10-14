#!/bin/bash


id=1

for repid in `seq 1 10`
do

statenum=`nl ../../prepare/assign/detail/assign_10kb_rep${repid}.dat | tail -n 1 | awk '{print $1}'`

for j in `seq 1 50`
do
ls ../../simulation/replica/${repid}/simulation_nocompart/${j}/output_10.dat
done > filelist.dat




while read line
do
tail -n ${statenum} ${line} > ${id}.dat

cd mergepoint
echo -e ../${id}.dat '\n' ../../../prepare/assign/detail/assign_10kb_rep${repid}.dat '\n' ${id}.dat | ./mergecoor.o
cd ..
rm ${id}.dat

id=$[id+1]
done < filelist.dat



done

