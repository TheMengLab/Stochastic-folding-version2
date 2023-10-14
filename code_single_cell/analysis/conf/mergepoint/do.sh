#!/bin/bash

repid=$1
cp ../../../../../assign/detail/assign_10kb_rep${repid}.dat assign.dat

rm [1-9]*.dat
confnum=`ls ../[0-9]*.dat | nl | tail -n 1 | awk '{print $1}'`


for i in `seq 1 $confnum`
do
echo -e ../${i}.dat '\n' assign.dat '\n' ${i}.dat | ./mergecoor.o
done
