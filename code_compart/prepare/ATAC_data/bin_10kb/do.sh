#!/bin/bash

chrnum=$1

for i in `seq 1 $chrnum`
do
cp ../chr${i}.dat test.dat
#awk '{print $1,$2,1}' ../chr${i}.dat > test.dat
echo -e test.dat '\n' density${i}.dat | ./getweightdensityfrombed.o
done
