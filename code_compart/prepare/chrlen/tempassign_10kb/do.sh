#!/bin/bash
chrnum=$1

for i in `seq 1 ${chrnum}`
do
chrlen=`awk '{printf("%d\n",($1-1)/10000+1)}' ../chrlen${i}.dat`

echo $chrlen
for j in `seq 1 $chrlen`
do
echo 1
done > assign${i}.dat


done
