#!/bin/bash
chrnum=$1

for i in `seq 1 $chrnum`
do
nohup bash temp.sh ${i} &
done




