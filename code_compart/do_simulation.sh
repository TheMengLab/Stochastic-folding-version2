#!/bin/bash


for i in `seq 1 10`
do
cd simulation/replica/${i}/
for j in `seq 1 50`
do
cd simulation/replica/${i}/simulation_nocompart/${j}/
bash doall.sh
cd ../../
done
cd ../../..
done
