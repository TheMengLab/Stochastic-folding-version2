#!/bin/bash


for i in `seq 1 600`
do
cd simulation/replica/${i}
bash doall.sh
cd ../../..
done

