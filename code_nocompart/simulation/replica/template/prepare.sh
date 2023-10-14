#!/bin/bash
repid=$1
chrnum=$2

cd initconf/
bash do.sh $repid $chrnum
cd ..

cd simulation_nocompart
bash copy.sh $repid
cd ..


