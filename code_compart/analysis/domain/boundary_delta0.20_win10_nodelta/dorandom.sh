#!/bin/bash

delta=$1
winsize=$2

#for i in `seq 1 23`
#do
#cd ${i}
#bash do.sh $i $delta $winsize
#cd ../
#done

cat */random/dist_sort*.dat | sort -g | grep -v 10000 > test_random.dat
statenum=`nl test_random.dat | tail -n 1 | awk '{print $1}'`
targetcount=`cat test_random.dat | awk '{if($1<10.5) count++} END {print count}'`
totalcount=`cat /media/group/2b2676b5-76e1-42af-864b-2152cf8da195/20201210/code/c/Nucleosome_model/segment_sample/HiC/simulation/Hi-C_highresolution/K562/10kb_resolution_intrachromosomal/boundary/delta/${delta}/win${winsize}_nodelta/*/targetpeakid.dat | nl | tail -n 1 | awk '{print $1}'`
echo $targetcount $statenum $totalcount | awk '{print $1,$2,$3}' > statis_random.dat


cat statis_random.dat | awk '{print $1/$2}' > statis_random_ratio.dat
