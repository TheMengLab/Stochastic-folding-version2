#!/bin/bash

repid=$1
confid=${repid}
res=5120

#cp ../../../data/medcutoff/chr_all/assignall_1280kb.dat .

#awk '{print 1}' ../../../prepare/assignment/multi_chr/allgenome/assign_${res}kb.dat > sudoassign.dat
awk '{print 1}' ../../../prepare/assign/detail/assign_${res}kb_rep${repid}.dat > sudoassign.dat

cp sudoassign.dat sudodensity.dat
cp sudoassign.dat sudomass.dat
for i in `seq $confid $confid`
do
#cp ../../initconf/conf${i}.dat  initconf${res}.dat
cp ../../initconf/conf${i}.dat  initconf${res}.dat
time echo -e initconf${res}.dat '\n' sudoassign.dat '\n' sudodensity.dat '\n' 1 '\n' sudomass.dat '\n' ../../../prepare/assign/detail/beadcompartassign_${res}kb_rep${repid}.dat '\n' 1 '\n' 100000 '\n' 100000 '\n' 100000 '\n' ${confid} '\n' output_${res}.dat | ./runMD_feature_v25_continue_v2_density_mass_compartment.o
#Input filename for original data:
#Input the assignment for chromatin point:
#Input the value for bond lenth:
#Input the step number for iteration:
#Input the interval for saving:
#Input the value for random number generator:
#Input filename for output:

done
