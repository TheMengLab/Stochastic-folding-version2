#!/bin/bash

repid=$1

#cp ../../../data/medcutoff/chr_all/assignall_1280kb.dat .

for i in `seq $repid $repid`
do
cp ../../initconf/conf${i}.dat  initconf.dat
#time echo -e initconf.dat '\n' ../../assignment/sudoassign.dat '\n' ../../assignment/sudodensity.dat '\n' 1 '\n' ../../assignment/sudomass.dat '\n' ../../assignment/compartmentasssign.dat '\n' 1 '\n' 1000000 '\n' 1000000 '\n' 5000 '\n' ${i} '\n' output_1.dat | ./runMD_feature_v25_continue_v2_density_mass_compartment.o
time echo -e initconf.dat '\n' ../../assignment/sudoassign.dat '\n' ../../assignment/sudodensity.dat '\n' 1 '\n' ../../assignment/sudomass.dat '\n' 1 '\n' 1000000 '\n' 1000000 '\n' 5000 '\n' ${i} '\n' output_1.dat | ./runMD_feature_v25_continue_v2_density_mass.o
#Input filename for original data:
#Input the assignment for chromatin point:
#Input the value for bond lenth:
#Input the step number for iteration:
#Input the interval for saving:
#Input the value for random number generator:
#Input filename for output:

done
