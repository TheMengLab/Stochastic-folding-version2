#!/bin/bash

repid=$1

id1=640
id2=320

#pointnum=`nl ../../../prepare/assignment/multi_chr/allgenome/assign_${id1}kb.dat | tail -n 1 | awk '{print $1}'`
pointnum=`nl ../../../prepare/assign/detail/assign_${id1}kb_rep${repid}.dat | tail -n 1 | awk '{print $1}'`
tail -n ${pointnum} output_${id1}.dat > conf${id1}.dat


#echo -e conf${id1}.dat '\n' ../../../prepare/assignment/multi_chr/allgenome/assign_${id1}kb.dat '\n' ../../../prepare/assignment/multi_chr/allgenome/compartassign_${id1}kb.dat '\n' ../../../prepare/assignment/multi_chr/allgenome/assign_${id2}kb.dat '\n' ../../../prepare/assignment/multi_chr/allgenome/compartassign_${id2}kb.dat '\n' initconf${id2}.dat | ./insert_contract_assign_seg.o 
echo -e conf${id1}.dat '\n' ../../../prepare/assign/detail/assign_${id1}kb_rep${repid}.dat '\n' ../../../prepare/assign/detail/compartassign_${id1}kb.dat '\n' ../../../prepare/assign/detail/assign_${id2}kb_rep${repid}.dat '\n' ../../../prepare/assign/detail/compartassign_${id2}kb.dat '\n' initconf${id2}.dat | ./insert_contract_assign_seg.o 
#Input the filename for the conformation:
#Input filename for first assignment:
#Input the filename for first chr assignment:
#Input filename for second assignment:
#Input the filename for second chr assignment:
#Input the filename for output:
#2

#awk '{print 1}' ../../../prepare/assignment/multi_chr/allgenome/assign_${id2}kb.dat > sudoassign.dat
awk '{print 1}' ../../../prepare/assign/detail/assign_${id2}kb_rep${repid}.dat > sudoassign.dat
cp sudoassign.dat sudodensity.dat
cp sudoassign.dat sudomass.dat

#time echo -e initconf${id2}.dat '\n' sudoassign.dat '\n' sudodensity.dat '\n' 1 '\n' sudomass.dat '\n' ../../../prepare/assignment/multi_chr/allgenome/beadcompartassign_${id2}kb.dat '\n' 1 '\n' 20000 '\n' 20000 '\n' 20000 '\n' ${repid} '\n' output_${id2}.dat | ./runMD_feature_v25_continue_v2_density_mass_compartment_continue.o
time echo -e initconf${id2}.dat '\n' sudoassign.dat '\n' sudodensity.dat '\n' 1 '\n' sudomass.dat '\n' ../../../prepare/assign/detail/beadcompartassign_${id2}kb_rep${repid}.dat '\n' 1 '\n' 20000 '\n' 20000 '\n' 20000 '\n' ${repid} '\n' output_${id2}.dat | ./runMD_feature_v25_continue_v2_density_mass_compartment_continue.o


#pointnum1=`nl ../../../prepare/assign/allgenome_hierarchical/assign${id1}.dat | tail -n 1 | awk '{print $1}'`
#tail -n $pointnum1 output_${id1}.dat > temp.dat
#
#echo -e temp.dat '\n' ../../../prepare/assign/allgenome_hierarchical/assign${id1}.dat '\n' ../../../prepare/assign/allgenome_hierarchical/assign${id2}.dat '\n' conf${id2}.dat | ./insert_contract_assign.o
#
#
#
#time echo -e conf${id2}.dat '\n' ../../../prepare/assign/allgenome_hierarchical/assign${id2}.dat '\n' ../../../prepare/assign/DNase/hierarchical/density${id2}.dat '\n' 0.015 '\n' 1 '\n' 5000 '\n' 5000 '\n' 5000 '\n' ${i} '\n' output_${id2}.dat | ./runMD_feature_v25_continue_v2_density.o
##time echo -e conf${id2}.dat '\n' ../../../data/medcutoff/chr_all/assignall_${id2}kb.dat '\n' 1 '\n' 10000 '\n' 10000 '\n' 1 '\n' output_${id2}.dat | ./runMD_feature_v25_continue.o
##Input filename for original data:
##Input the assignment for chromatin point:
##Input the value for bond lenth:
##Input the step number for iteration:
##Input the interval for saving:
##Input the value for random number generator:
##Input filename for output:
#
