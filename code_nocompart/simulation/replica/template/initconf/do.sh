#!/bin/bash

repid=$1
startid=1
chrnum=$2
num=50

nl ../../../../prepare/assign/detail/beadcompartassign_5120kb_rep${repid}.dat | awk '{print $1-1,$2,$3}' > assign.dat
#nl ../../../assign/detail/beadcompartassign_5120kb_rep${repid}.dat | awk '{print $1-1,$2,$3}' > assign.dat

rm conf*.dat

for i in `seq 1 ${num}`
do

echo -e ${chrnum} '\n' $[i+repid*num] '\n' testconf.dat | ./genrandomconf_nounit.o 
#Input the number of beads in the system:
#Input the initial seed for random number generation:
#Input the filename for output:



echo -e testconf.dat '\n' assign.dat '\n' $[i+repid*num] '\n' conf${i}.dat | timeout 10  ./genrandomconf_nounit_startpoint.o


#check bad conf
count=`ls conf${i}.dat | nl | tail -n 1 | awk '{print $1+0}'`
if [ -f conf${i}.dat ]
then
count=1
else
count=0
fi
echo $count

tempid=${i}
#while [ -f conf${i}.dat ]
while [ $count -ne 1 ]
do
echo check conf $i
tempid=$[tempid+num]

echo -e ${chrnum} '\n' $[tempid+3*repid*num] '\n' testconf.dat | ./genrandomconf_nounit.o 
#Input the number of beads in the system:
#Input the initial seed for random number generation:
#Input the filename for output:

echo -e testconf.dat '\n' assign.dat '\n' ${i} '\n' conf${i}.dat | timeout 10  ./genrandomconf_nounit_startpoint.o

#echo -e testconf.dat '\n'  chrassign.dat '\n' chrcontact.dat '\n' 10000 '\n' 10000 '\n' ${tempid} '\n' output${i}.dat | ./runMD_Morse_list_block_all_noghost_prob_list_cell_v34_continue.o
#tail -n ${chrnum} output${i}.dat | awk '{print $1*40,$2*40,$3*40}'> startcoor${i}.dat
#echo -e startcoor${i}.dat '\n' ../../../prepare/assignall${id}.dat '\n' ${i} '\n' conf${i}.dat | timeout 10 ./genrandomconf_nounit_startpoint.o

if [ -f conf${i}.dat ]
then
count=1
else
count=0
fi

#count=`ls conf${i}.dat | nl | tail -n 1 | awk '{print $1+0}'`

done 
done

