#!/bin/bash

#chrid=1
repnum=600


for i in `seq 1 ${chrnum}`
do
awk '{print 1}' ../../chrlen/tempassign_10kb/assign${i}.dat > assign${i}_10kb_compbk.dat
done



for i in `seq 1 ${repnum}`
do


for j in `seq 1 $chrnum`
do
cat assign${j}_10kb_compbk.dat
done > assign_10kb_compbk.dat

cp ../len${i}.dat len_10kb_rep${i}.dat

for j in `seq 1 23`
do
awk -v CID=${j} '{print CID-1}' assign${j}_10kb_compbk.dat
done > chrassign_10kb.dat

paste chrassign_10kb.dat assign_10kb_compbk.dat | awk '{print $1,$2}' > compartassign_10kb.dat



list="
10
20
40
80
160
320
640
1280
2560
"

for res in $list
do
echo -e len_${res}kb_rep${i}.dat '\n' compartassign_${res}kb.dat '\n' 2 '\n' len_$[res*2]kb_rep${i}.dat '\n' assign_$[res*2]kb_rep${i}.dat '\n' compartassign_$[res*2]kb.dat | ./gendensity_detailbead_mass_coarse_chrassign.o
#echo -e pointlen_${res}kb.dat '\n' 2 '\n' pointlen_$[res*2]kb.dat '\n' assign_$[res*2]kb.dat | ./gendensity_detailbead_mass_coarse.o 
#Input the filename for original pointlen:
#Input the ratio for coarse-graining:
#Input the filename for output pointlen:
#Input the filename for output assignment:

done


res=10
echo -e len_${res}kb_rep${i}.dat '\n' compartassign_${res}kb.dat '\n' 1 '\n' len_$[res]kb_rep${i}_check.dat '\n' assign_$[res]kb_rep${i}.dat '\n' compartassign_$[res]kb_bk.dat | ./gendensity_detailbead_mass_coarse_chrassign.o

for res in $list
do
echo -e assign_${res}kb_rep${i}.dat '\n' compartassign_${res}kb.dat '\n' beadcompartassign_${res}kb_rep${i}.dat | ./micro2macromap_chain.o
done


res=5120

 echo -e assign_${res}kb_rep${i}.dat '\n' compartassign_${res}kb.dat '\n' beadcompartassign_${res}kb_rep${i}.dat | ./micro2macromap_chain.o





done









#
#for i in `seq 1 23`
#do
#cat ../density${i}.dat
#done > densityall.dat
#
#for i in `seq 1 23`
#do
#cat assign${i}_10kb_compbk.dat
#done > assign_10kb_compbk.dat
#
#
#for i in `seq 1 23`
#do
#awk -v CID=${i} '{print CID-1}' ../density${i}.dat
#done > chrassign_10kb.dat
#
#
#paste chrassign_10kb.dat assign_10kb_compbk.dat | awk '{print $1,$2}' > compartassign_10kb.dat
##cp assign_10kb_compbk.dat compartassign_10kb.dat
#
#sort -g densityall.dat | grep -v 0$ > test.dat
#linenum=`nl test.dat | tail -n 1 | awk '{print $1}'`
#basevalue=`sort -g test.dat | head -n $[linenum/2] | tail -n 1 | awk '{print $1*3/4}'`
#
#echo $basevalue
#echo -e densityall.dat '\n' ${basevalue} '\n' 8 '\n' assign_10kb.dat '\n' pointlen_10kb.dat '\n' mass_10kb.dat | ./gendensity_detailbead_mass.o
##echo $basevalue
##echo -e density_chr1_20M_50M.dat '\n' 0.014 '\n' 6 '\n' assign.dat '\n' pointlen.dat| ./gendensity_detailbead.o 
##echo -e densityall.dat '\n'  8 '\n' assign.dat '\n' pointlen.dat '\n' mass.dat | ./gendensity_detailbead_mass_uplimit.o 
##echo -e density_chr1_20M_50M.dat '\n' 0.014 '\n' 6 '\n' assign.dat '\n' pointlen.dat '\n' mass3.dat | ./gendensity_detailbead_mass3.o 
##Input filename for density data:
##Input filename for interval:
##Input the maxlen for each bin:
##Input filename for output:
#
#list="
#10
#20
#40
#80
#160
#320
#640
#1280
#2560
#"
#
#for res in $list
#do
#echo -e pointlen_${res}kb.dat '\n' compartassign_${res}kb.dat '\n' 2 '\n' pointlen_$[res*2]kb.dat '\n' assign_$[res*2]kb.dat '\n' compartassign_$[res*2]kb.dat | ./gendensity_detailbead_mass_coarse_chrassign.o
##echo -e pointlen_${res}kb.dat '\n' 2 '\n' pointlen_$[res*2]kb.dat '\n' assign_$[res*2]kb.dat | ./gendensity_detailbead_mass_coarse.o 
##Input the filename for original pointlen:
##Input the ratio for coarse-graining:
##Input the filename for output pointlen:
##Input the filename for output assignment:
#
#done
#
#
##for i in `seq 1 20`
##do
##cat assign${i}_10kb_compbk.dat
##done > compartassign_10kb.dat
#
#
#for res in $list
#do
#
## for chrid in `seq 1 20`
## do
## awk -v CID=${chrid} '{print CID-1,$1}' compartassign${chrid}_${res}kb.dat 
## done > temp.dat
##
## mv temp.dat compartassign_${res}kb.dat
# 
# 
# #awk '{print $2}' compartassign1_${res}kb.dat > compartassign1_${res}kb_1col.dat
# 
##awk '{print $2}' compartassign_$[res]kb.dat > temp.dat
# 
# echo -e assign_${res}kb.dat '\n' compartassign_${res}kb.dat '\n' beadcompartassign_${res}kb.dat | ./micro2macromap_chain.o
#
#
#done
#
#
#res=5120
#
# echo -e assign_${res}kb.dat '\n' compartassign_${res}kb.dat '\n' beadcompartassign_${res}kb.dat | ./micro2macromap_chain.o
#
##
##res=640
##
##awk -v CID=${chrid} '{print CID-1,$1}' compartassign${chrid}_${res}kb.dat > temp.dat
##mv temp.dat compartassign${chrid}_${res}kb.dat
##
##
###awk '{print $2}' compartassign1_${res}kb.dat > compartassign1_${res}kb_1col.dat
##
##
##echo -e assign${chrid}_${res}kb.dat '\n' compartassign${chrid}_${res}kb.dat '\n' beadcompartassign${chrid}_${res}kb.dat | ./micro2macromap_chain.o
##
##
