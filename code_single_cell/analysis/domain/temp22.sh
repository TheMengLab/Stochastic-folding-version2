#!/bin/bash

chrid=22

filenum=`ls ../[1-9]*.dat | nl | tail -n 1 | awk '{print $1}'`

for id in `seq ${chrid} ${chrid}`
do

nl ../../../assign/detail/chrassign_10kb.dat | awk '{print $1,$2}' | grep " "$[id-1]$ > temp_chr${id}.dat
#nl ../../../../../assignment/multi_chr/allgenome/chrassign_10kb.dat | awk '{print $1,$2}' | grep " "$[id-1]$ > temp_chr${id}.dat
sid=`head -n 1 temp_chr${id}.dat | awk '{print $1}'`
eid=`tail -n 1 temp_chr${id}.dat | awk '{print $1}'`
echo $id $sid $eid

for i in `seq 1 $filenum`
do
head -n $[eid] ../${i}.dat | tail -n $[eid-sid+1] > tempcoor/${i}_chr${id}.dat
done

ls tempcoor/[1-9]*_chr${id}.dat > filelist_chr${id}.dat

nl tempcoor/1_chr${id}.dat | awk '{print $1,1}' > targetid_chr${id}.dat

atomnum=`nl tempcoor/1_chr${id}.dat | tail -n 1 | awk '{print $1}'`

time echo -e filelist_chr${id}.dat '\n' ${atomnum} '\n' targetid_chr${id}.dat '\n' 4 '\n' 100 '\n' matrix4_chr${id}_norm_top100.dat | ./getcontactmatrix_dark_noneigh_norm_topline.o 
#Input the filename for filelist of original coordinate:
#Input the number of atoms in this system:
#Input filename for target atomlist(start from 1):
#Input the cutoff for contact dist(with radius of 10):
#Input the value for top linenum:
#Input filename for output:


#time echo -e filelist.dat '\n' $[eid-sid+1] '\n' targetid.dat '\n' 4 '\n' matrix4_norm.dat | ./getcontactmatrix_dark_noneigh_norm.o
#time echo -e filelist.dat '\n' $[eid-sid+1] '\n' targetid.dat '\n' 4 '\n' matrix4_norm_chr${id}.dat '\n' 1 | ./getcontactmatrix_dark_noneigh_norm_targetchrid_coarse_scale.o

#time echo -e filelist.dat '\n' $[eid-sid+1] '\n' 1 '\n' matrix${id}_10kb.dat | ./getdistmatrix_coor_normlog_median_ratio.o
#
#
##cp run.r_bk run.r
##sed -i "s/matrix10/matrix4_norm_chr${id}/g" run.r
##Rscript run.r
##mv domain_50.dat domain_50_matrix4_norm_chr${id}.dat
##rm matrix4_norm_chr${id}.dat
#rm tempcoor/[1-9]*_chr${id}.dat
done


