#!/bin/bash

chrid=$1
delta=$2
winsize=$3
#cp /media/group/2b2676b5-76e1-42af-864b-2152cf8da195/20201210/code/c/Nucleosome_model/segment_sample/HiC/simulation/Hi-C_highresolution/K562/10kb_resolution_intrachromosomal/compartment/compartment_assign/clustering_H3K4me3_assign/chr_10kb/assign${chrid}_10kb.dat .
cp ../../../../prepare/assign/detail/assign${chrid}_10kb_compbk.dat assign${chrid}_10kb.dat

statenum=`nl ../../matrix4_chr${chrid}_norm_top100.dat | tail -n 1 | awk '{print $1}'`
#echo -e ../test.dat '\n' $statenum '\n' test.dat '\n' score.dat | ./getarrowhead_maxscore_nodark.o
#echo -e ../test.dat '\n' ${statenum} '\n' assign${chrid}_10kb.dat '\n' test.dat '\n' score.dat | ./getarrowhead_maxscore_nodark_darkassign.o 
echo -e ../../matrix4_chr${chrid}_norm_top100.dat '\n' ${statenum} '\n' assign${chrid}_10kb.dat '\n' 100 '\n' test_topline.dat '\n' score_topline.dat | ./getarrowhead_maxscore_nodark_darkassign_topline.o
#Input filename for original data:
#Input the dimension of matrix:
#Input the filename for compartassign(1 col):
#Input filename for arrowhead output:
#Input the filename for score list:


#echo -e ../test.dat '\n' 5136 '\n' test.dat | ./getarrowhead.o 
##Input filename for original data:
##Input the dimension of matrix:
##Input filename for output:

#echo -e score.dat '\n' score_statis.dat | ./getaverage.o
echo -e score_topline.dat '\n' assign${chrid}_10kb.dat '\n' score_statis.dat | ./getaverage_darkassign.o

#cutoff=`cat score_statis.dat | awk '{print $1}'`
cutoff=`cat score_statis.dat | awk '{print $2*0.25}'`
cat score_statis.dat
cutoff=${delta}
echo ${cutoff}

#echo -e score_topline.dat '\n' 20 '\n' peakid.dat '\n' test.dat | ./gethighpeak.o
echo -e score_topline.dat '\n' ${winsize} '\n' peakid_delta.dat '\n' test.dat | ./gethighpeak_delta.o
#awk -v CUTOFF=${cutoff} '{if(($2>0)&&($2-$3>CUTOFF)) print $1,$2}' peakid_delta.dat > targetpeakid.dat
awk -v CUTOFF=${cutoff} '{if($2>CUTOFF) print $1,$2}' peakid_delta.dat > targetpeakid.dat


##echo -e targetpeakid.dat '\n' /media/group/2b2676b5-76e1-42af-864b-2152cf8da195/20201210/code/c/Nucleosome_model/segment_sample/HiC/simulation/Hi-C_highresolution/K562/10kb_resolution_intrachromosomal/boundary/delta/${delta}/win${winsize}/${chrid}/targetpeakid.dat '\n' dist.dat | ./getdomain_dist.o
#echo -e targetpeakid.dat '\n' /media/group/2b2676b5-76e1-42af-864b-2152cf8da195/20201210/code/c/Nucleosome_model/segment_sample/HiC/simulation/Hi-C_highresolution/K562/10kb_resolution_intrachromosomal/boundary/delta/${delta}/win${winsize}_nodelta/${chrid}/targetpeakid.dat '\n' dist.dat | ./getdomain_dist.o
#
#
#
#awk '{print $2}' dist.dat | sort -g > dist_sort.dat

rm test_topline.dat

