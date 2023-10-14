#!/bin/bash


#echo -e density22.dat '\n' 1 '\n' 1 '\n' 1 '\n' 10 '\n' testmatrix_rep10.dat '\n' pointlen_rep10.dat | ./getrandomcontact_norm_pointlen.o
#basevalue=`cat /media/group/2b2676b5-76e1-42af-864b-2152cf8da195/20201210/code/c/Nucleosome_model/segment_sample/HiC/simulation/single_cell/MD_attraction/compartment_individualcell/data/sc-ATAC-seq/K562/test_average/ave_median.dat`
#basevalue=`cat /media/group/2b2676b5-76e1-42af-864b-2152cf8da195/20201210/code/c/Nucleosome_model/segment_sample/HiC/simulation/single_cell/MD_attraction/compartment_individualcell/data/sc-ATAC-seq/K562/test_average/ave_statis.dat`
basevalue=`cat ../data_raw/ATAC_sc/Rep*/*/densityall.dat | nl | awk '{SUM+=$2} END {print SUM/$1}'`

for i in `seq 1 6`
do
for j in `seq 1 100`
do
cp ../data_raw/ATAC_sc/Rep${i}/${j}/densityall.dat .

echo -e densityall.dat '\n' ${basevalue} '\n'  len$[i*100+j-100].dat | ./getnorm_pointlen.o

done
done


