#!/bin/bash

for i in `seq 1 23`
do
#awk -v CID=${i} '{print CID-1,$1}' /media/group/2b2676b5-76e1-42af-864b-2152cf8da195/20201210/code/c/Nucleosome_model/segment_sample/HiC/simulation/single_cell/Steven/data/Population/compartment/compartment_assign/PCA_H3K4me3_assign/chr_10kb/assign${i}_10kb.dat > compartassign${i}_10kb.dat
cp /media/group/2b2676b5-76e1-42af-864b-2152cf8da195/20201210/code/c/Nucleosome_model/segment_sample/HiC/simulation/Hi-C_highresolution/K562/10kb_resolution_intrachromosomal/compartment/compartment_assign/clustering_H3K4me3_assign/chr_10kb/assign${i}_10kb.dat assign${i}_10kb_compbk.dat
#awk -v CID=${i} '{print$1}' /media/group/2b2676b5-76e1-42af-864b-2152cf8da195/20201210/code/c/Nucleosome_model/segment_sample/HiC/simulation/single_cell/Steven/data/Population/compartment/compartment_assign/PCA_H3K4me3_assign/chr_10kb/assign${i}_10kb.dat > compartassign${i}_10kb.dat
done
