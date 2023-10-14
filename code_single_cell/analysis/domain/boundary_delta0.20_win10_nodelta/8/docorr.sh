#!/bin/bash

chrid=$1
delta=$2
winsize=$3
cp /media/group/2b2676b5-76e1-42af-864b-2152cf8da195/20201210/code/c/Nucleosome_model/segment_sample/HiC/simulation/Hi-C_highresolution/K562/10kb_resolution_intrachromosomal/compartment/compartment_assign/clustering_H3K4me3_assign/chr_10kb/assign${chrid}_10kb.dat .


echo -e score_topline.dat '\n' /media/group/2b2676b5-76e1-42af-864b-2152cf8da195/20201210/code/c/Nucleosome_model/segment_sample/HiC/simulation/Hi-C_highresolution/K562/10kb_resolution_intrachromosomal/boundary/delta/${delta}/win${winsize}_nodelta/${chrid}/score_topline.dat '\n' assign${chrid}_10kb.dat '\n' 500 '\n' corr${chrid}.dat | ./getcorr_darkassign.o




