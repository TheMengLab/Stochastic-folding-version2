#!/bin/bash

chrnum=23	#The number of chromosome in the simulated cell, the value of 23 means 22 autosomes+ X chromosome



#Starting from ATAC/DNase data(genomic length of each chromosome is also needed) to prepare input files for performing simulation for stochastic folding
cd prepare/


#convert the information of chromosome size to the number of 10-kb segments. The information is deposited in a serise of file ./chrlen/tempassign_10kb/assign1.dat, ./chrlen/tempassign_10kb/assign2.dat, and etc. For example, the line number of ./chrlen/tempassign_10kb/assign1.dat corresponds the number of 10-kb segments in chromosome 1
cd chrlen/
bash doall.sh $chrnum
cd ..


#convert the DNA accessibility data to the signal for each 10-kb segments. These 10-kb segments signals are deposited in ./ATAC_data/bin_10kb/further/density1.dat, ./ATAC_data/bin_10kb/further/density2.dat, and etc. For example, the value of the first row in ./ATAC_data/bin_10kb/further/density2.dat represents the DNA accessibility signal value of the first 10-kb segments(0-10kb) of chromosome 2
cd ATAC_data/
bash doall.sh $chrnum
cd ..


#determin the number of beads (or k values in the manuscript) for each 10-kb segments. The outputs of the command are files ./assign/detail/len_10kb_rep1.dat, ./assign/detail/len_10kb_rep2.dat and etc(10 files in total). 10 file corresponds to 10 polymer models. The difference between two polymer models lies in the varied k values (algorithm of determining k values are described in manuscript). In each file, the k value of i-th 10-kb segment are listed as the i-th row.

cd assign/
bash doall.sh $chrnum	#chrnum is the number of chromosome in the reference genome. Since we only involve 22 autosomes and one X chromosome in this study, thus the we take 23 as the value of chrnum
cd ..

cd ../



#generating initial structures of low resolution level(5120kb resolution in this work) for simulation. Also, modify some parameter the scripts for the random folding conformation of 10-kb resolution level.
cd simulation/replica/
bash prepare_simulation.sh $chrnum
cd ../..

