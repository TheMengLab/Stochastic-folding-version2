Overview

1.Introduction
2.Repo Contents
3.System Requirements
4.Installation Guide
5.Demo
6.Results



1.Introduction
The repository contains the codes and script for simulating 3D chromosome structures from single-cell DNA accessibility data based on stochastic folding model. The code is written by C and the script is written in bash
Both the source code and the compiled code are presented in the repository.
The input data for the codes contains two parts: (1)the single-cell DNA accessibility data(we use ATAC-seq of K562 as a demonstration here, deposited in the folder ./prepare/data_raw/ATAC_sc/ and (2) the sequence length of each chromosome in the reference genome (we use hg19 genome, deposited as ./prepare/data_raw/hg19.chrom.sizes)
The output of the codes contains the simulated conformation ensemble, together with the futher analyzed results, such as the separation score of each chromosome segment, the genomic positions of predicted domain boundaries, and etc.





2.Repo Contents
prepare: repository to convert DNA accessibility data to the input data prepared for simulating 3D chromatin structures 
simulation: repository for calculating folded chromatin structures in confined space based on stochastic folding.
analysis: repository to identify domain boundary from simulated chromatin structure ensemble.

In the respository ./prepare/data_raw/, there are two types of data obtained from experimental: 
(1)single-cell ATAC data

(2)hg19.chrom.sizes: file containing the length of reference genome.


The calculation can be roughly divided into 3 steps:
(1) preparation step. In this step, we construct 10 polymer models from DNA accessibility data and generate low-resolution initial structures as starting points for structural simulation
command: bash prepareall.sh
The scripte prepareall.sh contains multiple commands and each command has its own functions. Detail descriptions for commands in the script are documented in the script instead of here.


(2) simulation step: In this step, we first simulate the chromatin structures of low-resolution model based on stochastic folding and then convert the simulated low-resolution structures to high-resolution ones. We repeat this simulation and structural convert until 10-kb resolution structures are achieved. The script perform the calculation fo simulate the random folded polymer structure hierarchically from low resolution (5120-kb) level to high resolution level (10-kb). The outputs the script are a series of files, deposited as ./simulation/replica/1/simulation_nocompart/1/output_10.dat, ./simulation/replica/1/simulation_nocompart/2/output_10.dat, ..., ./simulation/replica/2/simulation_nocompart/1/output_10.dat, ./simulation/replica/2/simulation_nocompart/2/output_10.dat, and etc. There are in total of 500 files are generated.
You should be noted that the scripte do_simulation.sh is designed for signle cpu. If you just run it, it will take very very long time. You need to modify the script if you want it running on multiple cpus. The modification is documented in the script.
command: bash do_simulation.sh




(3) structure analysis: In this step, we first convert the files containing the coordinates of each polymer beads to the files containing the coordinates of each 10-kb segments. With the converted files, we then calculate the separation score for each 10-kb segments, based on which we can identify the genomic positions of domain boundaries.
Details of each commands are documented in the script file do_analysis.sh
command: bash do_analysis.sh



3.System Requirements
The code is tested on Linux operating systems with the system Ubuntu 16.04
It requires only a standard computer with enough RAM to support the calculation. We recommend a computer with the following specs:

RAM: 8+ GB
CPU: 4+ cores, 3.2+ GHz/core



4.Installation Guide
awk and mysql are required before running the script, which can be installed by 

sudo apt-get install mawk
sudo apt-get install mysql-server

it may take serveral minutes to finish the installation



5.Demo
Just run the three scripts one by one
bash prepareall.sh	#this step is fast and it may take several minutes to finish
bash do_simulation.sh	#this script is designed using one core and it may take several months to finish this. you can modify the script to use multiple cores or cpus to save time
bash do_analysis.sh	#it may take several one or two days to finish. you can also use modify the script to use multiple cores or cpus to save time.



6.Results
After finishing these 3 steps, the files containing coordinates of simulated chromatin regions are deposed in the directory ./analysis/conf/mergepoint/ with filename as 1.dat, 2.dat,... In these coordinate files, each row corresponds to one 10-kb segment and these 10-kb segments are listed in the sequential order. The chromosome id of each row is listed in the file ./prepare/assign/detail/chrassign_10kb.dat
The separation score values of all 10-kb segments are deposited as ./analysis/domain/boundary_delta0.20_win10_nodelta/1/score_topline.dat (chromosome 1), ./analysis/domain/boundary_delta0.20_win10_nodelta/2/score_topline.dat (chromosome 2), and etc. The genomic positions of domain boundaries are deposited as ./analysis/domain/boundary_delta0.20_win10_nodelta/1/targetpeakid.dat (chromosome 1), ./analysis/domain/boundary_delta0.20_win10_nodelta/2/targetpeakid.dat (chromosome 2), and etc. For example, the first line of ./analysis/domain/boundary_delta0.20_win10_nodelta/1/targetpeakid.dat is "77 0.653527". It means the 78-th 10-kb segment(chr 1:760kb-770kb) is identified as domain boundary with the separation score 0.653527.

