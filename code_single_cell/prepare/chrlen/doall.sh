#!/bin/bash
#wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
chrnum=$1
datafile=../data_raw/hg19.chrom.sizes

for i in `seq 1 $[chrnum-1]`
do
grep -P chr${i}"\t" $datafile | awk '{print $2}' > chrlen${i}.dat
done 


grep -P chrX"\t" $datafile | awk '{print $2}' > chrlen${chrnum}.dat



cd tempassign_10kb
bash do.sh $chrnum
cd ..

