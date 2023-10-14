#!/bin/bash
chrnum=$1
datafile=../data_raw/ATAC_data.dat



for i in `seq 1 $[chrnum-1]`
do
grep -P chr${i}"\t" $datafile | awk '{print $2,$3,$4}' > chr${i}.dat
done

grep -P chrX"\t" $datafile | awk '{print $2,$3,$4}' > chr${chrnum}.dat



cd bin_10kb
bash do.sh $chrnum
cd further/
bash do.sh $chrnum
cd ../
cd ../

