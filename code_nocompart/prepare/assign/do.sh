#!/bin/bash

cp ../ATAC_data/bin_10kb/further/densityall.dat .

#echo -e density22.dat '\n' 1 '\n' 1 '\n' 1 '\n' 10 '\n' testmatrix_rep10.dat '\n' pointlen_rep10.dat | ./getrandomcontact_norm_pointlen.o
echo -e densityall.dat '\n' 1 '\n' 1 '\n' 1 '\n' 10 '\n' testmatrix_rep10.dat '\n' pointlen_rep10.dat | ./getrandomcontact_norm_pointlen.o

for i in `seq 1 10`
do
head -n ${i} pointlen_rep10.dat | tail -n 1 > len${i}.dat

echo -e len${i}.dat '\n' test.dat | ./row2col.o
mv test.dat len${i}.dat
done

