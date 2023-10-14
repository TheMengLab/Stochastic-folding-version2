#!/bin/bash
chrnum=$1

for i in `seq 1 10`
do
rm -r $i
cp -r template ${i}
#cp -r 1 ${i}
cd ${i}
bash prepare.sh $i $chrnum
cd ../
done
