#!/bin/bash

repid=$1

#sed -i "s/sh 1 1/sh ${repid} 1/g" -- 1/doall.sh


for i in `seq 1 50`
do
rm -r $i
#cp 1/do*.sh ${i}/
cp -r bk ${i}
done



for i in `seq 1 50`
do
sed -i "s/sh 1 1/sh ${repid} ${i}/g" -- ${i}/doall.sh
sed -i "s/sh 1$/sh ${repid}/g" -- ${i}/doall.sh
done
