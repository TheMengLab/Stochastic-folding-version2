#!/bin/bash
chrnum=$1


for i in `seq 1 600`
do
rm -r $i
cp -r bk ${i}
sed -i "s/ 1/ ${i}/g" ${i}/doall.sh
done

cd ../initconf/
bash doall.sh
cd ../replica/

