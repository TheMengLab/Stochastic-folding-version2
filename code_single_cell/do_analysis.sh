#!/bin/bash
chrnum=23

cd analysis/conf/
bash do_mergepoint.sh
cd ../..


cd analysis/domain/
bash do.sh $chrnum
cd boundary_delta0.20_win10_nodelta/
bash doall.sh
cd ../../..



