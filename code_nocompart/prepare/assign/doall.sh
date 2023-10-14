#!/bin/bash
chrnum=$1

bash do.sh

cd detail
bash do.sh $chrnum
cd ..
