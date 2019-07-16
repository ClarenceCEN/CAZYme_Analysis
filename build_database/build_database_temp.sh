#!/bin/bash

for file in /project/flatiron2/cen/cazy_database/output/*.tax;
do
name=`basename $file .tax`
mkdir /project/flatiron2/cen/burst_output/$name
cd /project/flatiron2/cen/burst_output/$name
/project/flatiron2/cen/burst-v0.99.8-linux-64/burst15 -r /project/flatiron2/cen/cazy_database/output/$name.fasta -a $name.acx -o $name.edx -d DNA -s
done