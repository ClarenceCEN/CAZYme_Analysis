#!/bin/bash


for f in ./*.txt
do
name=`basename $f .txt`
biom convert -i $f -o $name.biom --to-json --table-type "OTU table" --process-obs-metadata=taxonomy
done

