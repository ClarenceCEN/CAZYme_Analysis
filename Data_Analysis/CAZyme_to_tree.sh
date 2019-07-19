#!/bin/bash


for f in ./*.txt
do
name=`basename $f .txt`
biom convert -i $f -o $name.biom --to-json --table-type "OTU table" --process-obs-metadata=taxonomy
done

merge_otu_tables.py -i MCTs01.biom,MCTs03.biom -o merged_otu_table.biom

for f in ./*.biom
do
name=`basename $f .biom`
if [$name != 'MCTs01' && $name != 'MCTs03'];then
echo ${name}
merge_otu_tables.py -i merged_otu_table.biom,$name.biom -o merged_otu_table.biom
fi
done
