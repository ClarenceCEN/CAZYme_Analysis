#!/bin/bash


for f in ./*.txt
do
name=`basename $f .txt`
biom convert -i $f -o $name.biom --to-json --table-type "OTU table" --process-obs-metadata=taxonomy
done

merge_otu_tables.py -i MCTs01.biom,MCTs03.biom -o merged_otu_table.biom
merge_otu_tables.py -i merged_otu_table.biom,MCTs04.biom -o merged_otu_table.biom
merge_otu_tables.py -i merged_otu_table.biom,MCTs05.biom -o merged_otu_table.biom
merge_otu_tables.py -i merged_otu_table.biom,MCTs06.biom -o merged_otu_table.biom
merge_otu_tables.py -i merged_otu_table.biom,MCTs07.biom -o merged_otu_table.biom
merge_otu_tables.py -i merged_otu_table.biom,MCTs08.biom -o merged_otu_table.biom
merge_otu_tables.py -i merged_otu_table.biom,MCTs09.biom -o merged_otu_table.biom
merge_otu_tables.py -i merged_otu_table.biom,MCTs10.biom -o merged_otu_table.biom
merge_otu_tables.py -i merged_otu_table.biom,MCTs11.biom -o merged_otu_table.biom
merge_otu_tables.py -i merged_otu_table.biom,MCTs12.biom -o merged_otu_table.biom
merge_otu_tables.py -i merged_otu_table.biom,MCTs13.biom -o merged_otu_table.biom
merge_otu_tables.py -i merged_otu_table.biom,MCTs14.biom -o merged_otu_table.biom
merge_otu_tables.py -i merged_otu_table.biom,MCTs15.biom -o merged_otu_table.biom
merge_otu_tables.py -i merged_otu_table.biom,MCTs16.biom -o merged_otu_table.biom
merge_otu_tables.py -i merged_otu_table.biom,MCTs18.biom -o merged_otu_table.biom
merge_otu_tables.py -i merged_otu_table.biom,MCTs19.biom -o merged_otu_table.biom
merge_otu_tables.py -i merged_otu_table.biom,MCTs20.biom -o merged_otu_table.biom
merge_otu_tables.py -i merged_otu_table.biom,MCTs21.biom -o merged_otu_table.biom
merge_otu_tables.py -i merged_otu_table.biom,MCTs22.biom -o merged_otu_table.biom
merge_otu_tables.py -i merged_otu_table.biom,MCTs23.biom -o merged_otu_table.biom
merge_otu_tables.py -i merged_otu_table.biom,MCTs24.biom -o merged_otu_table.biom
merge_otu_tables.py -i merged_otu_table.biom,MCTs25.biom -o merged_otu_table.biom
merge_otu_tables.py -i merged_otu_table.biom,MCTs26.biom -o merged_otu_table.biom
merge_otu_tables.py -i merged_otu_table.biom,MCTs27.biom -o merged_otu_table.biom
merge_otu_tables.py -i merged_otu_table.biom,MCTs28.biom -o merged_otu_table.biom
merge_otu_tables.py -i merged_otu_table.biom,MCTs29.biom -o merged_otu_table.biom
merge_otu_tables.py -i merged_otu_table.biom,MCTs32.biom -o merged_otu_table.biom
merge_otu_tables.py -i merged_otu_table.biom,MCTs33.biom -o merged_otu_table.biom
merge_otu_tables.py -i merged_otu_table.biom,MCTs34.biom -o merged_otu_table.biom
merge_otu_tables.py -i merged_otu_table.biom,MCTs35.biom -o merged_otu_table.biom
merge_otu_tables.py -i merged_otu_table.biom,MCTs36.biom -o merged_otu_table.biom
merge_otu_tables.py -i merged_otu_table.biom,MCTs37.biom -o merged_otu_table.biom


#for f in ./*.biom
#do
#name=`basename $f .biom`
#echo $name
#if [$name != 'MCTs01' && $name != 'MCTs03'];then
#echo ${name}
#merge_otu_tables.py -i merged_otu_table.biom,'$name'.biom -o merged_otu_table.biom
#fi
#done
