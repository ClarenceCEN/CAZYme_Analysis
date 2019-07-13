#!/bin/bash

python modify_header.py


for file in ./*_modified.faa
do
name=`basename $file`
output_name="${name%.*}.out"
if [ ! -f "$output_name" ]; then
	echo $output_name
	echo '----'${name%.*}'-----'
	/plus/work/soft/hmmer-3.1/binaries/hmmscan --cpu 1 --domtblout ${name%.*}.out /plus/work/soft/CAZY/dbCAN-fam-HMMs.txt ${name%.*}.faa >${name%.*}.out2
	rm ${name%.*}.out2
fi
bash hmmscan-parser.sh ${name%.*}.out > ${name%.*}.tab
done

python build_database.py -i ./ -o output/ -l cazy_level_tab.txt