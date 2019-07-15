#!/bin/bash


for f in *.fasta;
do
mv $f ${f/_scaffolds.fasta/.fasta}
done

#The protein prediction part.
#I used GeneMark.
for file in ./fasta
do
name=`basename $file`
perl /plus/work/soft/genemark_suite_linux_64/gmsuite/gmsn.pl --prok --format GFF --fnn --faa ${name%.*}.fasta
done


python modify_header.py


#download the hmm-parser.sh

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