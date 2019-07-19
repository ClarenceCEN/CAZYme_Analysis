#!/bin/bash

#go to the assembled fasta directory

cd /project/flatiron2/cen/cazy_database

for f in *.fasta;
do
mv $f ${f/_scaffolds.fasta/.fasta}
done

#The protein prediction part.
#I used GeneMark.
for file in ./.fasta
do
name=`basename $file .fasta`
perl /plus/work/soft/genemark_suite_linux_64/gmsuite/gmsn.pl --prok --format GFF --fnn --faa $name.fasta
done


python3 modify_header.py   #change the header of faa and fnn files.


#download the hmm-parser.sh

for file in ./*_modified.faa
do
name=`basename $file .faa`
output_name="$name.out"
if [ ! -f "$output_name" ]; then
	echo $output_name
	echo '----'$name'-----'
	/plus/work/soft/hmmer-3.1/binaries/hmmscan --cpu 1 --domtblout $name.out /plus/work/soft/CAZY/dbCAN-fam-HMMs.txt $name.faa >$name.out2
	rm $name.out2
fi
bash hmmscan-parser.sh $name.out > $name.tab
done

python3 build_database.py -i ./ -o output/ -l cazy_level_tab.txt

cd output/

for f in *.*;
do
mv $f ${f/_modified/}
done

mkdir /project/flatiron2/cen/burst_database
cd /project/flatiron2/cen/burst_database

for file in /project/flatiron2/cen/cazy_database/output/*.tax;
do
name=`basename $file .tax`
#mkdir /project/flatiron2/cen/burst_database/$name
#cd /project/flatiron2/cen/burst_database/$name
burst15 -r /project/flatiron2/cen/cazy_database/output/$name.fasta -a $name.acx -o $name.edx -d DNA -s
done

mkdir /project/flatiron2/cen/burst_output
cd /project/flatiron2/cen/

python3 burst_alignment.py -i /project/flatiron2/cen/dietstudy -o /project/flatiron2/cen/burst_output -d /project/flatiron2/cen/burst_database -t /project/flatiron2/cen/cazy_database/output/ -m /project/flatiron2/cen/dietstudy/food_map.txt -u UserName -id SampleID

mkdir /project/flatiron2/cen/diestudy_output
cd /project/flatiron2/cen

python3 count_cazymes.py -i burst_output/ -o Cazyme_output/