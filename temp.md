
###Step 1 Shi7

Sorry I lost my Shi7 code.

###Step 2 Concatenate sequence by subject.

See `concatenate_by_person_0701.py`. It takes 6 arguments: `-i`, the path of the fastq directory, `-o`: the output path to place the combined_seqs.fna (default: cwd), `-t`: 'FASTQ or FASTA', `-m`: the path of mapping file, `-u`: column name for usernames, `-id`: column name for samples.

###Step 3 Run Spades

To assemble the shi7'd and combined fasta files. See `run_spades.pbs`.

###Step 4 Protein prediction.

We were supposed to use Prokka but it did not work. So I used GeneMark.
The header of the result fasta file from GeneMark contains "|", which did not work well with regex expression so I used `modify_header.py` to change the headers. 


###Step 5 Annotate the sequences with Cazyme database
Again, I did this step on my university's server. The output file is similar with blast output file (.tab). See `CAZyme_Annotation/build_database/example/modified_example.tab`


###Step 6 Prepare files for building reference database for BURST

See `build_database.py`. It takes 4 arguments, `-i`: The input Cazyme annotation result tables and DNA sequences, they should have the same file name except for their extension, `-o`: The output directory of Cazyme databases, `-l`: The level tab file. See `CAZyme_Annotation/build_database/cazy_level.tab`, `--multi`: Whether keep multiple blast results? The default is false, which mean we only include the result with largest blast Coverage value.


###Step 7 Build BURST database file
I wrote a for loop for each subject's files. See `CAZyme_Annotation/build_database/cazy.sh` Step 7 section.


###Step 8 BURST
Use `burst_alignment.py`. It takes 7 arguments. `-i`: The path for each subject's fasta files, `-d`: The path for burst database files, `-t`: The path of tax files for burst. (Generated in step 6), `-o`:The output directory, `-m`: the path of mapping file, `-u`: the 
column name for usernames, `-id`: the column name for samples.


###Step 9 Transformation
In this step we count the Cazyme hits. Use `count_cazymes.py`. It takes 2 arguments. `-i`: The path for burst output file, which contains folders, `-o`: The path for output files. 

###Step 10 Merge the output files
Use QIIME/Python/R to merge all subjects' output file.