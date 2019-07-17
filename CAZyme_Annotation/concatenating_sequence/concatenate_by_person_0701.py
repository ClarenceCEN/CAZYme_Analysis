#!/usr/bin/env python
from __future__ import print_function, division
import argparse
import os
import re
import pandas as pd
import gzip
import shutil

TRUE_FALSE_DICT = {
    "True": True,
    "False": False,
    "true": True,
    "false": False,
    "T": True,
    "F": False,
    "t": True,
    "f": False,
    "1": True,
    "0": False,
    "TRUE": True,
    "FALSE": False
}

def make_arg_parser():
    parser = argparse.ArgumentParser(description='This is the commandline interface for combine_seqs',
                                     usage='combine_seqs v0.0.1 -i <input> -o <output> ...')
    parser.add_argument('-i', '--input', help='Set the directory path of the fastq directory', required=True)
    parser.add_argument('-o', '--output', help='Set the directory path of the output to place the combined_seqs.fna (default: cwd)', default=os.getcwd())
    parser.add_argument('-s', '--strip_underscore', help='Prune sample names after the first underscore (default: %(default)s)',default="False")
    parser.add_argument('-t', '--type', help='FASTQ or FASTA', choices=['FASTA', 'FASTQ','FASTA.GZ','FASTQ.GZ'], default='FASTA')
    parser.add_argument('--debug', help='Retain all intermediate files (default: Disabled)', dest='debug', action='store_true')
    parser.add_argument('--mapping_file', '-m', help='mapping file, required', required=True, type=str)
    parser.add_argument('--username', '-u', help='column name for usernames', required=True, type=str)
    parser.add_argument('--sampleid', '-id', help='column name for samples', required=True, type=str)
    return parser

def read_fasta(fh):
    """
    :return: tuples of (title, seq)
    """
    title = None
    data = None
    for line in fh:
        if line[0] == ">":
            if title:
                yield (title.strip(), data)
            title = line[1:]
            data = ''
        else:
            data += line.strip()
    if not title:
        yield None
    yield title.strip(), data


def read_fastq(fh):
    # Assume linear FASTQS
    while True:
        title = next(fh)
        #print(title)
        while title[0] != '@':
            title = next(fh)
        # Record begins
        if title[0] != '@':
            raise IOError('Malformed FASTQ files, verify they are linear and contain complete records.')
        title = title[1:].strip()
        sequence = next(fh).strip()
        garbage = next(fh).strip()
        if garbage[0] != '+':
            raise IOError('Malformed FASTQ files, verify they are linear and contain complete records.')
        qualities = next(fh).strip()
        if len(qualities) != len(sequence):
            raise IOError('Malformed FASTQ files, verify they are linear and contain complete records.')
        yield title, sequence


def format_basename(filename,run_type): # filename="/project/cs/MCT.001.S001.R1.fastq.gz"
    if run_type == 'FASTQ.GZ' or run_type == 'FASTA.GZ':
        i = -2
    else:
        i = -1

    if TRUE_FALSE_DICT[STRIP]:
        parts = os.path.basename(filename).split('_') #parts = ['MCT.001.S001.R1.fastq.gz']

        if len(parts) == 1:
            return re.sub('[^0-9a-zA-Z]+', '.', '.'.join(parts[0].split('.')[:i])) #substitude other char with '.'
        else:
            appendage = ''
            for section in parts[1:]:
                if section.find("R1") != -1:
                    appendage = 'R1'
                elif section.find("R2") != -1:
                    appendage = 'R2'
            return re.sub('[^0-9a-zA-Z]+', '.', parts[0])+appendage
    else:
        return re.sub('[^0-9a-zA-Z]+', '.', '.'.join(os.path.basename(filename).split('.')[:i]))


def convert_combine(inputs, run_type, output_path,username):
    output_filename = os.path.join(output_path, username+'.fasta')
    with open(output_filename, 'a') as outf_fasta:
        for path in inputs: #e.g. for path: /project/cs/MCT.001.S001.R1.fastq.gz in ...
            basename = format_basename(path, run_type)#delete the extension, basename = 'MCT.001.S001.R1'
            try:
                seqs = 0
                if run_type=='FASTQ.GZ' or run_type=='FASTA.GZ':
                    gz_path = path
                    file_path = os.path.splitext(path)[0]
                    print(gz_path)
                    g_file = gzip.GzipFile(gz_path,'rb')
                    with open(file_path, "wb") as f:
                        f.write(g_file.read())


                    with open(file_path) as inf:
                        if run_type == 'FASTA.GZ':
                            gen = read_fasta(inf)
                        else:
                            gen = read_fastq(inf)
                        for i, (title, seq) in enumerate(gen):
                            print('%s_%i' % (basename, i))
                            outf_fasta.write('>%s_%i %s\n%s\n' % (basename, i, title, seq))
                            seqs += 1
                else:
                    with open(path) as inf:
                        if run_type == "FASTA" or run_type == 'FASTA.GZ':
                            gen = read_fasta(inf)
                        else:
                            gen = read_fastq(inf)
                        for i, (title, seq) in enumerate(gen):
                            print('%s_%i' % (basename,i))
                            outf_fasta.write('>%s_%i %s\n%s\n' % (basename, i, title, seq))
                            seqs += 1
            except RuntimeError:
                print(basename+' from '+username+' concatenated!')
    return [output_filename]

def read_mapping_file(mapping_file_path,username_col,sampleid_col):

    mapping_file = pd.read_csv(mapping_file_path,sep='\t')
    UserNames = set(mapping_file[username_col])

    User_Sample_dict = {Username: mapping_file[mapping_file[username_col] == Username][sampleid_col].tolist() for Username in
                        UserNames}

    return User_Sample_dict

def combine_seqs(User_Sample_dict,output_path,run_type,paths):

    for username in User_Sample_dict.keys():
        sampleids = User_Sample_dict.get(username)  # e.g. MCT.f.0001, MCT.f.0002 etc
        if run_type == 'FASTQ.GZ' or run_type == 'FASTA.GZ':
            pattern = re.compile(r"^(.*)[._]S[0-9]+.*(R1|R2)\.[0-9]+\.f.*\.gz$")   #anything before 'S' or start with sample id?
        else:
            pattern = re.compile(r"^(.*)[._]S[0-9]+.*(R1|R2)\.[0-9]+\.f.*$")
        paths_for_user = []
        for path in paths:
            basename = format_basename(path, run_type)
            sampleid = re.search(pattern, basename)

            if sampleid and sampleid.group(1) in sampleids:
                print(sampleid.group(1))
                paths_for_user.append(path)

        convert_combine(paths_for_user, run_type, output_path, username)



def main():
    parser = make_arg_parser()
    args = parser.parse_args()

    global STRIP
    STRIP = args.strip_underscore
    global DEBUG
    DEBUG = args.debug

    # FIRST CHECK IF THE INPUT AND OUTPUT PATH EXIST. IF DO NOT, RAISE EXCEPTION AND EXIT
    if not os.path.exists(args.input):
        raise ValueError('Error: Input directory %s doesn\'t exist!' % args.input)

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    mode = args.type

    if mode == 'FASTA.GZ':
        paths = [os.path.join(args.input, f) for f in os.listdir(args.input) if f.endswith('fasta.gz') or f.endswith('fna.gz') or f.endswith('fa.gz') or f.endswith('fn.gz')]
    elif mode == 'FASTQ.GZ':
        paths = [os.path.join(args.input, f) for f in os.listdir(args.input) if f.endswith('fastq.gz') or f.endswith('fq.gz')]
    elif mode == 'FASTA':
        paths = [os.path.join(args.input, f) for f in os.listdir(args.input) if
                 f.endswith('fasta') or f.endswith('fna') or f.endswith('fa') or f.endswith('fn')]
    elif mode == 'FASTQ':
        paths = [os.path.join(args.input, f) for f in os.listdir(args.input) if
                 f.endswith('fastq') or f.endswith('fq')]

    user_sample_dict = read_mapping_file(args.mapping_file,args.username,args.sampleid)

    combine_seqs(user_sample_dict,args.output,mode,paths)

if __name__ == '__main__':
    main()