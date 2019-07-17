import os
import re
import argparse
import pandas as pd

def make_arg_parser():
    parser = argparse.ArgumentParser(description='This is the commandline interface for combine_seqs',
                                     usage='combine_seqs v0.0.1 -i <input> -o <output> ...')
    parser.add_argument('-i', '--input', help='Set the directory path of the fasta directory', required=True)
    parser.add_argument('-d', '--database', help='Set the directory path of the fasta directory', required=True)
    parser.add_argument('-o', '--output', help='Set the directory path of the output of b6 files', default=os.getcwd())
    parser.add_argument('--mapping_file', '-m', help='mapping file, required', required=True, type=str)
    parser.add_argument('--username', '-u', help='column name for usernames', required=True, type=str)
    parser.add_argument('--sampleid', '-id', help='column name for samples', required=True, type=str)
    return parser

def read_mapping_file(mapping_file_path,username_col,sampleid_col):

    mapping_file = pd.read_csv(mapping_file_path,sep='\t')
    UserNames = set(mapping_file[username_col])

    User_Sample_dict = {Username: mapping_file[mapping_file[username_col] == Username][sampleid_col].tolist() for Username in
                        UserNames}

    return User_Sample_dict


def align_seqs(User_Sample_dict,output_path,paths):

    for username in User_Sample_dict.keys():
        sampleids = User_Sample_dict.get(username)  # e.g. MCT.f.0001, MCT.f.0002 etc
        pattern = re.compile(r"^(.*)[._]S[0-9]+.*(R1|R2)\.[0-9]+\.")

        paths_for_user = []
        for path in paths:
            basename = os.path.splitext(os.path.split(path)[1])[0]
            sampleid = re.search(pattern, basename)
            if sampleid and sampleid.group(1) in sampleids:
                print(sampleid.group(1))
                paths_for_user.append(path)
        #print(paths_for_user)


def main():
    parser = make_arg_parser()
    args = parser.parse_args()

    if not os.path.exists(args.input):
        raise ValueError('Error: Input directory %s doesn\'t exist!' % args.input)

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    input_paths = [os.path.join(args.input, f) for f in os.listdir(args.input) if
             f.endswith('fasta') or f.endswith('fna') or f.endswith('fa') or f.endswith('fn')]

    user_sample_dict = read_mapping_file(args.mapping_file, args.username, args.sampleid)

    align_seqs(user_sample_dict,args.output,input_paths)

if __name__ == '__main__':
    main()