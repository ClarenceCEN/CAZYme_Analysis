import os
import re
import argparse
import subprocess
import pandas as pd


def make_arg_parser():
    parser = argparse.ArgumentParser(description='This is the commandline interface for burst_alignment',
                                     usage='burst_alignment v0.0.1 -i <input> -o <output> ...')
    parser.add_argument('-i', '--input', help='Set the path of the fasta directory', required=True)
    parser.add_argument('-d', '--database', help='Set the path of the database directory', required=True)
    parser.add_argument('-t', '--tax', help='Set the directory path of the taxonomy tab-limited files', required=True)
    parser.add_argument('-o', '--output', help='Set the path of the output of b6 files', default=os.getcwd())
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

def burst_command(query,a,r,tax,output):
    command = 'burst15 -q '+query+' -a '+a+' -r '+r+' -b '+tax+' -o ' + output
    print(command)
    ret = subprocess.run(command,shell=True)
    if ret.returncode==0:
        print('Success:',ret)
    else:
        print('Error:',ret)


def align_seqs(User_Sample_dict,output_path,paths,database_path,tax_path):

    for username in User_Sample_dict.keys():
        sampleids = User_Sample_dict.get(username)  # e.g. MCT.f.0001, MCT.f.0002 etc
        pattern = re.compile(r"^(.*)[._]S[0-9]+.*(R1|R2)\.[0-9]+\.")  # string before S
        a_path = os.path.join(database_path,username+'.acx')
        r_path = os.path.join(database_path,username+'.edx')
        t_path = os.path.join(tax_path,username+'.tax')
        paths_for_user = []
        for path in paths:
            basename = os.path.splitext(os.path.split(path)[1])[0]
            find_sampleid = re.search(pattern, basename)
            if find_sampleid and find_sampleid.group(1) in sampleids:
                print(find_sampleid.group(1))
                sampleid = find_sampleid.group(1)
                paths_for_user.append(path)

                if not os.path.exists(os.path.join(output_path,username)):
                    os.makedirs(os.path.join(output_path,username))

                o_path = os.path.join(output_path,username,sampleid+'.b6')
                print('BURST on '+path+' of '+username)
                burst_command(path,a_path,r_path,t_path,o_path)
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

    align_seqs(user_sample_dict,args.output,input_paths,args.database,args.tax)

if __name__ == '__main__':
    main()