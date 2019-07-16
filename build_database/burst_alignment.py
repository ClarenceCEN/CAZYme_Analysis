import os
import re

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
            pattern = re.compile(r"^(.*)[._]S[0-9]+.*(R1|R2)\.[0-9]+\.fast[aq]\.gz$")   #anything before 'S' or start with sample id?
        else:
            pattern = re.compile(r"^(.*)[._]S[0-9]+.*(R1|R2)\.[0-9]+\.fast[aq]$")
        paths_for_user = []
        for path in paths:
            basename = format_basename(path, run_type)
            sampleid = re.search(pattern, basename)

            if sampleid and sampleid.group(1) in sampleids:
                print(sampleid.group(1))
                paths_for_user.append(path)


def main():
