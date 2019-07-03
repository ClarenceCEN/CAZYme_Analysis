import re
import pandas as pd
import os
import argparse

def make_agr_parser():
    parser = argparse.ArgumentParser(description='This is the commandline interface for building Cazymes database')
    parser.add_argument('-i1', '--input_file1', help='The input dbCAN result table', required=True)
    parser.add_argument('-i2', '--input_file2', help='The input assembled,annotated protein sequence file', required=True)
    parser.add_argument('-o', '--output', help='The output directory of Cazyme databases', required=True,default=os.getcwd())
    return parser


def build_database(table,seq,output):
    tab_file = pd.read_csv(table, sep='\t', header=None)
    tab_file.columns = ['Family_HMM', 'HMM_length', 'Query_ID', 'Query_length', 'E_value',
                        'HMM_start', 'HMM_end', 'Query_start', 'Query_end', 'Coverage']

    with open(seq, 'r') as f:
        seq_file = f.read()

    CAZYme_family = {'GT': 'GlycosylTransferases',
                     'GH': 'GlycosideHydrolases',
                     'CE': 'CarbohydrateEsterases',
                     'PL': 'PolysaccharideLyases',
                     'AA': 'AuxiliaryActivities',
                     'CBM': 'Carbohydrate-BindingModules'}

    basename = os.path.splitext(os.path.basename(seq))[0]

    output_fasta = os.path.join(output,basename)+'.fasta'
    output_txt = os.path.join(output, basename) + '.txt'

    if os.path.exists(output_fasta):
        os.remove(output_fasta)

    if os.path.exists(output_txt):
        os.remove(output_txt)

    for i in range(0, tab_file.shape[0]):
        family = tab_file.loc[i, 'Family_HMM'].replace('.hmm', '')
        query_id = tab_file.loc[i, 'Query_ID']
        query_start = tab_file.loc[i, 'Query_start']
        query_end = tab_file.loc[i, 'Query_end']
        query_length = tab_file.loc[i, 'Query_length']
        print('%d: %s' % (i + 1, family))
        seq_pattern = re.compile(r'>' + query_id + '.*?\n(.*?)(>|\n$)', re.S)
        try:
            query_sequence = re.search(seq_pattern, seq_file).group(1).replace('\n', '')  # does |n count as char?
        except:
            print('Did not find %s in sequence!' % family)
            continue
        if (len(query_sequence) != query_length):
            print("The length of %s did not match! Please check the input seqs!" % family)
            continue
        print('Query Sequence: %s' % query_sequence)

        with open(output_fasta, 'a') as f:
            f.write('>cazy_%04d_%s\n%s\n' % (i, query_id, query_sequence))

        target_sequence = query_sequence[query_start - 1:query_end]
        print('Target Seuqunce: %s\n' % target_sequence)

        # Target_sequences.append(target_sequence)

        name_pattern = re.compile(r'(GH|GT|CBM|PL|AA|CE)(\d+)(_\d+)?')
        try:
            cazy_cat = re.search(name_pattern, family).group(1)
            cazy_num = re.search(name_pattern, family).group(2)
        except:
            print('Find %s!' % family)
            print("Interesting family! What is that?")

        level_1 = 'CAZYme'
        if cazy_cat in ['GH', 'GT', 'AA', 'CE', 'PL']:
            level_2 = 'Enzyme_Classes'
            level_3 = CAZYme_family[cazy_cat]
            level_4 = '%s%s' % (cazy_cat, cazy_num)
            if cazy_cat == 'GH' and cazy_num in [5, 13, 30, 43]:
                try:
                    cazy_sub_num = re.search(name_pattern, family).group(3)
                except:
                    print('%s%s had no sub family.' % (cazy_cat, cazy_num))
                    cazy_sub_num = ''
                level_5 = '%s%s%s' % (cazy_cat, cazy_num, cazy_sub_num)
                with open(output_txt, 'a') as f:
                    f.write('cazy_%04d_example\t%s;%s;%s;%s;%s\n' % (i, level_1, level_2, level_3, level_4, level_5))
            else:
                with open(output_txt, 'a') as f:
                    f.write('cazy_%04d_example\t%s;%s;%s;%s\n' % (i, level_1, level_2, level_3, level_4))
        elif cazy_cat == 'CBM':
            level_2 = 'Associated_Modules'
            level_3 = CAZYme_family[cazy_cat]
            level_4 = '%s%s' % (cazy_cat, cazy_num)
            with open(output_txt, 'a') as f:
                f.write('cazy_%04d_example\t%s;%s;%s;%s\n' % (i, level_1, level_2, level_3, level_4))


def main():
    parser = make_agr_parser()
    args = parser.parse_args()

    dbCan_table = args.input_file1
    pro_seq = args.input_file2
    output_path = args.output

    if not os.path.exists(output_path):
        os.makedirs(output_path)

    build_database(dbCan_table,pro_seq,output_path)

if __name__ == '__main__':
    main()