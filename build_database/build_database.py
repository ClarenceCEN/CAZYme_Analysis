import re
import pandas as pd
import os
import argparse

def make_agr_parser():
    parser = argparse.ArgumentParser(description='This is the commandline interface for building Cazymes database')
    parser.add_argument('-i', '--input', help='The input dbCAN result tables and protein sequences, they should have the same file name.', required=True)
    #parser.add_argument('-i2', '--input2', help='The input assembled,annotated protein sequence file', required=True)
    parser.add_argument('-o', '--output', help='The output directory of Cazyme databases', required=True,default=os.getcwd())
    parser.add_argument('-l','--level',help='The level tab file.')
    return parser

def find_value_recursive(cazy_dict,cazyme):
    global flag
    global levels
    flag = 0
    levels = cazyme
    for key,value in cazy_dict.items():
        if isinstance(value,dict):
            global n
            n += 1
            find_value_recursive(value, cazyme)
            if flag == 1:
                n -= 1
                levels =  'L'+str(n)+'_'+key+';'+levels
                return levels
            else:
                n -= 1
        elif cazyme == value:
            n += 1
            flag = 1
            levels =  'L'+str(n)+'_'+key+';'+'L'+str(n+1)+'_'+levels
            global n_level
            n_level = n+1
            return

def build_database(table,seq,output,level):
    tab_file = pd.read_csv(table, sep='\t', header=None)
    tab_file.columns = ['Family_HMM', 'HMM_length', 'Query_ID', 'Query_length', 'E_value',
                        'HMM_start', 'HMM_end', 'Query_start', 'Query_end', 'Coverage']

    with open(seq, 'r') as f:
        seq_file = f.read()

    CAZYme_family = {'CAZyme':{'Enzyme_Classes':{'GlycosylTransferases':'GT',
              'GlycosideHydrolases':'GH',
              'CarbohydrateEsterases':'CE',
              'PolysaccharideLyases':'PL',
              'AuxiliaryActivities':'AA'},
              'Associated_Modules':{'Carbohydrate-BindingModules':'CBM'}}}

    basename = os.path.splitext(os.path.basename(seq))[0]

    output_fasta = os.path.join(output,basename)+'.fasta'
    output_txt = os.path.join(output, basename) + '.txt'

    if os.path.exists(output_fasta):
        os.remove(output_fasta)

    if os.path.exists(output_txt):
        os.remove(output_txt)

    for i in range(0, tab_file.shape[0]):
        global n, n_level
        n = 0
        n_level = 0

        family = tab_file.loc[i, 'Family_HMM'].replace('.hmm', '')
        query_id = tab_file.loc[i, 'Query_ID']
        query_start = tab_file.loc[i, 'Query_start']
        query_end = tab_file.loc[i, 'Query_end']
        query_length = tab_file.loc[i, 'Query_length']
        print('%d: %s' % (i + 1, family))
        seq_pattern = re.compile(r'>' + query_id + '.*?\n(.*?)(>|\n$)', re.S)
        try:
            query_sequence = re.search(seq_pattern, seq_file).group(1).replace('\n', '')  # does \n count as char?
        except:
            print('Did not find %s in sequence!' % family)
            continue

        if (len(query_sequence) != query_length):
            print("The length of %s did not match! Please check the input seqs!" % family)
            continue
        print('Query Sequence: %s' % query_sequence)
        target_sequence = query_sequence[query_start - 1:query_end]
        print('Target Seuqunce: %s\n' % target_sequence)

        with open(output_fasta, 'a') as f:
            f.write('>cazy_%04d_%s\n%s\n' % (i, query_id, target_sequence))



        name_pattern = re.compile(r'(GH|GT|CBM|PL|AA|CE)(\d+)(_\d+)?')
        try:
            cazy_cat = re.search(name_pattern, family).group(1)
            cazy_num = re.search(name_pattern, family).group(2)
        except:
            print('Find %s!' % family)
            print("Interesting family! What is that?")
            continue

        '''
        # method 1    
        cazy_prefix = find_value_recursive(CAZYme_family,cazy_cat)
        cazy_tax = cazy_prefix+cazy_num
        if cazy_cat == 'GH' and int(cazy_num) in [5,13,30,43]:
            try:
                cazy_sub_num = re.search(name_pattern,family).group(3)
            except:
                print('%s%s had no sub family.'%(cazy_cat,cazy_num))
            cazy_tax = cazy_tax+';L'+str(n_level+1)+'_'+family
        '''

        # method 2
        cazy_tab = pd.read_csv(level,seq='\t')
        temp = cazy_tab.loc[cazy_tab.loc[:, 'L4'] == cazy_cat, :]
        cazy_tax = ''
        for i in range(0, cazy_tab.shape[1]):
            cazy_tax = cazy_tax + temp.columns[i] + '_' + temp.iloc[0, i] + ';'
        cazy_tax = cazy_tax.strip(';') + cazy_num
        if cazy_cat == 'GH' and int(cazy_num) in [5, 13, 30, 43]:
            try:
                cazy_sub_num = re.search(name_pattern, family).group(3)
            except:
                print('%s%s had no sub family.' % (cazy_cat, cazy_num))
            cazy_tax = cazy_tax + ';L5_' + family

        print(cazy_tax)
        print('\n')

        with open(output_txt, 'a') as f:
            f.write('cazy_%04d_%s\t%s\n' % (i, query_id, cazy_tax))

def main():
    parser = make_agr_parser()
    args = parser.parse_args()

    #dbCan_table = args.input_file1
    #pro_seq = args.input_file2
    input_path = args.input
    output_path = args.output
    level_path = args.level

    if not os.path.exists(output_path):
        os.makedirs(output_path)

    for f in os.listdir(input_path):
        if f.endswith("faa"):
            print(f)
            pro_seq = f
            dbCan_table = os.path.splitext(f)[0]+'.tab'
            #try:
            build_database(dbCan_table, pro_seq, output_path, level_path)
            #except:
            #    print(e)





if __name__ == '__main__':
    main()