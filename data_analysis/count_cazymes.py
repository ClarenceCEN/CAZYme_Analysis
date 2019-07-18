import pandas as pd
import os
import argparse

def make_arg_parser():
    parser = argparse.ArgumentParser(description='This is the commandline interface for cazyme_counts')
    parser.add_argument('-i', '--input', help="Set the directory path of the burst output directory, which contains folders.", required=True)
    parser.add_argument('-o', '--output', help='Set the directory path of the output of b6 files', default=os.getcwd())

    return parser

def main():
    parser = make_arg_parser()
    args = parser.parse_args()

    if not os.path.exists(args.input):
        raise ValueError('Error: Input directory %s doesn\'t exist!' % args.input)

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    input_path = os.path.abspath(args.input)
    output_path = os.path.abspath(args.output)

    #temp_mean = pd.DataFrame()

    for dir in os.listdir(input_path):
        #print(dirs)
        temp = pd.DataFrame()
        for f in os.listdir(os.path.join(input_path,dir)):
            try:
                if f.endswith(".b6"):
                    print(dir)
                    print(f)
                    f_table = pd.read_csv(os.path.join(input_path,dir,f), sep='\t', header=None)
                    f_table.loc[:, 13] = f[:-3]
                    f_table = f_table.loc[:, [12, 13]]
                    temp = pd.concat([temp,f_table],axis=0,ignore_index=True)
            except:
                continue
        temp.columns = ['CAZyme', 'SampleID']
        temp['count'] = 1
        print(temp)

        temp_wide = temp.pivot_table(index='SampleID', columns='CAZyme', values='count', aggfunc='sum')
        temp_wide[temp_wide.isnull()] = 0.0
        temp_wide.T.to_csv(os.path.join(output_path,dir)+'.txt',sep='\t')

        #temp_wide_mean = temp_wide.mean()
        #temp_wide_mean.columns = ['CAZyme',dir]

        #temp_mean = pd.concat([temp_mean, temp_wide_mean], axis=1, ignore_index=True)

    #temp_mean.to_csv(os.path.join(output_path,'CAzyme_mean.txt'))







if __name__ == '__main__':
    main()