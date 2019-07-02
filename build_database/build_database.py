import re
import pandas as pd


tab_file = pd.read_csv('./example.tab',sep='\t',header=None)


tab_file.columns = ['Family_HMM','HMM_length','Query_ID','Query_length','E_value',
                    'HMM_start','HMM_end','Query_start','Query_end','Coverage']



Query_ID = 3 - 1 # Query ID column
Query_start = 8 - 1 # Query start column
Query_end = 9 - 1 # Query start column