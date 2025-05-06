#!/usr/bin/env python3
import sys
import pandas as pd
import time
import pickle
from tqdm import tqdm

'''
This script provides the initial pickling for the excel files to be used in ir_pipeline.py

ir_pipeline.py WILL NOT WORK if you don't run this on your file first.

WILL override previous files of the same output.

Output: pickles/df_[TIMEPOINT].pkl for each sheet(timepoint) in input

ROUGH time estimate: 2 min assuming excel of 4 sheets and ~300k rows
'''



# ./pickler.py data/IRratio_Master.xlsx

# ./pickler.py data/All_Replicates_IRData.xlsx

def brine_mix(sys_file):
    input_file = sys_file[1]
    xls = pd.ExcelFile(input_file)
    
    # Define sheets and their corresponding output files
    sheets_to_pickle = [
        #('Master Sheet', 'df_master.pkl'),
        ('CTRL', 'pickles/test_ctrl.pkl'),
        ('1HR', 'pickles/test_1hr.pkl'),
        ('2HR', 'pickles/test_2hr.pkl'),
        ('4HR', 'pickles/test_4hr.pkl')
    ]
    
    for sheet_name, output_file in tqdm(sheets_to_pickle, desc="Pickling sheets"):
        df = pd.read_excel(xls, sheet_name=sheet_name)
        with open(output_file, 'wb') as f:
            pickle.dump(df, f)




if __name__ == '__main__':
    # This code is literally just for pickling
    start = time.time()
    brine_mix(sys.argv)
    print(f'Time to Pickle: {time.time() - start}')