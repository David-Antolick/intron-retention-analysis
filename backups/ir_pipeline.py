#!/usr/bin/env python3
import pandas as pd
import requests
import time
import sys
import pickle
from Bio.Seq import Seq
from tqdm import tqdm
import numpy as np

'''
This file is the combined IR data sorting and processing pipeline. Before you run this, make sure you've used
pickler.py to create the pickle files this file will use (otherwise it will break very fast).

Be careful, as the output from this WILL override any previous versions.

Don't mind any warnings from biopython about partial codons, its supposed to be that way.

Output: data/final_output.xlsx
'''




#########################
   # Sorting Stage #
#########################
def fix_ctrl_sheet(df):
    """Process individual dataset sheets by filtering based on differential IR ratio."""

    df[['Gene Name', 'Ensembl ID', 'Description']] = df['Name'].str.split('/', expand=True)

    df = df.drop(columns=['Name'], errors='ignore')

    return df

def process_sheet(df, time_label, filter_num=0.01):
    """Process dataset sheets by filtering based on differential IR ratio, dynamically handling time points."""

    df = df.drop(columns=['Description'], errors='ignore')  # Drop unused 'Description' column

    # Dynamically identify relevant columns based on the time point label
    valid_ir_columns = []

    for i in range(1, 4):  # Assuming up to 3 replicates
        ir_col = f'Irratio_{time_label}{i}'
        warn_col = f'Warnings_{time_label}{i}'

        if ir_col in df.columns and warn_col in df.columns:
            # Set affected IR values to NaN instead of dropping entire rows
            df.loc[df[warn_col].str.contains('LowCover', na=False), ir_col] = np.nan
            df.loc[df[warn_col].str.contains('LowSplicing', na=False), ir_col] = np.nan

            valid_ir_columns.append(ir_col)  # Track valid IR columns

    # Calculate the average IR ratio across valid columns
    df['Average IRratio'] = df[valid_ir_columns].mean(axis=1, skipna=True)
    df = df.dropna(subset=['Average IRratio'])


    # Apply significance filtering
    if time_label != 'CTRL' :
        df = df[abs(df['Average IRratio']) >= filter_num].copy()

    return df

def create_master(df_filt_list):
    """Merge sorted data into a master dataframe with differential IR ratios."""
    df_master = df_filt_list[0][['Gene Name', 'Ensembl ID', 'Start', 'End', 'Average IRratio']].rename(columns={'Average IRratio': 'Average IRratio CTRL'})
    df_master = df_master.merge(df_filt_list[1][['Gene Name', 'Ensembl ID', 'Start', 'End', 'Average IRratio']].rename(columns={'Average IRratio': 'Average IRratio 1HR'}), on=['Gene Name', 'Ensembl ID', 'Start', 'End'], how='outer')
    df_master = df_master.merge(df_filt_list[2][['Gene Name', 'Ensembl ID', 'Start', 'End', 'Average IRratio']].rename(columns={'Average IRratio': 'Average IRratio 2HR'}), on=['Gene Name', 'Ensembl ID', 'Start', 'End'], how='outer')
    df_master = df_master.merge(df_filt_list[3][['Gene Name', 'Ensembl ID', 'Start', 'End', 'Average IRratio']].rename(columns={'Average IRratio': 'Average IRratio 4HR'}), on=['Gene Name', 'Ensembl ID', 'Start', 'End'], how='outer')


    #df_master['Average IRratio CTRL'] = df_master['Average IRratio CTRL'].fillna(0)
    
    df_master['DIFF_IRratio_TNF_1HR'] = df_master['Average IRratio 1HR'] - df_master['Average IRratio CTRL']
    df_master['DIFF_IRratio_TNF_2HR'] = df_master['Average IRratio 2HR'] - df_master['Average IRratio CTRL']
    df_master['DIFF_IRratio_TNF_4HR'] = df_master['Average IRratio 4HR'] - df_master['Average IRratio CTRL']

    return df_master





#########################
  # Processing Stage #
#########################

# NEED to add the intron fulling feature
class EnsemblLookup:
    """Handles fetching full gene sequences from the Ensembl REST API."""
    def __init__(self, ensembl_ids):
        self.ensembl_ids = ensembl_ids
        self.gene_sequences = {}
    
    def _fetch_gene_sequences_batch(self, ensembl_ids, batch_number, total_batches):
        #  Fetch full gene sequences in batch from the Ensembl REST API using Ensembl IDs
        url = 'https://rest.ensembl.org/sequence/id'
        headers = {
            "Content-Type": "application/json",
            "Accept": "application/json",
            "User-Agent": "ir_processing/1.0"
        }
        data = {'ids': ensembl_ids, 'type': 'genomic'}

        try:
            response = requests.post(url, headers=headers, json=data, timeout=60)
            if response.status_code == 200:
                seq_data_list = response.json()
                result = {}

                for seq_data in seq_data_list:
                    ensembl_id = seq_data['id']
                    seq = seq_data['seq']
                    desc = seq_data.get('desc', '')
                    result[ensembl_id] = (seq, desc)
                return result
            
            else:
                print(f"Error fetching data for batch {batch_number}/{total_batches}: Status Code: {response.status_code}")
                print(f"Response content: {response.content.decode('utf-8')}")
                return {}
            
        except requests.exceptions.RequestException as e:
            print(f"RequestException for batch {batch_number}/{total_batches}: {e}")
            return {}

    def build_gene_sequence_cache(self):
        """Fetch Ensembl sequences in batches using TQDM progress bar."""
        unique_ids = list(set(self.ensembl_ids))
        batch_size = 50  # Ensembl API batch limit
        total_batches = (len(unique_ids) + batch_size - 1) // batch_size

        for batch_number, i in enumerate(tqdm(range(0, len(unique_ids), batch_size), desc="Fetching Ensembl Sequences", unit="batch"), start=1):
            batch_ids = unique_ids[i:i + batch_size]
            result = self._fetch_gene_sequences_batch(batch_ids, batch_number, total_batches)
            self.gene_sequences.update(result)
            time.sleep(1)  # Be gentle, reccomend 1

        return self.gene_sequences
    
def fetch_and_merge_ensembl(df_list):
    """Fetch Ensembl sequences once and merge into all DataFrames in df_list."""
    
    # Step 1: Collect all unique Ensembl IDs across all sheets
    all_ensembl_ids = set()
    for df in df_list:
        all_ensembl_ids.update(df['Ensembl ID'].dropna().astype(str))

    # Step 2: Run Ensembl lookup once for all unique IDs
    ensembl_lookup = EnsemblLookup(list(all_ensembl_ids))
    sequence_dict = ensembl_lookup.build_gene_sequence_cache()

    # Step 3: Merge Ensembl sequences into all sheets
    df_list = [merge_seqs(df, sequence_dict) for df in df_list]

    return df_list


def merge_seqs(df, seq_dict):
    # Creates df from sequence dict
    df_gene_sequences = pd.DataFrame.from_dict(
        seq_dict, orient='index', columns=['Gene Sequence', 'Desc']
    ).reset_index()

    df_gene_sequences.rename(columns={'index': 'Ensembl ID'}, inplace=True)

    # Merges them correlated by Ensembl OD
    df = df.merge(df_gene_sequences, on='Ensembl ID', how='left')

    return df


def _reverse_comp(row):
    # Used with .apply() in the the merge_seqs() func
    if row['strand'] == '-1':
        return str((Seq(row['Intron Sequence']).reverse_complement()))
    else: return row['Intron Sequence']


def _get_intron_seq(row):
    # Extract intron sequences  
        seq = row['Gene Sequence']
        start = row['intron_start']
        end = row['intron_end']
        if pd.isnull(seq) or pd.isnull(start) or pd.isnull(end):
            return None
        else:
            return str(seq[start+1:end+1])


def sort_df(df):

    # Handles all the column renaming and creation

    # Extract gene_start from 'Desc'
    # Example : 'chromosome:GRCh38:19:50415799:50428409:1
    df['ref_chr'] = df['Desc'].str.split(':').str[2]
    df['gene_start'] = df['Desc'].str.split(':').str[3].astype(float)
    df['gene_end'] = df['Desc'].str.split(':').str[4].astype(float)
    df['strand'] = df['Desc'].str.split(':').str[5].astype(float)
    df = df.drop('Desc', axis=1)


    # Compute intron start and end indices
    df['intron_start'] = (df['Start'] - df['gene_start']).clip(lower=0).astype(int)
    df['intron_end'] = (df['End'] - df['gene_start']).clip(lower=0).astype(int)


    # Applies and corrects intron sequence if it is on the negative strand (takes reverse compliment)
    df['Intron Sequence'] = df.apply(_get_intron_seq, axis=1)
    df['Intron Sequence'] = df.apply(_reverse_comp, axis=1)
    

    return df




#########################
# Post-Processing Stage #
#########################
def translate(df):
    """Translate intron sequences into amino acids."""
    df['Translation'] = df['Intron Sequence'].apply(lambda seq: str(Seq(str(seq)).translate()) if isinstance(seq, str) else None)
    return df

def analyze(df):
    """Analyze frame status and stop codon presence."""
    df['In Frame'] = df['Intron Sequence'].apply(lambda seq: len(seq.strip()) % 3 == 0 if isinstance(seq, str) else False).astype(int)
    df['Stop'] = df['Translation'].apply(lambda trans: '*' in trans if isinstance(trans, str) else False).astype(int)
    df['Out, No Stop'] = ((df['In Frame'] == 0) & (df['Stop'] == 0)).astype(int)
    df['In, No Stop'] = ((df['In Frame'] == 1) & (df['Stop'] == 0)).astype(int)
    df['In, In Stop'] = ((df['In Frame'] == 1) & (df['Stop'] == 1)).astype(int)
    return df





if __name__ == '__main__':
    # PARAMS
    filter_num = 0.05   

    output_path = 'data/final_output.xlsx'


    # Sorting Phase


    # assigns pre-pickled files to variables
    ctrl_df = pd.read_pickle('df_ctrl.pkl')
    tnf_1hr_df = pd.read_pickle('df_1hr.pkl')
    tnf_2hr_df = pd.read_pickle('df_2hr.pkl')
    tnf_4hr_df = pd.read_pickle('df_4hr.pkl')


    # fixes weird control sheet inconsistency
    ctrl_df = fix_ctrl_sheet(ctrl_df)

    df_list = [
        process_sheet(ctrl_df, 'CTRL', filter_num),
        process_sheet(tnf_1hr_df, '1HR', filter_num),
        process_sheet(tnf_2hr_df, '2HR', filter_num),
        process_sheet(tnf_4hr_df, '4HR', filter_num)
    ]
    df_master = create_master(df_list)

    df_master = df_master[df_master[['Average IRratio 1HR', 'Average IRratio 2HR', 'Average IRratio 4HR']].notna().any(axis=1)]

    df_list_with_seqs = fetch_and_merge_ensembl([df_master] + df_list)

    updated_dfs = []
    for df in df_list_with_seqs:
        df = sort_df(df)
        df = translate(df)
        df = analyze(df)
        updated_dfs.append(df)
  
    
updated_dfs[1] = updated_dfs[1][abs(updated_dfs[1]['Average IRratio']) >= filter_num].copy()



with pd.ExcelWriter(output_path) as writer:
    updated_dfs[0].to_excel(writer, sheet_name='Master', index=False)
    updated_dfs[1].to_excel(writer, sheet_name='CTRL', index=False)
    updated_dfs[2].to_excel(writer, sheet_name='1HR', index=False)
    updated_dfs[3].to_excel(writer, sheet_name='2HR', index=False)
    updated_dfs[4].to_excel(writer, sheet_name='4HR', index=False)