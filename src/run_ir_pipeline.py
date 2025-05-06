#!/usr/bin/env python3
import pandas as pd
from process_ir_data import fix_ctrl_sheet, process_sheet, create_master
from ensembl_lookup import EnsemblLookup
from process_sequences import SequenceProcessor

"""
This script orchestrates the IR data  pipeline by integrating modular components for data handling, 
sequence retrieval, and sequence analysis. 

### **How It Works**
1. **Data Processing** (`data_processing.py`)
   - Loads pre-pickled datasets (CTRL, 1HR, 2HR, 4HR) and normal counts (see pickler.py)
   - Cleans and filters datasets based on intron retention (IR) ratios (filter_num) and normalized counts (seperate .xlsx file).
   - Combines data into a master file with differential IR ratios.
   
2. **Ensembl Sequence Retrieval** (`ensembl_api.py`)
   - Fetches full gene sequences from the Ensembl REST API for all unique Ensembl IDs across datasets.
   - Merges retrieved sequences into all datasets before further processing.
   
3. **Sequence Processing** (`sequence_processing.py`)
   - Extracts intron sequences based on genomic positions.
   - Translates sequences into amino acids and analyzes reading frames.
   - Identifies whether introns are in-frame and whether they contain stop codons.
   
### **Important Notes**
- Ensure `pickler.py` has been used to generate the necessary pickle files before running this script.
- This script **overwrites previous output files** (specified in block below).
- Biopython warnings about partial codons are expected and can be ignored (not all introns are in frame).

### **Output**
Processed results are saved in output_path with separate sheets for:
- Master summary dataset
- Individual time-point datasets (CTRL, 1HR, 2HR, 4HR)
- Metadata
"""


if __name__ == '__main__':
    # PARAMETERS
    filter_num = 0.05
    normal_cnts = pd.read_excel('data/Normalized_Counts_All_Replicates.xlsx')
    output_path = 'data/final_output.xlsx'

    # Load pre-pickled data
    ctrl_df = pd.read_pickle('pickles/df_ctrl.pkl')
    tnf_1hr_df = pd.read_pickle('pickles/df_1hr.pkl')
    tnf_2hr_df = pd.read_pickle('pickles/df_2hr.pkl')
    tnf_4hr_df = pd.read_pickle('pickles/df_4hr.pkl')

    # Process Control Sheet
    ctrl_df = fix_ctrl_sheet(ctrl_df)

    # Process Experimental Sheets
    processed_df_list = [
        process_sheet(ctrl_df, 'CTRL', normal_cnts, filter_num),
        process_sheet(tnf_1hr_df, '1HR', normal_cnts, filter_num),
        process_sheet(tnf_2hr_df, '2HR', normal_cnts, filter_num),
        process_sheet(tnf_4hr_df, '4HR', normal_cnts, filter_num)
    ]

    df_master = create_master(processed_df_list)
    df_master = df_master[df_master[['Average IRratio 1HR', 'Average IRratio 2HR', 'Average IRratio 4HR']].notna().any(axis=1)]
    
    # Apply filtering to CTRL dataset after master dataset creation
    processed_df_list[0] = processed_df_list[0][abs(processed_df_list[0]['Average IRratio']) >= filter_num].copy()
    
    # Collect all unique Ensembl IDs from all datasets (Master + Time Points)
    all_ensembl_ids = set()
    for df in processed_df_list + [df_master]:  # Ensure we check all datasets
        all_ensembl_ids.update(df['Ensembl ID'].dropna().astype(str))

    # Fetch Ensembl sequences and exon data for ALL datasets
    ensembl_lookup = EnsemblLookup(list(all_ensembl_ids))
    sequence_dict = ensembl_lookup.build_gene_sequence_cache()
    exon_df = ensembl_lookup.run_exon_fetching()

    # Merge sequences into ALL datasets (Master + Time Points)
    df_list_with_seqs = [SequenceProcessor(df).merge_sequences(sequence_dict).sort_dataframe()
                        .translate_sequences().analyze_sequences().get_dataframe()
                        for df in [df_master] + processed_df_list]

    # Reorder columns for clarity
    column_order = {
        "CTRL": ["Gene Name", "Ensembl ID", "Irratio_CTRL1", "Irratio_CTRL2", "Irratio_CTRL3", "Average Irratio CTRL", "ref_chr", "gene_start", "gene_end", "strand", "intron_start", "intron_end", "Intron Sequence", "Translation", "In Frame", "Stop", "Out, No Stop", "In, No Stop", "In, In Stop"],
        "1HR": ["Gene Name", "Ensembl ID", "Irratio_1HR1", "Irratio_1HR2", "Irratio_1HR3", "Average Irratio 1HR", "ref_chr", "gene_start", "gene_end", "strand", "intron_start", "intron_end", "Intron Sequence", "Translation", "In Frame", "Stop", "Out, No Stop", "In, No Stop", "In, In Stop"],
        "2HR": ["Gene Name", "Ensembl ID", "Irratio_2HR1", "Irratio_2HR2", "Irratio_2HR3", "Average Irratio 2HR", "ref_chr", "gene_start", "gene_end", "strand", "intron_start", "intron_end", "Intron Sequence", "Translation", "In Frame", "Stop", "Out, No Stop", "In, No Stop", "In, In Stop"],
        "4HR": ["Gene Name", "Ensembl ID", "Irratio_4HR1", "Irratio_4HR2", "Irratio_4HR3", "Average Irratio 4HR", "ref_chr", "gene_start", "gene_end", "strand", "intron_start", "intron_end", "Intron Sequence", "Translation", "In Frame", "Stop", "Out, No Stop", "In, No Stop", "In, In Stop"]
    }

    for i, label in enumerate(["CTRL", "1HR", "2HR", "4HR"], start=1):
        df_list_with_seqs[i] = df_list_with_seqs[i][column_order[label]]
    
    # Save to Excel
    with pd.ExcelWriter(output_path) as writer:
        df_list_with_seqs[0].to_excel(writer, sheet_name='Master', index=False)
        df_list_with_seqs[1].to_excel(writer, sheet_name='CTRL', index=False)
        df_list_with_seqs[2].to_excel(writer, sheet_name='1HR', index=False)
        df_list_with_seqs[3].to_excel(writer, sheet_name='2HR', index=False)
        df_list_with_seqs[4].to_excel(writer, sheet_name='4HR', index=False)
        exon_df.to_excel(writer, sheet_name='Exon Data', index=False)