import pandas as pd
import numpy as np


def fix_ctrl_sheet(df):
    """
    Process the control dataset sheet by extracting relevant gene information from 'Name' column.
    Parameters:
    df (pd.DataFrame): Input DataFrame containing the control dataset.
    
    Returns:
    pd.DataFrame: Processed DataFrame with extracted gene information.
    """
    df[['Gene Name', 'Ensembl ID', 'Description']] = df['Name'].str.split('/', expand=True)
    df = df.drop(columns=['Name'], errors='ignore')
    return df


def process_sheet(df, time_label, normalized_counts, filter_num=0.05):
    """
    Process dataset sheets by handling warnings, computing IR ratios, and filtering based on splicing quality.

    Parameters:
    df (pd.DataFrame): Input dataset for a specific time point.
    time_label (str): Timepoint label (e.g., 'CTRL', '1HR', '2HR', '4HR').
    normalized_counts (pd.DataFrame): DataFrame containing normalized gene counts.
    filter_num (float): Threshold for filtering based on IR ratio.
    
    Returns:
    pd.DataFrame: Processed dataset after filtering and normalization.
    """
        
    df = df.drop(columns=['Description'], errors='ignore')  # Drop 'Description' column
    
    warnings_array = np.zeros((len(df), 3), dtype=int)
    valid_ir_columns = []
    
    for i in range(1, 4):  # 3 replicates
        ir_col = f'Irratio_{time_label}{i}'
        warn_col = f'Warnings_{time_label}{i}'
        
        if ir_col in df.columns and warn_col in df.columns:
            warn_mask = df[warn_col].str.contains('LowSplicing|LowCover', na=False)
            warnings_array[:, i-1] = warn_mask.astype(int)
            valid_ir_columns.append(ir_col)
    
    all_warnings_mask = warnings_array.sum(axis=1) == 3
    df.loc[all_warnings_mask, valid_ir_columns] = np.nan
    df['Average IRratio'] = df[valid_ir_columns].mean(axis=1, skipna=True)
    df = df.dropna(subset=['Average IRratio'])

    # Apply filtering to not CTRL (so differences in IR are more percise)
    if time_label != 'CTRL':
        df = df[abs(df['Average IRratio']) >= filter_num].copy()

    norm_col_map = {
        "CTRL": ["CTRL 1 Count", "CTRL 2 Count", "CTRL 3 Count"],
        "1HR": ["1HR 1 Count", "1HR 2 Count", "1HR 3 Count"],
        "2HR": ["2HR 1 Count", "2HR 2 Count", "2HR 3 Count"],
        "4HR": ["4HR 1 Count", "4HR 2 Count", "4HR 3 Count"]
    }

    if time_label in norm_col_map:
        norm_cols = ["Ensembl ID"] + norm_col_map[time_label]
        norm_counts_subset = normalized_counts[norm_cols]  # Subset normalized counts file
        norm_counts_subset = norm_counts_subset.rename(columns={col: col.strip() for col in norm_counts_subset.columns})

        df = df.merge(norm_counts_subset, on='Ensembl ID', how='left')
        norm_replicate_cols = norm_col_map[time_label]
        df = df[df[norm_replicate_cols].ge(1).any(axis=1)]  # Keep if at least one replicate is ≥ 1

    return df



def create_master(df_filt_list):
    """
    Merge processed datasets into a master DataFrame containing differential IR ratios.

    Parameters:
    df_filt_list (list of pd.DataFrame): List of processed DataFrames for different time points (CTRL, 1HR, 2HR, 4HR).
    
    Returns:
    pd.DataFrame: Master DataFrame containing IR ratios and differential IR ratios.
    """
    df_master = df_filt_list[0][['Gene Name', 'Ensembl ID', 'Start', 'End', 'Average IRratio']].rename(columns={'Average IRratio': 'Average IRratio CTRL'})
    df_master = df_master.merge(df_filt_list[1][['Gene Name', 'Ensembl ID', 'Start', 'End', 'Average IRratio']].rename(columns={'Average IRratio': 'Average IRratio 1HR'}), on=['Gene Name', 'Ensembl ID', 'Start', 'End'], how='outer')
    df_master = df_master.merge(df_filt_list[2][['Gene Name', 'Ensembl ID', 'Start', 'End', 'Average IRratio']].rename(columns={'Average IRratio': 'Average IRratio 2HR'}), on=['Gene Name', 'Ensembl ID', 'Start', 'End'], how='outer')
    df_master = df_master.merge(df_filt_list[3][['Gene Name', 'Ensembl ID', 'Start', 'End', 'Average IRratio']].rename(columns={'Average IRratio': 'Average IRratio 4HR'}), on=['Gene Name', 'Ensembl ID', 'Start', 'End'], how='outer')

    
    df_master['DIFF_IRratio_TNF_1HR'] = df_master['Average IRratio 1HR'] - df_master['Average IRratio CTRL']
    df_master['DIFF_IRratio_TNF_2HR'] = df_master['Average IRratio 2HR'] - df_master['Average IRratio CTRL']
    df_master['DIFF_IRratio_TNF_4HR'] = df_master['Average IRratio 4HR'] - df_master['Average IRratio CTRL']


    return df_master