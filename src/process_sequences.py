import pandas as pd
from Bio.Seq import Seq

class SequenceProcessor:
    """
    Handles sequence-related transformations, including merging Ensembl sequences,
    sorting introns, translating sequences, and analyzing frame status.
    
    Attributes:
        df (pd.DataFrame): DataFrame containing sequence data for processing.
    """
    
    def __init__(self, df):
        self.df = df
    

    def merge_sequences(self, seq_dict):
        """
        Merges Ensembl sequences into the DataFrame based on Ensembl ID.
        
        Parameters:
            seq_dict (dict): Mapping of Ensembl ID to (sequence, description) pairs.
        
        Returns:
            SequenceProcessor: The updated SequenceProcessor instance.
        """

        df_gene_sequences = pd.DataFrame.from_dict(
            seq_dict, orient='index', columns=['Gene Sequence', 'Desc']
        ).reset_index()
        df_gene_sequences.rename(columns={'index': 'Ensembl ID'}, inplace=True)
        self.df = self.df.merge(df_gene_sequences, on='Ensembl ID', how='left')
        return self
    

    def sort_dataframe(self):
        """
        Extracts and processes intron sequences based on genomic locations.
        
        Returns:
            SequenceProcessor: The updated SequenceProcessor instance.
        """
        self.df['ref_chr'] = self.df['Desc'].str.split(':').str[2]
        self.df['gene_start'] = self.df['Desc'].str.split(':').str[3].astype(float)
        self.df['gene_end'] = self.df['Desc'].str.split(':').str[4].astype(float)
        self.df['strand'] = self.df['Desc'].str.split(':').str[5].astype(float)
        self.df = self.df.drop('Desc', axis=1)

        self.df['intron_start'] = (self.df['Start'] - self.df['gene_start']).clip(lower=0).astype(int)
        self.df['intron_end'] = (self.df['End'] - self.df['gene_start']).clip(lower=0).astype(int)
        
        self.df['Intron Sequence'] = self.df.apply(self._get_intron_seq, axis=1)
        self.df['Intron Sequence'] = self.df.apply(self._reverse_comp, axis=1)
    
        return self

    def translate_sequences(self):
        """
        Translates intron sequences into amino acids.
        
        Returns:
            SequenceProcessor: The updated SequenceProcessor instance.
        """
        self.df['Translation'] = self.df['Intron Sequence'].apply(
            lambda seq: str(Seq(str(seq)).translate()) if isinstance(seq, str) else None
        )
        return self
    
    
    def analyze_sequences(self):
        """
        Determines if introns are in-frame and whether they contain stop codons.
        
        Returns:
            SequenceProcessor: The updated SequenceProcessor instance.
        """
        self.df['In Frame'] = self.df['Intron Sequence'].apply(
            lambda seq: len(seq.strip()) % 3 == 0 if isinstance(seq, str) else False
        ).astype(int)
        self.df['Stop'] = self.df['Translation'].apply(
            lambda trans: '*' in trans if isinstance(trans, str) else False
        ).astype(int)
        self.df['Out, No Stop'] = ((self.df['In Frame'] == 0) & (self.df['Stop'] == 0)).astype(int)
        self.df['In, No Stop'] = ((self.df['In Frame'] == 1) & (self.df['Stop'] == 0)).astype(int)
        self.df['In, In Stop'] = ((self.df['In Frame'] == 1) & (self.df['Stop'] == 1)).astype(int)
        return self
    

    def get_dataframe(self):
        """Return the processed DataFrame."""
        return self.df
    
    @staticmethod
    def _reverse_comp(row):
        """
        Reverse complements sequence if on the negative strand.
        
        Parameters:
            row (pd.Series): Row containing sequence and strand information.
        
        Returns:
            str: Reverse complemented sequence if on the negative strand, otherwise original sequence.
        """
        if row['strand'] == -1:
            return str(Seq(row['Intron Sequence']).reverse_complement())
        return row['Intron Sequence']
    
    
    @staticmethod
    def _get_intron_seq(row):
        """
        Extracts intron sequences from the full gene sequence.
        
        Parameters:
            row (pd.Series): Row containing gene sequence and intron start/end positions.
        
        Returns:
            str or None: Extracted intron sequence or None if missing data.
        """
        seq = row['Gene Sequence']
        start = row['intron_start']
        end = row['intron_end']
        if pd.isnull(seq) or pd.isnull(start) or pd.isnull(end):
            return None
        return str(seq[int(start)+1:int(end)+1])
