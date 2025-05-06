import requests
import time
import pandas as pd
from tqdm import tqdm
from pybiomart import Server

class EnsemblLookup:
    """Handles fetching full gene sequences and exon annotations from the Ensembl REST API."""
    
    def __init__(self, ensembl_ids):
        self.ensembl_ids = ensembl_ids
        self.gene_sequences = {}
        self.exon_data = pd.DataFrame()  # Change from dict to DataFrame
    
    def _fetch_gene_sequences_batch(self, ensembl_ids, batch_number, total_batches):
        """Fetch full gene sequences in batch from the Ensembl REST API using Ensembl IDs."""
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
                result = {seq_data['id']: (seq_data['seq'], seq_data.get('desc', '')) for seq_data in seq_data_list}
                return result
            else:
                print(f"Error fetching data for batch {batch_number}/{total_batches}: Status Code: {response.status_code}")
                print(f"Response content: {response.content.decode('utf-8')}")
                return {}
        except requests.exceptions.RequestException as e:
            print(f"RequestException for batch {batch_number}/{total_batches}: {e}")
            return {}
    
    def build_gene_sequence_cache(self):
        """Fetch Ensembl sequences in batches using a progress bar."""
        unique_ids = list(set(self.ensembl_ids))
        batch_size = 50  # Ensembl API batch limit
        total_batches = (len(unique_ids) + batch_size - 1) // batch_size

        for batch_number, i in enumerate(tqdm(range(0, len(unique_ids), batch_size), 
                                              desc="Fetching Ensembl Sequences", unit="batch"), start=1):
            batch_ids = unique_ids[i:i + batch_size]
            result = self._fetch_gene_sequences_batch(batch_ids, batch_number, total_batches)
            self.gene_sequences.update(result)
            time.sleep(1)  # Prevents excessive requests
        
        return self.gene_sequences

    def get_canonical_exon_coordinates(self, ensembl_ids, batch_size=100):
        """Fetch canonical exon annotations from Ensembl BioMart."""
        server = Server(host='http://oct2022.archive.ensembl.org')  # Use release 108

        try:
            mart = server.marts['ENSEMBL_MART_ENSEMBL'].datasets['hsapiens_gene_ensembl']
        except Exception as e:
            print(f"Error connecting to Ensembl release 108: {e}")
            return pd.DataFrame()

        results = []
        total_batches = (len(ensembl_ids) + batch_size - 1) // batch_size

        for i in tqdm(range(0, len(ensembl_ids), batch_size), 
                      desc="Fetching canonical exons", 
                      total=total_batches,
                      unit="batch"):
            batch = ensembl_ids[i:i+batch_size]

            try:
                response = mart.query(attributes=['ensembl_gene_id',
                                                  'ensembl_transcript_id',
                                                  'ensembl_exon_id',
                                                  'exon_chrom_start',
                                                  'exon_chrom_end',
                                                  'genomic_coding_start',  
                                                  'genomic_coding_end'], 
                                      filters={'link_ensembl_gene_id': batch,
                                               'transcript_is_canonical': True})

                results.append(response)
            except Exception as e:
                print(f"Error processing batch starting with {batch[0]}: {e}")
                print(f"Exception details: {str(e)}")

        if not results:
            return pd.DataFrame(columns=['ensembl_id', 'transcript_id', 'exon_id', 
                                         'exon_start', 'exon_end', 'cds_start', 'cds_end'])

        exon_df = pd.concat(results)
        exon_df.columns = ['ensembl_id', 'transcript_id', 'exon_id', 
                           'exon_start', 'exon_end', 'coding_start', 'coding_end']

        exon_df = exon_df.sort_values(['ensembl_id', 'exon_start'])

        return exon_df
    
    def run_exon_fetching(self):
        """
        Orchestrates the exon fetching using BioMart.
        Returns a DataFrame of exon annotations.
        """
        self.exon_data = self.get_canonical_exon_coordinates(self.ensembl_ids)
        return self.exon_data  # No need to call process_exon_data()

    def fetch_and_merge_ensembl(self, df_list):
        """
        - Gathers all Ensembl IDs from the DataFrames
        - Builds the gene sequence cache in batches (faster)
        - Fetches canonical exons from BioMart
        - Returns updated DF list (with sequences) + the exon DataFrame
        """
        all_ensembl_ids = set()
        for df in df_list:
            all_ensembl_ids.update(df['Ensembl ID'].dropna().astype(str))

        self.ensembl_ids = list(all_ensembl_ids)  # update internal list

        sequence_dict = self.build_gene_sequence_cache()

        exon_df = self.run_exon_fetching()

        updated_df_list = [self.merge_seqs(df, sequence_dict) for df in df_list]

        return updated_df_list, exon_df

    @staticmethod
    def merge_seqs(df, seq_dict):
        """Merge Ensembl sequences into the provided DataFrame."""
        df_gene_sequences = pd.DataFrame.from_dict(
            seq_dict, orient='index', columns=['Gene Sequence', 'Desc']
        ).reset_index()
        df_gene_sequences.rename(columns={'index': 'Ensembl ID'}, inplace=True)
        df = df.merge(df_gene_sequences, on='Ensembl ID', how='left')
        return df
