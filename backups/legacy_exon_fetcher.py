import requests
import pandas as pd
import time
from tqdm import tqdm

class EnsemblExonFetcher:
    """Class to fetch and process exon annotations from the Ensembl API."""
    def __init__(self, ensembl_ids):
        self.ensembl_ids = ensembl_ids
        self.exon_data = {}
    
    def get_exon_annotations(self, ensembl_id):
        """Fetch exon annotations for a given Ensembl gene ID."""
        url = f'https://rest.ensembl.org/overlap/id/{ensembl_id}?feature=exon'
        headers = {"Content-Type": "application/json", "Accept": "application/json"}
        
        try:
            response = requests.get(url, headers=headers, timeout=30)
            if response.status_code == 200:
                return response.json()
            else:
                print(f"Error fetching exons for {ensembl_id}: {response.status_code}")
                return None
        except requests.exceptions.RequestException as e:
            print(f"Request failed for {ensembl_id}: {e}")
            return None
    
    def fetch_all_exon_data(self):
        """Fetch exon annotations for all Ensembl IDs."""
        for ensembl_id in tqdm(self.ensembl_ids, desc="Fetching exon data"):
            exons = self.get_exon_annotations(ensembl_id)
            if exons:
                self.exon_data[ensembl_id] = exons
            time.sleep(1)  # Avoid overwhelming the API
        
    def process_exon_data(self):
        """Convert exon data into a structured DataFrame."""
        rows = []
        for ensembl_id, exons in self.exon_data.items():
            for exon in exons:
                rows.append({
                    "Ensembl ID": ensembl_id,
                    "Chromosome": exon.get("seq_region_name"),
                    "Exon Start": exon.get("start"),
                    "Exon End": exon.get("end"),
                    "Strand": exon.get("strand"),
                    "Exon ID": exon.get("id")
                })
        
        return pd.DataFrame(rows)
    
    def run(self):
        """Run the entire pipeline: fetch exon data and process it into a DataFrame."""
        self.fetch_all_exon_data()
        return self.process_exon_data()

# Example integration:
# unique_ensembl_ids = df_list[0]["Ensembl ID"].dropna().unique().tolist()
# exon_fetcher = EnsemblExonFetcher(unique_ensembl_ids)
# exon_df = exon_fetcher.run()
# exon_df.to_csv("data/exon_annotations.csv", index=False)
