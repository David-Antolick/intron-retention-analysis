# Intron Retention Analysis Pipeline

This repository contains a modular and extensible pipeline for analyzing **intron retention (IR)** events using RNA-seq data pre-processed by [IRFinder](https://github.com/williamritchie/IRFinder). The analysis was developed for and used in the **Lee Lab** at the **University of Pittsburgh** as part of intron-centric transcriptomic research.

The pipeline integrates raw IR ratio data with Ensembl gene annotations and performs sequence-level inference to characterize retained introns based on frame status, stop codon content, and positional information. Outputs include annotated datasets and comprehensive plots to guide biological interpretation and experimental prioritization.

---

## 🧪 Background

Intron retention plays a critical role in transcriptome regulation, often leading to **nonsense-mediated decay (NMD)** or altered protein products. Understanding which introns are retained, whether they disrupt coding frames, and how they vary across conditions can provide insight into stress responses, splicing regulation, and disease mechanisms.

This pipeline automates the downstream processing of IRFinder output to:
- Filter and normalize IR events by replicate quality.
- Integrate Ensembl exon and gene annotations.
- Reconstruct intron-retained sequences.
- Translate sequences to predict frame and stop codon effects.
- Visualize IR distributions, chromosomal breakdowns, and NMD classifications.

---

## 📁 Project Structure

```plaintext
ir_analysis/
├── src/                         # All modular Python scripts
│   ├── process_ir_data.py       # IR ratio filtering and master merging
│   ├── ensembl_lookup.py        # Fetches gene sequences and exon data via Ensembl REST/BioMart
│   ├── process_sequences.py     # Extracts, translates, and annotates retained introns
│   ├── run_ir_pipeline.py       # Main pipeline script (entry point)
│   ├── analyze_coding_impact.py # Computes distances to exon junctions, stop codon logic
│   ├── prepare_pickles.py       # Converts Excel sheets into pickled files for fast loading
│   ├── plot_exploratory.py      # General CG content, stop codon, intron length plots
│   ├── plot_intron_types.py     # NMD-based classification plots
│   ├── plot_qc_metrics.py       # IR ratio and missing data visualizations
├── data/                        # Input Excel files (excluded via .gitignore)
├── pickles/                     # Pickled intermediate data (excluded)
├── plots/                       # Output figures (excluded)
├── requirements.txt
├── .gitignore
└── README.md
````

---

## ⚙️ Installation

You can run this pipeline in a **Python 3.12** environment or within a **VS Code dev container**.

### Option 1: Native Environment

```bash
git clone https://github.com/David-Antolick/intron-retention-analysis.git
cd intron-retention-analysis
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### Option 2: Dev Container (VS Code)

* Open this folder in **VS Code**
* Accept the prompt: *“Reopen in Container”*
* It will auto-install all requirements

---

## ▶️ Running the Pipeline

1. **Prepare Pickled Data**

   ```bash
   python src/prepare_pickles.py data/IRratio_Master.xlsx
   ```

2. **Run the Main Pipeline**

   ```bash
   python run_ir_pipeline.py
   ```

3. Output will be saved to:

   ```
   data/final_output.xlsx
   ```

You may also run individual visualization modules to generate figures.

---

## 📊 Example Outputs

The pipeline will generate:

* **Excel files** with annotated IR events across all timepoints
* **Plots** showing:

  * IR ratio distributions
  * Intron frame/stop classifications
  * Chromosomal distribution of retained introns
  * CG content vs. IR changes
  * NMD likelihood via stop-to-exon distance

---

## 🔍 Dependencies

All Python packages are listed in `requirements.txt`, including:

* `pandas`, `numpy`
* `matplotlib`, `seaborn`
* `biopython`
* `tqdm`, `requests`
* `pybiomart`
* `openpyxl`

---

## 🧠 Author

Developed by **David Antolick** during research in the **Lee Lab** at the **University of Pittsburgh**, 2024–2025.
For questions or collaboration inquiries, feel free to contact or visit: [github.com/David-Antolick](https://github.com/David-Antolick)

---

## 📝 License

This project does not yet specify a license. Please contact the author for use in external research or applications.

Yes — it's important (and professional) to clarify that **data are excluded due to publication restrictions**. This prevents confusion for anyone trying to run the code and signals that you're respecting confidentiality and co-authorship boundaries.

Here’s a concise section you can add to the README under its own heading:

---

### 🔒 Data Availability

This repository **does not include any input data files** (e.g., `.xlsx`, `.pkl`) due to pending publication of the associated research. All scripts and pipeline components are provided for reproducibility and code review purposes only.

If you are interested in collaborating or accessing the data, please contact the author or the corresponding members of the [Lee Lab at the University of Pittsburgh](https://www.cell.com/iscience/fulltext/S2589-0042%2821%2901301-6) after publication.

