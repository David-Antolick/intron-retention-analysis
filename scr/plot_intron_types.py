#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt

# Thresholds (in amino acids) for deciding if a stop codon is likely to trigger NMD.
THRESHOLD_STOP_INTRON = 30
THRESHOLD_STOP_EXON   = 55

def classify_intron_new(row, thresh_stop_intron=THRESHOLD_STOP_INTRON, thresh_stop_exon=THRESHOLD_STOP_EXON):
    """
    Classify each intron (with the intron retained in the transcript) into one of three types:
      - "Outside coding sequence": No in-frame stop codon is found (i.e. the retained intron does not cause a stop).
      - "Non-NMD stop": An in-frame stop codon is present but the distances from the stop to the end of the intron 
                        (or to the following exon junction) are below the thresholds.
      - "Likely NMD": An in-frame stop codon is present and at least one of the distances is equal to or exceeds its threshold.
    
    This function expects the row to include:
      - "Stop": 1 if an in-frame stop codon is present, 0 otherwise.
      - "stop_to_intron_end": distance (in amino acids) from the first stop codon to the intron end.
      - "stop_to_exon_junction": distance (in amino acids) from the first stop codon to the next exon junction.
    
    If any values are missing or non-numeric, the row is classified as "Unclassified".
    """
    try:
        stop_present = int(row['Stop'])
    except (ValueError, TypeError, KeyError):
        return "Unclassified"
    
    if stop_present == 0:
        return "Outside coding sequence"
    elif stop_present == 1:
        # Retrieve distances; if missing, treat as 0 for safety.
        try:
            dist_intron = float(row['stop_to_intron_end']) if pd.notna(row['stop_to_intron_end']) else 0
        except (ValueError, TypeError):
            dist_intron = 0
        try:
            dist_exon = float(row['stop_to_exon_junction']) if pd.notna(row['stop_to_exon_junction']) else 0
        except (ValueError, TypeError):
            dist_exon = 0
        # If either distance meets or exceeds its threshold, flag as likely NMD.
        if dist_intron >= thresh_stop_intron or dist_exon >= thresh_stop_exon:
            return "Likely NMD"
        else:
            return "Non-NMD stop"
    else:
        return "Unclassified"

def main():
    # Read in the final Excel file (adjust filename if necessary)
    excel_file = "data/final_output_updated.xlsx"
    sheets = pd.read_excel(excel_file, sheet_name=None)
    
    # Specify the timepoint sheet names.
    timepoints = ["CTRL", "1HR", "2HR", "4HR"]
    type_counts = {}  # To store counts for each type per timepoint
    
    for tp in timepoints:
        if tp in sheets:
            df = sheets[tp].copy()
            # Classify each row using the new scheme.
            df['Type'] = df.apply(classify_intron_new, axis=1)
            counts = df['Type'].value_counts()
            type_counts[tp] = counts
            sheets[tp] = df  # (Optional: save the classified DataFrame back)
    
    # Build a DataFrame for the grouped bar graph.
    types = ["Outside coding sequence", "Non-NMD stop", "Likely NMD"]
    plot_df = pd.DataFrame(index=timepoints, columns=types).fillna(0)
    for tp in timepoints:
        for t in types:
            if tp in type_counts and t in type_counts[tp]:
                plot_df.loc[tp, t] = type_counts[tp][t]
    plot_df = plot_df.astype(int)
    
    # ---------------------------
    # Grouped Bar Graph (Counts)
    # ---------------------------
    ax = plot_df.plot(kind="bar", figsize=(10,6))
    ax.set_ylabel("Count")
    ax.set_title("Intron Types by Timepoint (Counts)")
    plt.xticks(rotation=0)
    plt.legend(title="Intron Type")
    plt.tight_layout()
    # Save the bar graph in the "plots" directory.
    plt.savefig("plots/intron_type_bar.png")
    plt.close()
    
    # ---------------------------
    # Pie Charts (one per timepoint)
    # ---------------------------
    fig, axes = plt.subplots(2, 2, figsize=(12,10))
    axes = axes.flatten()
    for i, tp in enumerate(timepoints):
        if tp in type_counts:
            counts = type_counts[tp].reindex(types, fill_value=0)
            axes[i].pie(counts, labels=types, autopct="%1.1f%%", startangle=90)
            axes[i].set_title(f"{tp} Intron Types")
    plt.tight_layout()
    plt.savefig("plots/intron_type_pie.png")
    plt.close()

if __name__ == "__main__":
    main()
