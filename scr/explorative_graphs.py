#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import os
import seaborn as sns
from collections import Counter

output_dir = "sachi_plots_nofilt"
os.makedirs(output_dir, exist_ok=True)
file_path = "data/final_output_sachi_nofilt.xlsx"
df = pd.read_excel(file_path)


# Ensure intron length is treated as an integer (not a fraction)
df['intron_length'] = (df['intron_end'] - df['intron_start']).astype(int)

time_points = ['1HR', '2HR', '4HR']
for time in time_points:
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=df[f'DIFF_IRratio_TNF_{time}'], y=df['intron_length'])
    plt.xlabel("Difference in IR Ratio")
    plt.ylabel("Intron Length (bp)")
    plt.title(f"Intron Length vs IR Ratio ({time})")
    plt.ylim(df['intron_length'].min(), df['intron_length'].max())  # Adjust y-axis range dynamically
    plt.grid(True)
    plt.savefig(os.path.join(output_dir, f"intron_length_vs_diff_ir_ratio_{time}.png"), dpi=300, bbox_inches='tight')
    plt.close()

# Violin plots for IR Ratio
plt.figure(figsize=(10, 6))
sns.violinplot(data=df[['Average IRratio CTRL', 'Average IRratio 1HR', 'Average IRratio 2HR', 'Average IRratio 4HR']])
plt.xticks(range(4), ['CTRL', '1HR', '2HR', '4HR'])
plt.xlabel("Condition")
plt.ylabel("IR Ratio")
plt.title("Distribution of IR Ratios Across Conditions")
plt.grid(True)
plt.savefig(os.path.join(output_dir, "violinplot_ir_ratio.png"), dpi=300, bbox_inches='tight')
plt.close()


# Chromosome Breakdown Normalized by Gene Count
# chromosome_gene_counts = df.groupby('ref_chr')['Gene Name'].nunique()
# plt.figure(figsize=(12, 6))
# plt.bar(chromosome_gene_counts.index, chromosome_gene_counts.values)
# plt.xlabel("Chromosome")
# plt.ylabel("Number of Unique Genes")
# plt.title("Chromosome Breakdown Normalized by Gene Count")
# plt.xticks(rotation=90)
# plt.grid(axis='y', linestyle='--', alpha=0.7)
# plt.savefig(os.path.join(output_dir, "chromosome_breakdown.png"), dpi=300, bbox_inches='tight')
# plt.close()


# Histogram of Stop Codon Positions
stop_positions = df['Translation'].str.find('*')
plt.figure(figsize=(10, 6))
plt.hist(stop_positions[stop_positions >= 0], bins=20, alpha=0.7, color='red', edgecolor='black')
plt.xlabel("Stop Codon Position")
plt.ylabel("Frequency")
plt.title("Distribution of First Stop Codon Positions")
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.savefig(os.path.join(output_dir, "stop_codon_positions.png"), dpi=300, bbox_inches='tight')
plt.close()

# Classification of In-Frame and Out-of-Frame Introns (Pie Chart)
classification_counts = df[['In, No Stop', 'In, In Stop', 'Out, No Stop']].sum()
plt.figure(figsize=(8, 6))
plt.pie(classification_counts, labels=classification_counts.index, autopct="%1.1f%%", colors=['green', 'red', 'orange'], startangle=140)
plt.title("Classification of Intron Translations")
plt.savefig(os.path.join(output_dir, "classification_intron_types_pie.png"), dpi=300, bbox_inches='tight')
plt.close()


# Compute intron position relative to gene start
df['relative_position'] = (df['intron_start']) / (df['gene_end'] - df['gene_start'])

# Plot histogram of intron retention positions across all genes
plt.figure(figsize=(10, 6))
plt.hist(df['relative_position'].dropna(), bins=20, alpha=0.7, color='blue', edgecolor='black')
plt.xlabel("Relative Intron Position (0 = Start, 1 = End)")
plt.ylabel("Count of Retained Introns")
plt.title("Distribution of Retained Intron Positions Within Genes")
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.savefig(os.path.join(output_dir, "intron_retention_positions.png"), dpi=300, bbox_inches='tight')
plt.close()

# Compute CG content for each intron
df['CG_content'] = df['Intron Sequence'].apply(lambda seq: (seq.count('C') + seq.count('G')) / len(seq) if isinstance(seq, str) and len(seq) > 0 else None)

# Scatter plots: CG content vs. differential IR ratio for each time point
time_points = ['1HR', '2HR', '4HR']
for time in time_points:
    plt.figure(figsize=(8, 6))
    plt.scatter(df['CG_content'], df[f'DIFF_IRratio_TNF_{time}'], alpha=0.5, color='blue')
    plt.axvline(df['CG_content'].mean(), color='red', linestyle='--', label='Mean CG Content')
    plt.xlabel("CG Content")
    plt.ylabel(f"Δ IR Ratio ({time})")
    plt.title(f"CG Content vs. Differential IR Ratio ({time})")
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.savefig(os.path.join(output_dir, f"CG content_IR_ratio{time}.png"), dpi=300, bbox_inches='tight')
    plt.close()


df['ref_chr'] = df['ref_chr'].astype(str)

time_points = ['CTRL', '1HR', '2HR', '4HR']

human_genes_per_chromosome = {
    "1": 2100, "2": 1400, "3": 1100, "4": 750, "5": 1000, "6": 1100, "7": 950,
    "8": 750, "9": 1000, "10": 900, "11": 1300, "12": 1300, "13": 350, "14": 700,
    "15": 600, "16": 900, "17": 1250, "18": 250, "19": 1500, "20": 600, "21": 200,
    "22": 550, "X": 900, "Y": 50
}

for time in time_points:
    chromosome_gene_counts = df.groupby('ref_chr')[f'Gene Name'].nunique()
    normalized_counts = chromosome_gene_counts / chromosome_gene_counts.index.map(human_genes_per_chromosome.get)
    plt.figure(figsize=(12, 6))
    plt.bar(normalized_counts.index, normalized_counts.values)
    plt.xlabel("Chromosome")
    plt.ylabel("Unique Genes / Total Genes in Chromosome")
    plt.title(f"Normalized Chromosome Breakdown ({time})")
    plt.xticks(rotation=90)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.savefig(os.path.join(output_dir, f"chromosome_breakdown_{time}.png"), dpi=300, bbox_inches='tight')
    plt.close()




# Define color scheme
palette_colors = {"IR < 0.05": "blue", "IR ≥ 0.05": "red"}

# Load all sheets from the Excel file
time_points = ['CTRL', '1HR', '2HR', '4HR']
df_sheets = pd.read_excel(file_path, sheet_name=None)

# Define function to compute CG content
def calculate_cg_content(sequence):
    if isinstance(sequence, str) and len(sequence) > 0:
        cg_count = sequence.count('C') + sequence.count('G')
        return (cg_count / len(sequence)) * 100  # Convert to percentage
    return None

# Create empty lists for plotting
cg_content = []
time_labels = []
ir_category = []

# Process each time point sheet
for time in time_points:

    df = df_sheets[time]

    # Compute CG content for each intron
    df['CG Content (%)'] = df['Intron Sequence'].apply(calculate_cg_content)

    # Separate into two categories: IR < 0.5 and IR ≥ 0.05
    low_ir = df[df['Average IRratio'] < 0.05]['CG Content (%)'].dropna()
    high_ir = df[df['Average IRratio'] >= 0.05]['CG Content (%)'].dropna()

    # Append data to lists
    cg_content.extend(low_ir)
    time_labels.extend([time] * len(low_ir))
    ir_category.extend(["IR < 0.05"] * len(low_ir))

    cg_content.extend(high_ir)
    time_labels.extend([time] * len(high_ir))
    ir_category.extend(["IR ≥ 0.05"] * len(high_ir))

# Convert to DataFrame for Seaborn
plot_df = pd.DataFrame({'CG Content (%)': cg_content, 'Timepoint': time_labels, 'IR Category': ir_category})

# Create violin plot with only 4 x-axis labels
plt.figure(figsize=(8, 6))
sns.violinplot(x='Timepoint', y='CG Content (%)', hue='IR Category', data=plot_df,
               palette=palette_colors, inner="quartile")

# Formatting
plt.xlabel("Timepoint")
plt.ylabel("CG Content (%)")
plt.title("CG Content Distribution Across Timepoints and IR Categories")
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.legend(title="IR Category")

# Save plot
plt.savefig(os.path.join(output_dir, "cg_content_violinplot.png"), dpi=300, bbox_inches="tight")
plt.close()




print("All plots saved in", output_dir)

