#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Load the data
file_path = "data/final_output_filtered.xlsx"  # Update path if necessary
df = pd.read_excel(file_path)

# Set style for plots
sns.set(style="whitegrid")

# Create output directory if it doesn't exist
output_dir = "qc_plots_filtered"
os.makedirs(output_dir, exist_ok=True)

### 1️⃣ Missing Values Heatmap ###
plt.figure(figsize=(12, 6))
sns.heatmap(df.isnull(), cbar=False, cmap='viridis')
plt.title("Missing Data Heatmap")
plt.savefig(os.path.join(output_dir, "missing_data_heatmap.png"), dpi=300, bbox_inches="tight")
plt.close()

### 2️⃣ Distribution of IR Ratios ###
ir_columns = ["Average IRratio CTRL", "Average IRratio 1HR", "Average IRratio 2HR", "Average IRratio 4HR"]

plt.figure(figsize=(12, 6))
for col in ir_columns:
    sns.histplot(df[col].dropna(), kde=True, label=col, bins=50, alpha=0.6)
plt.legend()
plt.xlabel("IR Ratio")
plt.ylabel("Frequency")
plt.title("Distribution of IR Ratios")
plt.savefig(os.path.join(output_dir, "ir_ratio_distribution.png"), dpi=300, bbox_inches="tight")
plt.close()

### 3️⃣ Boxplot of IR Ratios (Outlier Detection) ###
plt.figure(figsize=(12, 6))
sns.boxplot(data=df[ir_columns])
plt.xticks(rotation=45)
plt.title("Boxplot of IR Ratios (Outliers Detection)")
plt.savefig(os.path.join(output_dir, "ir_ratio_boxplot.png"), dpi=300, bbox_inches="tight")
plt.close()

### 4️⃣ In-Frame vs Stop Codon Classification ###
classification_counts = df[['Out, No Stop', 'In, No Stop', 'In, In Stop']].sum()

plt.figure(figsize=(8, 6))
classification_counts.plot(kind='bar', color=['skyblue', 'orange', 'red'])
plt.title("Classification of Intron Translations")
plt.ylabel("Count")
plt.xticks(rotation=45)
plt.savefig(os.path.join(output_dir, "intron_classification.png"), dpi=300, bbox_inches="tight")
plt.close()

print(f"Plots saved to '{output_dir}/'")