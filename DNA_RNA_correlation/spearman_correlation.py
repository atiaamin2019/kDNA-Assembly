import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr

# Load the datasets
dna_seq_df = pd.read_csv("Formatted_Dataset_DNAseq.csv")
rna_seq_df = pd.read_csv("Formatted_Dataset_RNAseq.csv")

# Extract the column names for the two groups
s_cols_dna = [col for col in dna_seq_df.columns if col.endswith('_S')]
r_cols_dna = [col for col in dna_seq_df.columns if col.endswith('_R')]

s_cols_rna = [col for col in rna_seq_df.columns if col.endswith('_S')]
r_cols_rna = [col for col in rna_seq_df.columns if col.endswith('_R')]

# Compute mean coverage for each minicircle for both groups
dna_seq_df["Mean_S_DNA"] = dna_seq_df[s_cols_dna].mean(axis=1)
dna_seq_df["Mean_R_DNA"] = dna_seq_df[r_cols_dna].mean(axis=1)

rna_seq_df["Mean_S_RNA"] = rna_seq_df[s_cols_rna].mean(axis=1)
rna_seq_df["Mean_R_RNA"] = rna_seq_df[r_cols_rna].mean(axis=1)

# Merge the computed mean values from both datasets based on minicircle names
merged_df = dna_seq_df[["Minicircles", "Mean_S_DNA", "Mean_R_DNA"]].merge(
    rna_seq_df[["Minicircles", "Mean_S_RNA", "Mean_R_RNA"]], on="Minicircles"
)

# Save merged dataframe for record
merged_df.to_csv("merged_mean_DNA_RNA_seq_cov.csv", index=False)

# ==== NEW PART ====
# Combine both groups into one series each (stack S and R together)
combined_dna = pd.concat([merged_df["Mean_S_DNA"], merged_df["Mean_R_DNA"]])
combined_rna = pd.concat([merged_df["Mean_S_RNA"], merged_df["Mean_R_RNA"]])

# Compute Spearman correlation across whole dataset
corr_all, p_value_all = spearmanr(combined_dna, combined_rna)

# Store results in a dataframe
correlation_results = pd.DataFrame({
    "Group": ["All"],
    "Correlation": [corr_all],
    "P-value": [p_value_all]
})

# Display correlation results
print("Correlation Results (Whole Dataset):")
print(correlation_results)

# ==== PLOTTING ====
# Visualization with unique colors for each minicircle
unique_minicircles = merged_df["Minicircles"].unique()
colors = sns.color_palette("hsv", len(unique_minicircles))
color_map = {minicircle: colors[i] for i, minicircle in enumerate(unique_minicircles)}

# Prepare combined dataframe for plotting
plot_df = pd.concat([
    merged_df[["Minicircles", "Mean_S_DNA", "Mean_S_RNA"]].rename(
        columns={"Mean_S_DNA": "Mean_DNA", "Mean_S_RNA": "Mean_RNA"}
    ),
    merged_df[["Minicircles", "Mean_R_DNA", "Mean_R_RNA"]].rename(
        columns={"Mean_R_DNA": "Mean_DNA", "Mean_R_RNA": "Mean_RNA"}
    )
])

# Create scatter plot
plt.figure(figsize=(8, 6))
scatter_handles = {}

# Scatter plot with different colors for each minicircle
for minicircle in unique_minicircles:
    subset = plot_df[plot_df["Minicircles"] == minicircle]
    scatter = plt.scatter(subset["Mean_DNA"], subset["Mean_RNA"], alpha=0.8,
                          edgecolors='black', marker='o', color=color_map[minicircle], label=minicircle)
    
    if minicircle not in scatter_handles:
        scatter_handles[minicircle] = scatter

# Fit a linear regression line (trend line)
m, b = np.polyfit(plot_df["Mean_DNA"], plot_df["Mean_RNA"], 1)
plt.plot(plot_df["Mean_DNA"], m * plot_df["Mean_DNA"] + b, color='red', linestyle='--', linewidth=2, label="Trend Line")

# Labels and title
plt.xlabel("Mean DNAseq Coverage", fontsize=12)
plt.ylabel("Mean RNAseq Coverage", fontsize=12)
plt.title(f"Spearman Correlation: r={corr_all:.2f}, p={p_value_all:.2e}", fontsize=14)
plt.grid(True, linestyle="--", alpha=0.5)

# Add legend
plt.legend(handles=scatter_handles.values(), labels=scatter_handles.keys(),
           loc="center left", bbox_to_anchor=(1, 0.5), title="Minicircles", fontsize=9)

# Save figure
plt.savefig("combined_correlation_analysis.pdf", format="pdf", bbox_inches="tight")

# Show plot
plt.tight_layout()
plt.show()