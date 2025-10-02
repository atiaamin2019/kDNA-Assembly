# Import necessary libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kruskal


# Load the dataset
file_path = "normalized_mean_coverage_final_mini_v2.csv"  # Update with the correct file path
#file_path = "combined_normalized_mean_coverage.csv"  # Update with the correct file path
df = pd.read_csv(file_path)

# Extract unique minicircles
minicircles = df["Minicircles"].unique()
# print(minicircles)

# Identify the four groups based on sample names
group_suffixes = ["_P0", "_1X", "_2X", "_4X"]
group_order = ["P0", "1X", "2X", "4X"]
group_colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]  # Blue, Orange, Green, Red

# Dictionary to store results
kruskal_results = []

# Perform Kruskal-Wallis test
significant_minicircles = {}

for minicircle in minicircles:
    subset = df[df["Minicircles"] == minicircle]
    groups = []
    for suffix in group_suffixes:  
        group_data = subset[subset["Sample"].str.strip().str.endswith(suffix)]["Normalized_Mean_Coverage"].dropna().values
        if len(group_data) > 0:
            groups.append(group_data)

    # Perform Kruskal-Wallis test only if all four groups have data
    if len(groups) == 4:
        stat, p_value = kruskal(*groups)
        
        kruskal_results.append([minicircle, p_value])
        
        if p_value < 0.05:  # If significant, store result
            significant_minicircles[minicircle] = p_value

print("significant_minicircles")
print(significant_minicircles)
# Convert results to DataFrame
kruskal_df = pd.DataFrame(kruskal_results, columns=["Minicircle", "P-Value"])

# Filter for significantly different minicircles
significant_df = kruskal_df[kruskal_df["P-Value"] < 0.05]

# Set seaborn style for publication-quality visuals
sns.set_style("ticks")
sns.set_context("notebook")

# Generate individual box plots with filled points for each significant minicircle
for minicircle, p_value in significant_minicircles.items():
    # Subset data for the specific minicircle
    subset = df[df["Minicircles"] == minicircle].copy()

    # Map sample names to ordered group categories
    subset["Group"] = subset["Sample"].str.extract(f"({'|'.join(group_suffixes)})")[0].str.strip("_")
    subset["Group"] = pd.Categorical(subset["Group"], categories=group_order, ordered=True)

    # Create figure
    plt.figure(figsize=(8, 6))

    # Box plot with transparency (alpha) and black edges
    ax = sns.boxplot(
        x="Group", y="Normalized_Mean_Coverage", data=subset, 
        showfliers=False, linewidth=2.5, boxprops=dict(facecolor='none', edgecolor="black", linewidth=2.5)
    )

    # Scatter plot with **filled circular dots** (marker='o') and distinct colors for each group
    sns.stripplot(
        x="Group", y="Normalized_Mean_Coverage", data=subset, size=9, alpha=0.8, 
        palette=group_colors, edgecolor="black", linewidth=0.8, marker='o'  # 'o' ensures filled circular dots
    )

    # Title and labels with improved formatting
    plt.title(f"Kruskal-Wallis Test: {minicircle}\np-value = {p_value:.3e}", fontsize=14, fontweight="bold", pad=15)
    plt.xlabel("Group", fontsize=18, fontweight="bold", color = "black")
    plt.ylabel("Normalized Mean Coverage", fontsize=18, fontweight="bold", color = "black")
    plt.yticks(fontsize=14, fontweight="bold")  # Increase font size of y-axis values
    plt.xticks(fontsize=14, fontweight="bold")  # Increase font size of y-axis values


    # Remove background grid lines
    plt.grid(False)

    # Adjusting grid and frame settings
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Show the plot
    plt.show()
    
    

