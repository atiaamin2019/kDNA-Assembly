import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

# Load the data
file_path = 'normalized_coverage_with_new_minicircle_ids.csv'
df = pd.read_csv(file_path)

# Group labeling
df['Group'] = df['Sample'].apply(lambda x: 'Sensitive' if x.endswith('_S') else ('Resistant' if x.endswith('_R') else 'Unknown'))

# Run Mann–Whitney U test per minicircle
results = []
for mc in df['Minicircle'].unique():
    group_data = df[df['Minicircle'] == mc]
    sens = group_data[group_data['Group'] == 'Sensitive']['Normalized_Mean_Coverage']
    res = group_data[group_data['Group'] == 'Resistant']['Normalized_Mean_Coverage']
    if len(sens) > 0 and len(res) > 0:
        stat, p = mannwhitneyu(sens, res, alternative='two-sided')
        results.append({'Minicircle': mc, 'p_value': p})

results_df = pd.DataFrame(results)
significant_minicircles = results_df[results_df['p_value'] < 0.05]['Minicircle'].tolist()
filtered_df = df[df['Minicircle'].isin(significant_minicircles)]

# Manual custom order
special = ['LINF min 6', 'LINF min 17', 'LINF min 41']
others = sorted([m for m in significant_minicircles if m not in special])
ordered_minicircles = others + special

# Plot manually
plt.figure(figsize=(max(18, len(ordered_minicircles) * 0.8), 10))

# Loop through and plot each box manually to control position
for i, mc in enumerate(ordered_minicircles):
    data = filtered_df[filtered_df['Minicircle'] == mc]
    
    sens = data[data['Group'] == 'Sensitive']['Normalized_Mean_Coverage']
    res = data[data['Group'] == 'Resistant']['Normalized_Mean_Coverage']
    
    # Boxplot for Sensitive
    b1 = plt.boxplot(sens, positions=[i - 0.2], widths=0.35,
                     patch_artist=True, boxprops=dict(facecolor='none', color='black'),
                     medianprops=dict(color='black'), capprops=dict(color='black'),
                     whiskerprops=dict(color='black'), flierprops=dict(marker=''))

    # Boxplot for Resistant
    b2 = plt.boxplot(res, positions=[i + 0.2], widths=0.35,
                     patch_artist=True, boxprops=dict(facecolor='none', color='black'),
                     medianprops=dict(color='black'), capprops=dict(color='black'),
                     whiskerprops=dict(color='black'), flierprops=dict(marker=''))

    # Scatter points
    plt.scatter([i - 0.2] * len(sens), sens, color='blue', alpha=0.8, s=40, label='Sensitive' if i == 0 else "")
    plt.scatter([i + 0.2] * len(res), res, color='red', alpha=0.8, s=40, label='Resistant' if i == 0 else "")

# Aesthetic settings
plt.xticks(ticks=range(len(ordered_minicircles)), labels=ordered_minicircles, rotation=90, fontsize=30)
plt.yticks(fontsize=30)
plt.xlabel("Minicircle", fontsize=30)
plt.ylabel("Normalized Mean Coverage", fontsize=30)
plt.title("Significant Differences in Minicircle Coverage (p < 0.05)", fontsize=30)

plt.legend(title='Group', fontsize=12, title_fontsize=13)
sns.despine()
plt.grid(False)

# Save
plt.tight_layout()
plt.savefig("significant_minicircles_custom_order.pdf", format='pdf', bbox_inches='tight')
plt.savefig("significant_minicircles_custom_order.png", format='png', dpi=300, bbox_inches='tight')