import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# Load the data from the CSV file
output_dir = os.path.dirname(os.path.abspath(__file__))
csv_file = os.path.join(output_dir, 'bp_summary_table.csv')
df = pd.read_csv(csv_file)

# Set up the figure with a larger size for better visibility
plt.figure(figsize=(12, 8))

# Set width of bars
barWidth = 0.25

# Set positions of the bars on X axis
r1 = np.arange(len(df) - 1)  # Exclude the 'Total' row for the bar chart
r2 = [x + barWidth for x in r1]
r3 = [x + barWidth for x in r2]

# Filter out the 'Total' row for plotting
df_plot = df[df['Biological Process'] != 'Total']

# Create bars
plt.bar(r1, df_plot['Seed Genes'], width=barWidth, edgecolor='grey', label='Seed Genes')
plt.bar(r2, df_plot['PPIs'], width=barWidth, edgecolor='grey', label='PPIs')
plt.bar(r3, df_plot['Human UniProt IDs'], width=barWidth, edgecolor='grey', label='Human UniProt IDs')

# Add labels and title
plt.xlabel('Biological Process', fontweight='bold', fontsize=12)
plt.ylabel('Count', fontweight='bold', fontsize=12)
plt.title('HSV1-2 Seed Genes, PPIs, and Human UniProt IDs by Biological Process', fontweight='bold', fontsize=14)

# Add xticks on the middle of the group bars
plt.xticks([r + barWidth for r in range(len(df_plot))], df_plot['Biological Process'], rotation=45, ha='right')

# Create legend & Show graphic
plt.legend()
plt.tight_layout()

# Save the figure
output_file = os.path.join(output_dir, 'bp_summary_chart.png')
plt.savefig(output_file, dpi=300)
plt.close()

print(f"\nChart saved to: {output_file}")

# Create a stacked bar chart to show proportions
plt.figure(figsize=(12, 8))

# Calculate percentages for each category
df_percent = df_plot.copy()
df_percent['Total Interactions'] = df_percent['PPIs']
df_percent['Unique Human Proteins'] = df_percent['Human UniProt IDs'] / df_percent['PPIs'] * 100
df_percent['Redundant Interactions'] = 100 - df_percent['Unique Human Proteins']

# Create stacked bars
bars = plt.bar(df_percent['Biological Process'], df_percent['PPIs'], edgecolor='grey')

# Add data labels on top of bars
for bar in bars:
    height = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2., height + 5,
            f'{height:.0f}', ha='center', va='bottom', fontsize=10)

# Add labels and title
plt.xlabel('Biological Process', fontweight='bold', fontsize=12)
plt.ylabel('Number of PPIs', fontweight='bold', fontsize=12)
plt.title('Total Protein-Protein Interactions by Biological Process', fontweight='bold', fontsize=14)
plt.xticks(rotation=45, ha='right')

# Add a text annotation with the total numbers
total_seeds = df[df['Biological Process'] == 'Total']['Seed Genes'].values[0]
total_ppis = df[df['Biological Process'] == 'Total']['PPIs'].values[0]
total_human = df[df['Biological Process'] == 'Total']['Human UniProt IDs'].values[0]

plt.figtext(0.5, 0.01, 
           f'Total Seed Genes: {total_seeds} | Total PPIs: {total_ppis} | Total Human UniProt IDs: {total_human}',
           ha='center', fontsize=12, bbox={'facecolor':'lightgrey', 'alpha':0.5, 'pad':5})

plt.tight_layout(pad=3)

# Save the figure
output_file = os.path.join(output_dir, 'bp_ppi_chart.png')
plt.savefig(output_file, dpi=300)
plt.close()

print(f"PPI chart saved to: {output_file}")