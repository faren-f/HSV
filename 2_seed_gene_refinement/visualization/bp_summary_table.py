import pandas as pd
import os

# Create a dictionary with the data from the log output
bp_data = {
    'Biological Process': [
        'Egress and Envelopment',
        'Entry and Uncoating',
        'Uncharacterized',
        'Replication and Transcription',
        'Assembly and Packaging',
        'Immune Evasion',
        'Total'
    ],
    'Seed Genes': [19, 6, 6, 19, 11, 6, 67],  # Sum of all seed genes
    'PPIs': [236, 69, 27, 190, 51, 204, 777],  # Sum of all PPIs
    'Human UniProt IDs': [205, 63, 23, 162, 44, 168, 579]  # From background control
}

# Create DataFrame
df = pd.DataFrame(bp_data)

# Set the Biological Process as the index for better display
df.set_index('Biological Process', inplace=True)

# Display the table
print("\nSummary Table of HSV1-2 Seed Genes by Biological Process Category\n")
print(df)

# Save the table to a CSV file
output_dir = os.path.dirname(os.path.abspath(__file__))
output_file = os.path.join(output_dir, 'bp_summary_table.csv')
df.to_csv(output_file)
print(f"\nTable saved to: {output_file}")

# Create a more visually appealing HTML table
html_file = os.path.join(output_dir, 'bp_summary_table.html')
html_content = """
<!DOCTYPE html>
<html>
<head>
    <title>HSV1-2 Seed Genes Summary</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        table { border-collapse: collapse; width: 80%; margin: 20px auto; }
        th, td { padding: 8px 12px; text-align: center; border: 1px solid #ddd; }
        th { background-color: #4CAF50; color: white; }
        tr:nth-child(even) { background-color: #f2f2f2; }
        tr:hover { background-color: #ddd; }
        tr:last-child { font-weight: bold; background-color: #e6f7ff; }
        h2 { text-align: center; color: #333; }
    </style>
</head>
<body>
    <h2>Summary Table of HSV1-2 Seed Genes by Biological Process Category</h2>
    <table>
        <tr>
            <th>Biological Process</th>
            <th>Seed Genes</th>
            <th>PPIs</th>
            <th>Human UniProt IDs</th>
        </tr>
"""

# Add rows for each BP category
for i, bp in enumerate(bp_data['Biological Process']):
    row_class = ' class="total-row"' if bp == 'Total' else ''
    html_content += f"""
        <tr{row_class}>
            <td>{bp}</td>
            <td>{bp_data['Seed Genes'][i]}</td>
            <td>{bp_data['PPIs'][i]}</td>
            <td>{bp_data['Human UniProt IDs'][i]}</td>
        </tr>"""

html_content += """
    </table>
</body>
</html>
"""

with open(html_file, 'w') as f:
    f.write(html_content)

print(f"HTML table saved to: {html_file}")