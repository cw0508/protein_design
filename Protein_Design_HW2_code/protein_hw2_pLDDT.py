import json
import matplotlib.pyplot as plt
import numpy as np

def extract_plddt_from_json(json_file):
    # Load the JSON file
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    # Extract pLDDT scores from the 'atom_plddts' field
    plddt_scores = data.get('atom_plddts', [])
    
    if not plddt_scores:
        print(f"No pLDDT scores found in {json_file}.")
        return None
    
    return plddt_scores

def plot_plddt(all_scores, labels):
    # Plot the pLDDT scores for multiple data sets
    plt.figure(figsize=(12, 8))
    
    for scores, label in zip(all_scores, labels):
        plt.plot(np.arange(len(scores)), scores, marker='o', linestyle='-', markersize=3, label=label)

    plt.title('pLDDT Scores across Residues (Multiple Files)', fontsize=14)
    plt.xlabel('Residue Index', fontsize=12)
    plt.ylabel('pLDDT Score', fontsize=12)

    # Add confidence threshold lines
    plt.axhline(y=90, color='green', linestyle='--', label='Very high confidence (>90)')
    plt.axhline(y=70, color='orange', linestyle='--', label='Medium confidence (70-90)')
    plt.axhline(y=50, color='red', linestyle='--', label='Low confidence (<70)')

    plt.legend()
    plt.grid(True)
    plt.show()

# List of file paths
json_files = [
    'fold_2024_10_27_16_37/fold_2024_10_27_16_37_full_data_0.json',
    'fold_2024_10_27_16_37/fold_2024_10_27_16_37_full_data_1.json',
    'fold_2024_10_27_16_37/fold_2024_10_27_16_37_full_data_2.json',
    'fold_2024_10_27_16_37/fold_2024_10_27_16_37_full_data_3.json'
]

json_files = [
    'fold_2024_10_27_16_37-2/fold_2024_10_27_16_37_full_data_0.json',
    'fold_2024_10_27_16_37-2/fold_2024_10_27_16_37_full_data_1.json',
    'fold_2024_10_27_16_37-2/fold_2024_10_27_16_37_full_data_2.json',
    'fold_2024_10_27_16_37-2/fold_2024_10_27_16_37_full_data_3.json'
]

json_files = [
    'fold_2024_10_27_16_37-3/fold_2024_10_27_16_37_full_data_0.json',
    'fold_2024_10_27_16_37-3/fold_2024_10_27_16_37_full_data_1.json',
    'fold_2024_10_27_16_37-3/fold_2024_10_27_16_37_full_data_2.json',
    'fold_2024_10_27_16_37-3/fold_2024_10_27_16_37_full_data_3.json'
]

json_files = [
    'fold_2024_10_27_16_37-4/fold_2024_10_27_16_37_full_data_0.json',
    'fold_2024_10_27_16_37-4/fold_2024_10_27_16_37_full_data_1.json',
    'fold_2024_10_27_16_37-4/fold_2024_10_27_16_37_full_data_2.json',
    'fold_2024_10_27_16_37-4/fold_2024_10_27_16_37_full_data_3.json'
]

all_plddt_scores = []
labels = ['Dataset 0', 'Dataset 1', 'Dataset 2', 'Dataset 3']

# Extract pLDDT scores for all files
for json_file in json_files:
    plddt_scores = extract_plddt_from_json(json_file)
    if plddt_scores:
        all_plddt_scores.append(plddt_scores)

# Plot the results
if all_plddt_scores:
    plot_plddt(all_plddt_scores, labels)