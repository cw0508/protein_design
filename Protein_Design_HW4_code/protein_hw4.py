# ---------------------------------------------------------------------
import os
from Bio.PDB import PDBParser, DSSP, PDBIO, Select
import matplotlib.pyplot as plt
import numpy as np

# Class to clean a PDB file
class CleanATOM(Select):
    def accept_atom(self, atom):
        return atom.is_disordered() == 0  # Remove disordered atoms

def clean_pdb(input_pdb, output_pdb):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", input_pdb)
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb, select=CleanATOM())

def calculate_alpha_fraction(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    model = structure[0]  # Assuming first model
    dssp = DSSP(model, pdb_file)

    alpha_helix_count = sum(1 for residue in dssp if residue[2] == 'H')
    total_residues = len(dssp)
    return alpha_helix_count / total_residues if total_residues > 0 else 0

# Input and output directories
pdb_dir = "./"  # Input PDB files directory
output_dir = "./output_results"  # Output directory for results

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Process PDB files
pdb_files = [f for f in os.listdir(pdb_dir) if f.endswith(".pdb")]
alpha_fractions = []

# File to save results
results_file = os.path.join(output_dir, "alpha_helix_fractions.txt")

with open(results_file, "w") as f:
    for pdb_file in pdb_files:
        try:
            # Clean the PDB file
            clean_file = os.path.join(output_dir, pdb_file.replace(".pdb", "_clean.pdb"))
            clean_pdb(os.path.join(pdb_dir, pdb_file), clean_file)

            # Run DSSP on the cleaned file
            alpha_fraction = calculate_alpha_fraction(clean_file)
            alpha_fractions.append(alpha_fraction)

            # Write results to file
            f.write(f"{pdb_file}: {alpha_fraction:.4f}\n")
        except Exception as e:
            print(f"Error processing {pdb_file}: {e}")

print("Alpha Helix Fractions:", alpha_fractions)

# Save the histogram plot
plt.hist(alpha_fractions, bins=np.linspace(0, 1, 20), edgecolor="black")
plt.title("Histogram of Alpha Helix Fractions per Design")
plt.xlabel("Fraction of Alpha Helices")
plt.ylabel("Frequency")
histogram_path = os.path.join(output_dir, "alpha_helix_histogram.png")
plt.savefig(histogram_path)
plt.show()



# ---------------------------------------------------------------------
from Bio.PDB import PDBParser, PPBuilder

# Load the PDB file
pdb_file = "/users/chessbunny/desktop/2d1s_sa_5.pdb" 
parser = PDBParser()
structure = parser.get_structure("protein", pdb_file)

# Extract the amino acid sequence
ppb = PPBuilder()
for i, pp in enumerate(ppb.build_peptides(structure)):
    print(f"Sequence {i + 1}: {pp.get_sequence()}")



# ---------------------------------------------------------------------
from pyrosetta import *
from pyrosetta.teaching import *

# Initialize PyRosetta (no extra options if ligand params are not needed)
init()

# Load the complex PDB file
pose1 = pose_from_pdb("ligand_protein_complex1.pdb")
pose2 = pose_from_pdb("ligand_protein_complex2.pdb")
pose3 = pose_from_pdb("ligand_protein_complex3.pdb")

# Get the standard Rosetta scoring function
scorefxn = get_fa_scorefxn()

# Calculate the total energy of the complex
total_energy1 = scorefxn(pose1)
total_energy2 = scorefxn(pose2)
total_energy3 = scorefxn(pose3)

# Print the results
print("Total energy of the complex1:", total_energy1)
print("Total energy of the complex1:", total_energy2)
print("Total energy of the complex1:", total_energy3)

