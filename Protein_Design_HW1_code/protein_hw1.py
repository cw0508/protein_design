import os 
import warnings
from Bio import PDB
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from prody import *
import time
import numpy as np
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser


# Suppress PDB warnings generically
warnings.filterwarnings("ignore", message="WARNING:.*discontinuous.*")

#-----------------------------------------------------------------------------------
def get_pdb_id(pdb_file):
    """
    Extracts and returns the PDB ID (filename without the extension) from a PDB file path.
    """
    return os.path.splitext(os.path.basename(pdb_file))[0]

#-----------------------------------------------------------------------------------
def extract_sequence_from_pdb(pdb_file):
    """
    Extracts the sequence of the first chain from the given PDB file using PPBuilder.
    """
    parser = PDB.PDBParser()
    structure = parser.get_structure("protein", pdb_file)
    
    # Use the PPBuilder to extract the sequence
    ppb = PDB.PPBuilder()
    
    # Get the sequence from the first polypeptide found in the first model
    for pp in ppb.build_peptides(structure):
        sequence = pp.get_sequence() 
        return str(sequence)

#-----------------------------------------------------------------------------------
def sequence_distance(pdb1, pdb2):
    """
    Calculates sequence distance (based on percentage of non-identical residues)
    between two proteins based on their PDB files.
    """
    seq1 = extract_sequence_from_pdb(pdb1)
    seq2 = extract_sequence_from_pdb(pdb2)

    # Use global pairwise alignment (Needleman-Wunsch) to align sequences
    alignments = pairwise2.align.globalxx(seq1, seq2)
    
    # Get the alignment with the highest score (optimal alignment)
    best_alignment = alignments[0]
    
    # Calculate the percentage of identical residues in the alignment
    identical_residues = sum(1 for a, b in zip(best_alignment[0], best_alignment[1]) if a == b)
    alignment_length = len(best_alignment[0])
    
    # Calculate percentage identity
    percentage_identity = (identical_residues / alignment_length) * 100
    return 100 - percentage_identity 

#-----------------------------------------------------------------------------------
def calculate_rmsd(pdb1, pdb2):
    """
    Calculates the RMSD (Root Mean Square Deviation) between two PDB files based on their backbone alignment.
    """
    # Parse the PDB files
    structure1 = parsePDB(pdb1)
    structure2 = parsePDB(pdb2)

    # Select backbone atoms for alignment
    backbone1 = structure1.select('backbone')
    backbone2 = structure2.select('backbone')

    if backbone1 is None or backbone2 is None:
        print(f"Error: Backbone atoms not found in {pdb1} or {pdb2}")
        return None

    # Check if both backbones have the same number of atoms
    if backbone1.numAtoms() != backbone2.numAtoms():
        print(f"Skipping RMSD calculation for {pdb1} and {pdb2} due to mismatched atom counts.")
        return None

    # Superimpose the two structures
    transformation = calcTransformation(backbone2, backbone1)
    transformation.apply(backbone2)

    # Calculate RMSD
    rmsd_value = calcRMSD(backbone1, backbone2)
    return rmsd_value

#-----------------------------------------------------------------------------------
def benchmark_proteins(pdb_set_1, pdb_set_2):
    """
    Benchmarks protein pairs for both sequence distance and structure distance (RMSD).
    """
    for pdb1 in pdb_set_1:
        for pdb2 in pdb_set_2:
            start_time = time.time()  
            
            try:
                seq_distance = sequence_distance(pdb1, pdb2)
                struct_distance = calculate_rmsd(pdb1, pdb2)
                
                pdb_id1 = get_pdb_id(pdb1) 
                pdb_id2 = get_pdb_id(pdb2)  

                print(f"Comparing {pdb_id1} and {pdb_id2}:")
                print(f"  Sequence distance: {seq_distance}")
                
                if struct_distance is not None:
                    print(f"  Structure distance (RMSD): {struct_distance}")
                else:
                    print(f"  Structure distance (RMSD): Skipped due to mismatched atom counts or missing backbone.")
                    
            except Exception as e:
                print(f"Error comparing {pdb1} and {pdb2}: {e}")
                
            finally:
                print(f"Time taken: {time.time() - start_time:.2f} seconds\n")


#-----------------------------------------------------------------------------------
# load the files: 
sequence_similar_set = [
    "2ptn.pdb",
    "1aq7.pdb",
    "1auj.pdb"
]

structure_similar_set = [
    "2ptn.pdb",
    "3ptb.pdb",
    "5mor.pdb"
]

random_pdb_ids = [
    "1xyz.pdb",
    "4hhb.pdb",
    "3def.pdb"
]

#-----------------------------------------------------------------------------------
print("------------------- START --------------------")
# Perform the benchmark for sequence-similar, structure-similar, and random sets
print("Benchmarking Sequence-Similar Set:")
benchmark_proteins(sequence_similar_set, sequence_similar_set)
print("----------------------------------------------")

print("\nBenchmarking Structure-Similar Set:")
benchmark_proteins(structure_similar_set, structure_similar_set)
print("----------------------------------------------")

print("\nBenchmarking Random Proteins Set:")
benchmark_proteins(random_pdb_ids, random_pdb_ids)
print("-------------------- END ----------------------")



#-----------------------------------------------------------------------------------
def compute_distance_map(pdb_file):
    """
    Computes the intra-protein residue-wise distance map for the given PDB file.
    The distances are calculated between the alpha carbon (Cα) atoms of all residues.
    """
    parser = PDBParser()
    structure = parser.get_structure("protein", pdb_file)
    
    # Extract alpha carbon (Cα) atoms
    ca_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:  # Check if residue has an alpha carbon (Cα)
                    ca_atoms.append(residue['CA'].get_coord())
    
    # Convert the list of Cα atom coordinates into a NumPy array
    ca_coords = np.array(ca_atoms)
    
    # Initialize an empty distance matrix
    n = len(ca_coords)
    distance_map = np.zeros((n, n))
    
    # Compute pairwise distances
    for i in range(n):
        for j in range(i, n):
            distance = np.linalg.norm(ca_coords[i] - ca_coords[j])
            distance_map[i, j] = distance_map[j, i] = distance  
    
    return distance_map

def plot_distance_map(distance_map, title="Intra-Protein Residue-Wise Distance Map"):
    """
    Plots the distance map using Matplotlib.
    """
    plt.imshow(distance_map, cmap='viridis')
    plt.colorbar(label='Distance (Å)')
    plt.title(title)
    plt.xlabel("Residue Index")
    plt.ylabel("Residue Index")
    plt.show()


pdb_file = "2ptn.pdb" 
distance_map = compute_distance_map(pdb_file)
plot_distance_map(distance_map)
