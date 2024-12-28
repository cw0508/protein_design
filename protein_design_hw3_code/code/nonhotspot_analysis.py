from pyrosetta import init, pose_from_pdb, get_fa_scorefxn
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.core.scoring import CA_rmsd
import pandas as pd
import os


def calculate_rmsd(reference_pdb_path, target_pdb_path):
    """
    Calculate RMSD between two PDB files using PyRosetta.
    """
    ref_pose = pose_from_pdb(reference_pdb_path)
    target_pose = pose_from_pdb(target_pdb_path)
    return CA_rmsd(ref_pose, target_pose)


def calculate_plddt_approximation(pose):
    """
    Approximate pLDDT using residue-level energy scores from PyRosetta.
    Ensures energies are calculated before accessing them.
    """
    score_function = get_fa_scorefxn()
    score_function(pose)  # Ensure energies are updated
    residue_scores = []

    for i in range(1, pose.total_residue() + 1):  # Residue indices in PyRosetta are 1-based
        residue_scores.append(pose.energies().residue_total_energy(i))

    # Normalize scores to approximate confidence
    normalized_scores = [1 / (1 + abs(score)) for score in residue_scores]
    avg_score = sum(normalized_scores) / len(normalized_scores)
    return avg_score


def analyze_binder(binder_pdb_path, target_pdb_path, output_pdb_path):
    """
    Analyzes the binding affinity between a binder and a target protein.
    """
    try:
        # Load the binder and target PDB structures
        binder_pose = pose_from_pdb(binder_pdb_path)
        target_pose = pose_from_pdb(target_pdb_path)
        
        # Combine binder and target into a single pose
        combined_pose = binder_pose.clone()
        combined_pose.append_pose_by_jump(target_pose, 1)

        # Save combined pose as a PDB for debugging/visualization
        combined_pose.dump_pdb(output_pdb_path)

        # Use the full-atom scoring function for affinity analysis
        score_function = get_fa_scorefxn()

        # Interface analysis to estimate binding energy and affinity
        interface_analyzer = InterfaceAnalyzerMover(1)
        interface_analyzer.apply(combined_pose)

        # Retrieve interface metrics
        binding_energy = interface_analyzer.get_interface_dG()
        binding_sasa = interface_analyzer.get_interface_delta_sasa()

        # RMSD between the binder and target PDBs
        rmsd_value = calculate_rmsd(binder_pdb_path, target_pdb_path)

        # Approximate pLDDT for the combined pose
        plddt_value = calculate_plddt_approximation(combined_pose)

        # Return results as a dictionary
        return {
            "binder_pdb": os.path.basename(binder_pdb_path),
            "binding_energy": binding_energy,
            "binding_sasa": binding_sasa,
            "rmsd": rmsd_value,
            "approx_plddt": plddt_value
        }

    except Exception as e:
        print(f"Error analyzing {binder_pdb_path}: {str(e)}")
        return {
            "binder_pdb": os.path.basename(binder_pdb_path),
            "binding_energy": None,
            "binding_sasa": None,
            "rmsd": None,
            "approx_plddt": None,
            "error": str(e)
        }


def analyze_multiple_binders(binders, target_pdb, output_dir="output"):
    """
    Analyze multiple binders and return results as a Pandas DataFrame.
    """
    results = []
    for binder_pdb in binders:
        output_pdb_path = os.path.join(output_dir, f"complex_{os.path.basename(binder_pdb)}")
        result = analyze_binder(binder_pdb, target_pdb, output_pdb_path)
        results.append(result)

    # Convert results to a Pandas DataFrame
    return pd.DataFrame(results)


# Main script
if __name__ == "__main__":
    # Initialize PyRosetta
    init("-mute all")

    # Define the list of binder PDB files
    binders = [
        "nca_nhs.pdb",
        "nca_nhs_m.pdb"
    ]

    # Target protein PDB file
    target_pdb = "2wh6.pdb"

    # Output directory
    output_dir = "output"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Analyze all binders
    results_df = analyze_multiple_binders(binders, target_pdb, output_dir)

    # Display the DataFrame
    print(results_df)

    # Save results to a CSV file
    output_csv = os.path.join(output_dir, "binding_analysis_results_nonhotspot.csv")
    results_df.to_csv(output_csv, index=False)
    print(f"Results saved to {output_csv}")
