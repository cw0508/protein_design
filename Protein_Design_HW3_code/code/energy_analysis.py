import pyrosetta
from pyrosetta import pose_from_pdb, get_fa_scorefxn
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover

# Initialize PyRosetta
pyrosetta.init()

# Load the binder poses
pose1 = pose_from_pdb("nca_bindcraft_1.pdb")  # First binder
pose2 = pose_from_pdb("chai_1.pdb")  # Second binder

# Load the target protein
target_pose = pose_from_pdb("2wh6.pdb")  # Replace with the actual target protein PDB file

# Use the full-atom scoring function
scorefxn = get_fa_scorefxn()

# Compute individual energy scores (stability)
energy1 = scorefxn(pose1)
energy2 = scorefxn(pose2)

# Print individual binder stability scores
print(f"Energy Score for Bindcraft Binder: {energy1}")
print(f"Energy Score for Chai-1 Binder: {energy2}")

# Combine target and binders into single poses
combined_pose1 = target_pose.clone()
combined_pose1.append_pose_by_jump(pose1, 1)  # Append binder1 to target protein

combined_pose2 = target_pose.clone()
combined_pose2.append_pose_by_jump(pose2, 1)  # Append binder2 to target protein

# Save combined poses for verification
combined_pose1.dump_pdb("combined_pose1.pdb")
combined_pose2.dump_pdb("combined_pose2.pdb")

# Analyze protein-protein interactions
interface_analyzer1 = InterfaceAnalyzerMover(1)  # Assume binder1 is in Chain B
interface_analyzer1.apply(combined_pose1)
binding_energy1 = interface_analyzer1.get_interface_dG()
binding_sasa1 = interface_analyzer1.get_interface_delta_sasa()

interface_analyzer2 = InterfaceAnalyzerMover(1)  # Assume binder2 is in Chain B
interface_analyzer2.apply(combined_pose2)
binding_energy2 = interface_analyzer2.get_interface_dG()
binding_sasa2 = interface_analyzer2.get_interface_delta_sasa()

# Use the full-atom scoring function
scorefxn = get_fa_scorefxn()

# Compute energy scores
energy1 = scorefxn(pose1)
energy2 = scorefxn(pose2)

# Print results
print(f"Bindcraft Binder Energy Score: {energy1}")
print(f"Chai-1 Energy Score Binder: {energy2}")

# Print interaction metrics
print(f"Bindcraft Binder Interaction Energy (ΔG): {binding_energy1} REU")
print(f"Bindcraft Binder Buried SASA (ΔSASA): {binding_sasa1} Å²")

print(f"Chai-1 Binder Interaction Energy (ΔG): {binding_energy2} REU")
print(f"Chai-1 Binder Buried SASA (ΔSASA): {binding_sasa2} Å²")

