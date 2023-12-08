import mdtraj as md
import numpy as np

# Load trajectory and topology
trajectory = md.load('trajectory.dcd', top='topology.pdb')

# Define the ligand selection using segment IDs
ligand_selection = 'segid PROB'

# Select protein atoms using segment ID
protein_selection = 'segid PROA'

# Calculate contacts between GPCR and ligand
contacts = md.compute_contacts(trajectory, contacts=[(ligand_selection, protein_selection)])

# Print the percentage of frames where the ligand interacts with the GPCR
interaction_percentage = np.mean(contacts[:, 0] < 0.4) * 100
print(f'Ligand interacts with GPCR in {interaction_percentage:.2f}% of frames.')
