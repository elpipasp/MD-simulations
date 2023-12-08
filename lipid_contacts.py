import mdtraj as md

# Load trajectory and topology
trajectory = md.load('trajectory.dcd', top='topology.pdb')

# Select lipid atoms in the membrane
lipid_selection = 'resname POPC'

# Calculate contacts between GPCR and lipids
lipid_contacts = md.compute_contacts(trajectory, contacts=[(lipid_selection, 'protein')])

# Print the percentage of frames where the GPCR interacts with lipids
interaction_percentage = np.mean(lipid_contacts[:, 0] < 0.4) * 100
print(f'GPCR interacts with lipids in {interaction_percentage:.2f}% of frames.')
