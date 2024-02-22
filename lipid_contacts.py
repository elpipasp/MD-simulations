import mdtraj as md

trajectory = md.load('t.dcd', top='t.pdb')
lipid_selection = 'resname POPC'
lipid_contacts = md.compute_contacts(trajectory, contacts=[(lipid_selection, 'protein')])
interaction_percentage = np.mean(lipid_contacts[:, 0] < 0.4) * 100
print(f'GPCR interacts with lipids in {interaction_percentage:.2f}% of frames.')
