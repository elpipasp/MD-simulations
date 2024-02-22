import mdtraj as md
import numpy as np

trajectory = md.load('t.dcd', top='t.pdb')
ligand_selection = 'segid PROB'
protein_selection = 'segid PROA'd
contacts = md.compute_contacts(trajectory, contacts=[(ligand_selection, protein_selection)])
interaction_percentage = np.mean(contacts[:, 0] < 0.4) * 100
print(f'Ligand interacts with GPCR in {interaction_percentage:.2f}% of frames.')
