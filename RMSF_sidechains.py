import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

trajectory = md.load('t.dcd', top='t.pdb')
atom_indices = trajectory.top.select('sidechain and resid 12 to 328')
rmsf = md.rmsf(trajectory, atom_indices=atom_indices)
plt.plot(rmsf, label='RMSF')
plt.xlabel('Residue Index')
plt.ylabel('RMSF (nm)')
plt.title('Root Mean Square Fluctuation (RMSF) for Side-Chain Atoms (Residues 12 to 328)')
plt.legend()
plt.show()
