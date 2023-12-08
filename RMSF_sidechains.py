import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

# Load trajectory and topology
trajectory = md.load('trajectory.dcd', top='topology.pdb')

# Select side-chain atoms for residues 12 to 328
atom_indices = trajectory.top.select('sidechain and resid 12 to 328')

# Calculate RMSF
rmsf = md.rmsf(trajectory, atom_indices=atom_indices)

# Plot RMSF
plt.plot(rmsf, label='RMSF')
plt.xlabel('Residue Index')
plt.ylabel('RMSF (nm)')
plt.title('Root Mean Square Fluctuation (RMSF) for Side-Chain Atoms (Residues 12 to 328)')
plt.legend()
plt.show()
