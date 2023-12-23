import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy.stats import gaussian_kde

# Step 1: Load MD trajectory
traj_files = ['D1.dcd', 'D2.dcd', 'D3.dcd', 'D4.dcd', 'D5.dcd', 'D6.dcd', 'D7.dcd', 'D8.dcd', 'D9.dcd', 'D10.dcd']
traj = md.load(traj_files, top='No_Nter_Gnrh_1_SC.pdb')

# Step 2: Specify residue range for the peptide
peptide_residues = range(1, 10)  # Replace with the actual residue range

# Step 3: Extract peptide coordinates with a specific segment name
selection_string = f'segname PROB and (residue {peptide_residues[0]} to {peptide_residues[-1]})'
peptide_traj = traj.atom_slice(traj.topology.select(selection_string))

# Step 4: Perform PCA
peptide_coords = peptide_traj.xyz.reshape(peptide_traj.n_frames, -1)
pca = PCA(n_components=2)
peptide_pcs = pca.fit_transform(peptide_coords)

# Step 5: Create a 2D KDE
x = peptide_pcs[:, 0]
y = peptide_pcs[:, 1]
kde = gaussian_kde(np.vstack([x, y]))

# Step 6: Evaluate the KDE on a grid
x_grid, y_grid = np.mgrid[x.min():x.max():100j, y.min():y.max():100j]
z = kde(np.vstack([x_grid.ravel(), y_grid.ravel()]))

# Step 7: Plot the density plot
plt.pcolormesh(x_grid, y_grid, z.reshape(x_grid.shape), shading='auto', cmap='viridis')
plt.colorbar(label='Density')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title('Density Plot in PC Space')
plt.show()
