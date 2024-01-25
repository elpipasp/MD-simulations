import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy.stats import gaussian_kde
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Step 1: Load MD trajectory
traj_files = ['D5.dcd', 'D6.dcd', 'D7.dcd', 'D8.dcd', 'D9.dcd', 'D10.dcd', 'D11.dcd', 'D12.dcd', 'D13.dcd']
traj = md.load(traj_files, top='ROS_1_SC.pdb')

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

# Plot the density plot
fig, ax = plt.subplots()
im = ax.pcolormesh(x_grid, y_grid, z.reshape(x_grid.shape), shading='auto', cmap='viridis')
ax.set_xlabel('PCA 1', fontsize=14, fontweight='bold')
ax.set_ylabel('PCA 2', fontsize=14, fontweight='bold')
ax.tick_params(axis='both', which='both', length=0)

# Create colorbar
divider = make_axes_locatable(ax)
cax = divider.append_axes("bottom", size="5%", pad=0.5)
cbar = plt.colorbar(im, cax=cax, orientation='horizontal', ticks=[z.min(), z.max()])
cbar.ax.set_xticklabels(['Low', 'High'], fontsize=10)

# Set Population Density label in bold and adjust pad
cbar.ax.set_xlabel('Population Density', fontsize=10, fontweight='bold', labelpad=-10)

# Adjust the spacing above the suptitle
plt.subplots_adjust(top=0.92)

# Title
plt.suptitle('PCA ROS-1 (1)', fontsize=14, fontweight='bold')
plt.show()
