import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

# Step 2: Load Trajectory and Select Atom Indices
trajectory = md.load('trajectory.dcd', top='topology.pdb')
atom_indices = trajectory.top.select('alpha and resid 1 to 100')

# Step 3: Prepare Data
data = trajectory.xyz[:, atom_indices, :]

# Step 4: Data Standardization
data_mean = data.mean(axis=0)
data_std = data.std(axis=0)
data_standardized = (data - data_mean) / data_std

# Step 5: Perform PCA
pca = md.decompose.principal_component_analysis(data_standardized)

# Step 6: Access PCA Results
eigenvalues = pca.eigenvalues
eigenvectors = pca.eigenvectors

# Step 7: Visualize Explained Variance
explained_variance_ratio = eigenvalues / eigenvalues.sum()
cumulative_variance_ratio = explained_variance_ratio.cumsum()

plt.plot(cumulative_variance_ratio)
plt.xlabel('Number of Principal Components')
plt.ylabel('Cumulative Explained Variance')
plt.show()

# Step 8: Project Trajectory onto Principal Components
projected_data = pca.transform(data_standardized)

