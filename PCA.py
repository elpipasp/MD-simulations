import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

# Step 1: Define a list of trajectory file names
trajectory_files = ['D1.dcd', 'D2.dcd', 'D3.dcd', 'D4.dcd']  # Add more file names as needed

# Step 2: Load Trajectories and Select Atom Indices
all_data = []

for traj_file in trajectory_files:
    trajectory = md.load(traj_file, top='No_Nter_Gnrh_1_SC.pdb')
    atom_indices = trajectory.top.select("name CA and resid 12 to 328")
    data = trajectory.xyz[:, atom_indices, :]
    all_data.append(data)

# Step 3: Concatenate Data from Multiple Trajectories
all_data = np.concatenate(all_data, axis=0)

# Step 4: Data Standardization
data_mean = all_data.mean(axis=0)
data_std = all_data.std(axis=0)
data_standardized = (all_data - data_mean) / data_std

# Step 5: Perform PCA using scikit-learn
num_components = min(all_data.shape[1], all_data.shape[2])
pca = PCA(n_components=num_components)
projected_data = pca.fit_transform(data_standardized.reshape(data_standardized.shape[0], -1))

# Step 6: Visualize Explained Variance
explained_variance_ratio = pca.explained_variance_ratio_
cumulative_variance_ratio = explained_variance_ratio.cumsum()

plt.plot(cumulative_variance_ratio)
plt.xlabel('Number of Principal Components')
plt.ylabel('Cumulative Explained Variance')
plt.show()

