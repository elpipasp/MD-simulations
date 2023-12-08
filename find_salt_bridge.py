import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

# Load trajectory and topology
trajectory = md.load('your_trajectory.dcd', top='your_topology.pdb')

# Define positively charged and negatively charged residues
positive_residues = ['LYS', 'ARG']
negative_residues = ['ASP', 'GLU']

# Define the distance threshold for salt bridges (3.2 Å)
distance_threshold = 0.32  # Convert to nanometers (1 Å = 0.1 nm)

# Find salt bridges and calculate their duration
salt_bridge_frames = []
salt_bridge_durations = []

for frame in range(len(trajectory)):
    positive_atoms = trajectory.top.select(f'resid {":".join(positive_residues)} and charge > 0')
    negative_atoms = trajectory.top.select(f'resid {":".join(negative_residues)} and charge < 0')

    distances = md.compute_distances(trajectory[frame], list(zip(positive_atoms, negative_atoms)))[0]

    # Check for salt bridges
    salt_bridges = np.where(distances < distance_threshold)[0]

    if salt_bridges.size > 0:
        # Record frames and durations of salt bridges
        for bridge in salt_bridges:
            salt_bridge_frames.append(frame)
            duration = 1
            while bridge + duration in salt_bridges:
                duration += 1
            salt_bridge_durations.append(duration)

# Plot salt bridge duration over the trajectory
plt.plot(salt_bridge_frames, salt_bridge_durations, marker='o', linestyle='-', color='b')
plt.xlabel('Frame')
plt.ylabel('Duration')
plt.title('Salt Bridge Duration Over Trajectory')
plt.show()
