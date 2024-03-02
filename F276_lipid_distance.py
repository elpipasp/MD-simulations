import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

trajectory_files = ["D19.dcd", "D20.dcd", "D21.dcd", "D22.dcd"] #inactive to active
protein_residue = "resname PHE and resid 276"
POPC_residue = "resname POPC and resid 243 and name C218"
distances = []

for traj_file in trajectory_files:
    u = mda.Universe("ROS_1_SC.pdb", traj_file)
    protein_atoms = u.select_atoms(protein_residue)
    POPC_atoms = u.select_atoms(POPC_residue)

    for ts in u.trajectory:
        #calculate the distance between the com of 276 and C212 atom
        distance = np.linalg.norm(protein_atoms.center_of_mass() - POPC_atoms.positions)
        distances.append(distance)
time_axis = np.arange(0, len(distances)) / 50.0

#plot
plt.figure(figsize=(10, 6))
plt.plot(time_axis, distances, label='Distance Protein Residue 276 - POPC Residue 243 C212', linestyle='-', linewidth=1, color='deepskyblue')
plt.xlabel('Time (ns)', fontsize=16, fontweight='bold')
plt.ylabel('Distance (Å)', fontsize=16, fontweight='bold')
plt.title('F276 CoM and Lipid aryl group distance', fontsize=16)
plt.xlim(0, 100)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid(True)
np.savetxt('distance_F276_lipid_2ROS1.txt', np.column_stack((time_axis, distances)), header='Time (ns) Distance (Å)', fmt='%1.3f %1.3f')
plt.show()
