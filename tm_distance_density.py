import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

dcd_files = ['D3.dcd','D19.dcd', 'D20.dcd', 'D21.dcd', 'D23.dcd']
psf_file = 'ROS_1_SC.psf'
u = mda.Universe(psf_file)

#define residues and atoms
tm3_selection = u.select_atoms("segid PROA and resid 139 and name CA")
tm6_selection = u.select_atoms("segid PROA and resid 265 and name CA")
tm7_selection = u.select_atoms("segid PROA and resid 323 and name CA")

tm3_tm6_distances = []
tm3_tm7_distances = []

for dcd_file in dcd_files:
    u.load_new(dcd_file)
    
    for ts in u.trajectory:
        distance_tm3_tm6 = np.linalg.norm(tm3_selection.positions - tm6_selection.positions)
        tm3_tm6_distances.append(distance_tm3_tm6)
        distance_tm3_tm7 = np.linalg.norm(tm3_selection.positions - tm7_selection.positions)
        tm3_tm7_distances.append(distance_tm3_tm7)

sns.kdeplot(x=tm3_tm7_distances, y=tm3_tm6_distances, cmap='viridis', fill=True, thresh=0.1, levels=100, center=0, cbar=True, cbar_kws={'ticks':[], 'orientation': 'vertical'})

#AF active for comparison
pdb_file = 'ClassA_gnrhr_human_Active_AF_2022-08-16_GPCRdb.pdb'
pdb_data = mda.Universe(pdb_file)

#extract coor of AF active receptor
pdb_tm3_position = pdb_data.select_atoms("segid A and resid 139 and name CA").positions
pdb_tm6_position = pdb_data.select_atoms("segid A and resid 265 and name CA").positions
pdb_tm7_position = pdb_data.select_atoms("segid A and resid 323 and name CA").positions

#calculate distances
pdb_tm3_tm6_distance = np.linalg.norm(pdb_tm3_position - pdb_tm6_position)
pdb_tm3_tm7_distance = np.linalg.norm(pdb_tm3_position - pdb_tm7_position)

print("TM3-TM6 Distance:", pdb_tm3_tm6_distance)
print("TM3-TM7 Distance:", pdb_tm3_tm7_distance)

plt.scatter(pdb_tm3_tm7_distance, pdb_tm3_tm6_distance, color='red', s=20, label='AF-Active GnRH1R')
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.xlabel('TM3-TM7 (Å)', fontsize=12, fontweight='bold')
plt.ylabel('TM3-TM6 (Å)', fontsize=12, fontweight='bold')
plt.title('ROS-1 (1)', fontsize=14, fontweight='bold')
plt.legend(loc='upper right', fontsize=10)
plt.show()
