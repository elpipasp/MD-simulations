import MDAnalysis as mda
from MDAnalysis.analysis.dihedrals import Ramachandran
import matplotlib.pyplot as plt

pdb_file = '2ROS2_600ns.pdb'
u = mda.Universe(pdb_file)
r = u.select_atoms("protein")
R = Ramachandran(r).run()

fig, ax = plt.subplots(figsize=(10, 8)) 
R.plot(ax=ax, color='k', marker='o', ref=True)
ax.set_xlabel('Psi (°)', fontsize=14, weight='bold')
ax.set_ylabel('Phi (°)', fontsize=14, weight='bold')
ax.set_title('ROS-2 (2)', fontsize=16, weight='bold')
ax.tick_params(axis='both', which='major', labelsize=14)
plt.grid(True)
plt.show()
