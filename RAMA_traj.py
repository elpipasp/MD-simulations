import MDAnalysis as mda
from MDAnalysis.analysis.dihedrals import Ramachandran
import matplotlib.pyplot as plt

psf_file = 'ROS_1_SC.psf'
dcd_file = 'D22.dcd'
u = mda.Universe(psf_file, dcd_file)
selection = u.select_atoms("segid PROA")
R = Ramachandran(selection).run()
fig, ax = plt.subplots(figsize=(10, 8))
R.plot(ax=ax, color='k', marker='o', alpha=0.3, ref=True) 
ax.set_xlabel('Phi (°)', fontsize=24, weight='bold')
ax.set_ylabel('Psi (°)', fontsize=24, weight='bold')
ax.set_title('GnRH1R in ROS-1', fontsize=22, weight='bold')
ax.tick_params(axis='both', which='major', labelsize=22)
plt.grid(True)
plt.tight_layout()
plt.show()
