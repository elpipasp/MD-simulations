import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
import matplotlib.pyplot as plt
import warnings
import numpy as np

warnings.filterwarnings('ignore')

adk_topology = 'NTER_FREE_SC.pdb'
adk_trajectory = 'D21.dcd'

u = mda.Universe(adk_topology, adk_trajectory)

protein_sidechain = u.select_atoms('protein and ((resname GLY and name CA) or (not resname GLY and name CB))')

average = align.AverageStructure(u, select='protein and name CA and (resid 1:328)', ref_frame=0).run()
ref = average.universe
aligner = align.AlignTraj(u, ref, select='protein and name CA and (resid 1:328)', in_memory=True).run()

R = rms.RMSF(protein_sidechain).run()
data = np.column_stack((protein_sidechain.resids, R.rmsf))
np.savetxt('rmsf_CB_1inactive.txt', data, header='Residue_number RMSF', fmt='%d %.4f') 

plt.plot(protein_sidechain.resids, R.rmsf)
plt.xlabel('Residue number', fontsize=14, weight='bold')
plt.ylabel('RMSF ($\AA$)', fontsize=14, weight='bold')
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
plt.xlim(12, 328)
plt.xticks([12, 33, 67, 74, 93, 110, 145, 152, 177, 204, 244, 257, 293, 302, 328], rotation=90)
midpoints = [(33 + 67) / 2, (74 + 93) / 2, (110 + 145) / 2, (152 + 177) / 2, (204 + 244) / 2, (257 + 293) / 2, (302 + 327) / 2]

for midpoint, label in zip(midpoints, ['TM1', 'TM2', 'TM3', 'TM4', 'TM5', 'TM6', 'TM7']):
    plt.text(midpoint, 4.5, label, ha='center', va='center', color='black', fontsize=14, weight='bold')

plt.axvspan(33, 67, zorder=0, alpha=0.2, color='grey', label='TM1')
plt.axvspan(74, 93, zorder=0, alpha=0.2, color='grey', label='TM2')
plt.axvspan(110, 145, zorder=0, alpha=0.2, color='grey', label='TM3')
plt.axvspan(152, 177, zorder=0, alpha=0.2, color='grey', label='TM4')
plt.axvspan(204, 244, zorder=0, alpha=0.2, color='grey', label='TM5')
plt.axvspan(257, 293, zorder=0, alpha=0.2, color='grey', label='TM6')
plt.axvspan(302, 327, zorder=0, alpha=0.2, color='grey', label='TM7')
plt.legend(fontsize=12)
plt.show()
