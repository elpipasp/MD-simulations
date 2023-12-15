import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
import matplotlib.pyplot as plt
import warnings

# suppress some MDAnalysis warnings about writing PDB files
warnings.filterwarnings('ignore')

# Replace these placeholders with the actual paths or filenames
adk_topology = 'No_Nter_Gnrh_1_SC.pdb'
adk_trajectory = 'D1.dcd'

u = mda.Universe(adk_topology, adk_trajectory)

# Define the selection for PROA (residues 12 to 328)
protein_proa = u.select_atoms('protein and name CA and (resid 12:328)')

average = align.AverageStructure(u, select='protein and name CA and (resid 12:328)', ref_frame=0).run()
ref = average.universe

aligner = align.AlignTraj(u, ref, select='protein and name CA and (resid 12:328)', in_memory=True).run()

c_alphas = protein_proa  # Use the selection for PROA
R = rms.RMSF(c_alphas).run()

plt.plot(c_alphas.resids, R.rmsf)
plt.xlabel('Residue number', fontsize=14, weight='bold')
plt.ylabel('RMSF ($\AA$)', fontsize=14, weight='bold')
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
# Set the x-axis limits and ticks
plt.xlim(12, 328)
plt.xticks([12, 20, 67, 74, 93, 110, 145, 152, 177, 204, 244, 257, 293, 302, 328], [12, 20, 67, 74, 93, 110, 145, 152, 177, 204, 244, 257, 293, 302, 328], rotation=90)

# Calculate the midpoints of each shaded region
midpoints = [(20 + 67) / 2, (74 + 93) / 2, (110 + 145) / 2, (152 + 177) / 2, (204 + 244) / 2, (257 + 293) / 2, (302 + 327) / 2]

# Add labels at a fixed y-coordinate on top of each box
for midpoint, label in zip(midpoints, ['TM1', 'TM2', 'TM3', 'TM4', 'TM5', 'TM6', 'TM7']):
    plt.text(midpoint, 4.5, label, ha='center', va='center', color='black', fontsize=14, weight='bold')

# Add shaded regions for each TM segment
plt.axvspan(20, 67, zorder=0, alpha=0.2, color='grey', label='TM1')
plt.axvspan(74, 93, zorder=0, alpha=0.2, color='grey', label='TM2')
plt.axvspan(110, 145, zorder=0, alpha=0.2, color='grey', label='TM3')
plt.axvspan(152, 177, zorder=0, alpha=0.2, color='grey', label='TM4')
plt.axvspan(204, 244, zorder=0, alpha=0.2, color='grey', label='TM5')
plt.axvspan(257, 293, zorder=0, alpha=0.2, color='grey', label='TM6')
plt.axvspan(302, 327, zorder=0, alpha=0.2, color='grey', label='TM7')

plt.show()
