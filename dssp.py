import MDAnalysis as mda
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#;oad trajectory and reference structure
u = mda.Universe("NTER_FREE_SC.psf", "D26.dcd")
protein_PROA = u.select_atoms("segid PROA")
#convert MDAnalysis Universe to mdtraj trajectory
traj = md.load("D26.dcd", top="NTER_FREE_SC.psf", stride=500)
protein_indices = protein_PROA.resids - 1
dssp = md.compute_dssp(traj, simplified=False)
protein_dssp = dssp[:, protein_indices]
structure_types, structure_counts = np.unique(protein_dssp, return_counts=True)
total_residues = len(protein_indices) * protein_dssp.shape[0]
structure_percentages = {s: np.sum(protein_dssp == s) / total_residues * 100 for s in structure_types}

colors = sns.color_palette("Set3", n_colors=len(structure_types))
plt.rcParams.update({'axes.labelsize': 16, 'xtick.labelsize': 14, 'ytick.labelsize': 14, 'legend.fontsize': 10})
plt.figure(figsize=(10, 6))
bars = plt.bar(structure_percentages.keys(), structure_percentages.values(), color=colors)
plt.xlabel('DSSP term', weight='bold')
plt.ylabel('Percentage (%)', weight='bold')
plt.title('Apo-GnRH1R', weight='bold', fontsize=16)

legend_labels = {'H': 'Alpha helix',
                 'B': 'Beta bridge',
                 'E': 'Extended strand',
                 'G': '3-helix (3/10 helix)',
                 'I': '5 helix (pi helix)',
                 'T': 'Hydrogen bonded turn',
                 'S': 'Bend',
                 ' ': 'Loops and irregular elements'}

legend_handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=colors[i], markersize=10) for i in range(len(structure_types))]
legend_texts = [f'{key}: {legend_labels[key]}' for key in structure_types]

plt.legend(legend_handles, legend_texts, title='DSSP Assignments', loc='upper right')

for bar in bars:
    yval = bar.get_height()
    plt.text(bar.get_x() + bar.get_width() / 2, yval, round(yval, 2), ha='center', va='bottom', fontsize=12)

plt.show()
