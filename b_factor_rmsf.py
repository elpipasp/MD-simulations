import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

rmsf_file_path = '/users/nkb19202/CORRECT/FREE/REP_1/namd/rmsf_inactive1_50ns.txt' #md rmsf of inactive
rmsf_data = np.loadtxt(rmsf_file_path, usecols=(0, 1))
pdb_file_path = '/users/nkb19202/CORRECT/FREE/REP_1/namd/7br3_clean.pdb' #b-factor pdb

#extract B-factor values from PDB for each residue
def extract_b_factors(pdb_file):
    b_factors = {}
    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith('ATOM'):
                residue_id = int(line[22:26].strip())
                b_factor = float(line[60:66].strip())
                if residue_id not in b_factors:
                    b_factors[residue_id] = []
                b_factors[residue_id].append(b_factor)
    return b_factors
b_factors_per_residue = extract_b_factors(pdb_file_path)

#align RMSF data with residue numbers from the PDB file
aligned_rmsf_data = np.zeros((len(b_factors_per_residue), 2))
for i, residue_id in enumerate(b_factors_per_residue.keys()):
    aligned_rmsf_data[i, 0] = residue_id
    aligned_rmsf_data[i, 1] = rmsf_data[rmsf_data[:, 0] == residue_id, 1]

#normalize both RMSF and B-factor values
normalized_rmsf = (aligned_rmsf_data[:, 1] - np.mean(aligned_rmsf_data[:, 1])) / np.std(aligned_rmsf_data[:, 1])
normalized_b_factors = np.array([np.mean(b_factors) for b_factors in b_factors_per_residue.values()])
normalized_b_factors = (normalized_b_factors - np.mean(normalized_b_factors)) / np.std(normalized_b_factors)

#plot
plt.figure(figsize=(12, 8))
plt.plot(aligned_rmsf_data[:, 0], normalized_b_factors, label='B-Factors Inactive Crystal')
plt.plot(aligned_rmsf_data[:, 0], normalized_rmsf, label='RMSF Apo-GnRH1R')
plt.xlabel('Residue', fontsize=22, fontweight='bold', labelpad=10)
plt.ylabel('Normalised Value', fontsize=22, fontweight='bold')
plt.xlim(18, 328)

#quantitative comparison (Pearson correlation coefficient)
correlation_coefficient, p_value = pearsonr(normalized_b_factors, normalized_rmsf)
correlation_coefficient_rounded = round(correlation_coefficient, 2)
p_value_scientific = f"{p_value:.2e}"

legend_text = f'Pearson: {correlation_coefficient_rounded}\np-value: {p_value_scientific}'
plt.legend(title=legend_text, loc='upper right', fontsize=16, title_fontsize=16)
plt.xticks([18, 33, 66, 74, 104, 110, 145, 153, 178, 204, 244, 257, 293, 301, 328], rotation=90, fontsize=20)
midpoints = [(33 + 66) / 2, (74 + 104) / 2, (110 + 145) / 2, (153 + 178) / 2,
             (204 + 244) / 2,(257 + 293) / 2,(301 +328) /2]

for midpoint,label in zip(midpoints,['TM1','TM2','TM3','TM4','TM5','TM6','TM7']):
    plt.text(midpoint,max(normalized_b_factors),label,
             ha='center',va='center',color='black',fontsize=18)
plt.yticks(fontsize=18)
plt.axvspan(33 ,67,zorder=0,alpha=0.2,color='grey',label='TM1')
plt.axvspan(74 ,93,zorder=0,alpha=0.2,color='grey',label='TM2')
plt.axvspan(110 ,145,zorder=0,alpha=0.2,color='grey',label='TM3')
plt.axvspan(152 ,177,zorder=0,alpha=0.2,color='grey',label='TM4')
plt.axvspan(204 ,244,zorder=0,alpha=0.2,color='grey',label='TM5')
plt.axvspan(257 ,293,zorder=0,alpha=0.2,color='grey',label='TM6')
plt.axvspan(302 ,327,zorder=0,alpha=0.2,color='grey',label='TM7')
plt.tight_layout()
plt.show()
