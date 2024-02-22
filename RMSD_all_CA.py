import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

df = pd.read_csv('rep_1_RMSD_450ns.dat', delim_whitespace=True)
df['frame'] /= 50
plt.plot(df['frame'], df['mol0'], linestyle='-', linewidth=2)
plt.xlim([-1, max(df['frame'])])
plt.ylim([0, 5])  
plt.xlabel('Time (ns)', fontsize=14, fontweight='bold')
plt.ylabel('RMSD (Å)', fontsize=14, fontweight='bold')
plt.title('C-alpha RMSD over Time', fontsize=16, fontweight='bold')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.gca().xaxis.set_major_locator(MultipleLocator(50))
plt.grid(True)
plt.show()
