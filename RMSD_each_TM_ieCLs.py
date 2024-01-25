import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

file_names = [
    'TM1_rmsd.dat', 'TM2_rmsd.dat', 'TM3_rmsd.dat',
    'TM4_rmsd.dat', 'TM5_rmsd.dat', 'TM6_rmsd.dat',
    'TM7_rmsd.dat'
]

plt.figure(figsize=(10, 6))
colors = ['red', 'blue', 'green', 'purple', 'orange', 'cyan', 'brown']
max_frame = 0

for i, file_name in enumerate(file_names):
    df = pd.read_csv(file_name, delim_whitespace=True)
    df['frame'] /= 50
    plt.plot(df['frame'], df['mol0'], label=file_name.split('_')[0], linewidth=1, color=colors[i])
    max_frame = max(max_frame, max(df['frame']))  # Update the maximum frame value

plt.xlim([0, 600])  
plt.ylim([0, 7])    
plt.xlabel('Time (ns)', fontsize=14, fontweight='bold')
plt.ylabel('RMSD (Å)', fontsize=14, fontweight='bold')
plt.title('7TM RMSD over Time (ROS-2 (1))', fontsize=16, fontweight='bold')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.gca().xaxis.set_major_locator(MultipleLocator(50))
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.grid(True)
plt.show()
