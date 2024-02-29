import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

file_names = [
    'rmsd_1to18_D1_D22.dat', 'rmsd_18to33_D1_D22.dat', 'TM1_rmsd.dat',
    'TM2_rmsd.dat', 'TM3_rmsd.dat', 'TM4_rmsd.dat',
    'TM5_rmsd.dat', 'TM6_rmsd.dat', 'TM7_rmsd.dat'
]

plt.figure(figsize=(10, 6))
colors = ['red', 'blue']
custom_colors = ['green', 'purple', 'orange', 'cyan', 'brown', 'magenta', 'teal']

max_frame = 0
for i, file_name in enumerate(file_names):
    df = pd.read_csv(file_name, delim_whitespace=True)
    df['frame'] /= 50
    if i == 0: 
        label = '1-17 N-ter'
    elif i == 1: 
        label = '18-33 N-ter'
    else:
        label = file_name.split('_')[0]
    if i < len(colors):
        plt.plot(df['frame'], df['mol0'], label=label, linewidth=1, color=colors[i])
    else:
        plt.plot(df['frame'], df['mol0'], label=label, linewidth=1, color=custom_colors[i-len(colors)])
    max_frame = max(max_frame, max(df['frame'])) 

plt.xlim([0, 1100])  
plt.ylim([0.1, 15])    
plt.xlabel('Time (ns)', fontsize=14, fontweight='bold')
plt.ylabel('RMSD (Å)', fontsize=14, fontweight='bold')
plt.title('ROS-1 #2', fontsize=16, fontweight='bold')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.gca().xaxis.set_major_locator(MultipleLocator(100))
plt.legend(loc='upper right')  
plt.grid(False)  
plt.show()
