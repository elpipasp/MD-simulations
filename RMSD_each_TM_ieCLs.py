import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

# List of file names
file_names = [
    'rep_1_RMSD_450ns_ICL1.dat', 'rep_1_RMSD_450ns_TM1.dat', 'rep_1_RMSD_450ns_NTER.dat',
    'rep_1_RMSD_450ns_TM2.dat', 'rep_1_RMSD_450ns_ECL1.dat', 'rep_1_RMSD_450ns_TM3.dat',
    'rep_1_RMSD_450ns_ICL2.dat', 'rep_1_RMSD_450ns_TM4.dat', 'rep_1_RMSD_450ns_ECL2.dat',
    'rep_1_RMSD_450ns_TM5.dat', 'rep_1_RMSD_450ns_ICL3.dat', 'rep_1_RMSD_450ns_TM6.dat',
    'rep_1_RMSD_450ns_ECL3.dat', 'rep_1_RMSD_450ns_TM7.dat'
]

# Create a new figure
plt.figure(figsize=(10, 6))

# List of distinct colors
colors = ['red', 'blue', 'green', 'purple', 'orange', 'cyan', 'brown', 'pink', 'gray', 'olive', 'lime', 'teal', 'magenta', 'yellow']

# Loop through each file and plot with a different color
for i, file_name in enumerate(file_names):
    df = pd.read_csv(file_name, delim_whitespace=True)
    df['frame'] /= 50
    plt.plot(df['frame'], df['mol0'], label=file_name.split('_')[-1].split('.')[0], linewidth=2, color=colors[i])

# Set custom limits for x and y axes
plt.xlim([-1, max(df['frame'])])  # Customize the x-axis limits
plt.ylim([0, 7])    # Customize the y-axis limits

# Set font size and make titles bold
plt.xlabel('Time (ns)', fontsize=14, fontweight='bold')
plt.ylabel('RMSD (Å)', fontsize=14, fontweight='bold')
plt.title('C-alpha RMSD over Time', fontsize=16, fontweight='bold')

# Increase tick label font size
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# Set x-axis ticks every 50 units
plt.gca().xaxis.set_major_locator(MultipleLocator(50))

# Add legend outside the right side of the plot
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

plt.grid(True)
plt.show()
