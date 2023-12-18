import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

# Read the data from the 'rep_1_RMSD_450ns.dat' file into a DataFrame
df = pd.read_csv('rep_1_RMSD_450ns.dat', delim_whitespace=True)

# Divide values in column 'frame' by 50
df['frame'] /= 50

# Plotting the 'frame' against 'mol0' with a simple line
plt.plot(df['frame'], df['mol0'], linestyle='-', linewidth=2)

# Set custom limits for x and y axes
plt.xlim([-1, max(df['frame'])])  # Customize the x-axis limits
plt.ylim([0, 5])    # Customize the y-axis limits

# Set font size and make titles bold
plt.xlabel('Time (ns)', fontsize=14, fontweight='bold')
plt.ylabel('RMSD (Å)', fontsize=14, fontweight='bold')
plt.title('C-alpha RMSD over Time', fontsize=16, fontweight='bold')

# Increase tick label font size
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# Set x-axis ticks every 50 units
plt.gca().xaxis.set_major_locator(MultipleLocator(50))

plt.grid(True)
plt.show()
