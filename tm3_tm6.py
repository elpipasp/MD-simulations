import pandas as pd
import matplotlib.pyplot as plt

csv_file_path = '/Users/nikipaspali/Desktop/MD_5/distances2.csv' 
df = pd.read_csv(csv_file_path, skiprows=2)

plt.figure(figsize=(6, 8))

plt.subplot(2, 1, 1)
plt.plot(df.iloc[:, 0], df.iloc[:, 3], label='ROS-1 #1', color='darkblue')
plt.plot(df.iloc[:, 0], df.iloc[:, 2], label='ROS-1 #2', color='blue')
plt.plot(df.iloc[:, 0], df.iloc[:, 1], label='ROS-1 #3', color='deepskyblue')
plt.axhline(y=7.5, color='black', linestyle='--')  # Line at y=8 is inactive crystal
plt.axhline(y=11, color='green', linestyle='--')  # Line at y=11 is active AF 
plt.xlim(850, 1100)
plt.xticks(range(850, 1101, 50))
plt.yticks(range(6, 15)) 
plt.tick_params(axis='both', which='major', labelsize=12)
plt.axvline(x=1000, color='black', linestyle='-')
plt.axvline(x=850, color='black', linestyle='-')

plt.xlabel('Time (ns)', fontsize=14, fontweight='bold')
plt.ylabel('$\mathbf{R139^{3.50}-T265^{6.33}}$', fontsize=14, fontweight='bold')
plt.legend(loc='upper left')

plt.subplot(2, 1, 2)
plt.plot(df.iloc[:, 0], df.iloc[:, 8], label='Free-GnRH1R #1', color='k')
plt.plot(df.iloc[:, 0], df.iloc[:, 9], label='Free-GnRH1R #2', color='grey')
plt.axhline(y=7.5, color='black', linestyle='--')  # Line at y=8 is inactive crystal
plt.axhline(y=11, color='green', linestyle='--')  # Line at y=11 is active AF 
plt.xlim(850, 1100)
plt.xticks(range(850, 1101, 50))
plt.yticks(range(6, 15)) 
plt.tick_params(axis='both', which='major', labelsize=12)
plt.xlabel('Time (ns)', fontsize=14, fontweight='bold')
plt.ylabel('$\mathbf{R139^{3.50}-T265^{6.33}}$', fontsize=14, fontweight='bold')
plt.legend(loc='upper left')

plt.tight_layout()

plt.savefig('/Users/nikipaspali/Desktop/MD_5/tms.pdf')
plt.show()
