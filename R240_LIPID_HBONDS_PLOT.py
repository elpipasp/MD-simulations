import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt('R240_hbonds.dat')#from vmd D19-D22.dcd
time = data[:, 0] / 50.0 + 900.0  #to ns and add 900 because D19 is 900ns
num_hydrogen_bonds = data[:, 1]

plt.figure(figsize=(10, 6))
plt.plot(time, num_hydrogen_bonds, linestyle='-', markersize=3, color='grey')
plt.xlabel('Time (ns)', fontsize=16, fontweight='bold')
plt.ylabel('Hydrogen Bonds', fontsize=16, fontweight='bold')
plt.title('R240-Lipid Hydrogen Bonds (ROS-1 #2)', fontsize=16)
plt.grid(True)
plt.xlim(900, 1100)
plt.ylim(0, 2)
plt.yticks([1, 2], fontsize=14)
plt.xticks(fontsize=14)
plt.show()
