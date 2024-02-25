import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

dcd_files = ['D18.dcd', 'D19.dcd', 'D20.dcd', 'D21.dcd', 'D22.dcd']

psf_file = 'ROS_1_SC.psf'
u = mda.Universe(psf_file)

#define selections for CA atoms of residues 139 (TM3), 265 (TM6), and 323 (TM7) in segname PROA
tm3_selection = u.select_atoms("segid PROA and resid 139 and name CA")
tm6_selection = u.select_atoms("segid PROA and resid 265 and name CA")
tm7_selection = u.select_atoms("segid PROA and resid 323 and name CA")

tm3_tm6_distances = []
tm3_tm7_distances = []
time_values = []

previous_time = 0  #initialise time from previous trajectory
for dcd_file in dcd_files:
    u.load_new(dcd_file)
  
    for ts in u.trajectory:
        #calculate distances between CA atoms of TM3 and TM6
        distance_tm3_tm6 = np.linalg.norm(tm3_selection.positions - tm6_selection.positions)
        tm3_tm6_distances.append(distance_tm3_tm6)
        time_values.append(ts.time + previous_time)
    previous_time = time_values[-1]
time_values = np.array(time_values)
tm3_tm6_distances = np.array(tm3_tm6_distances)

output_file = "tm3_tm6_distances.txt"
np.savetxt(output_file, np.column_stack((time_values, tm3_tm6_distances)), fmt='%10.5f', header='Time (ps)\tTM3-TM6 Distance (angstrom)', comments='')

# Plot
plt.plot(time_values, tm3_tm6_distances, color='b', linestyle='-')
plt.axhline(y=8, color='k', linestyle='--')
plt.xlabel('Time (ps)')
plt.ylabel('TM3-TM6 Distance (angstrom)')
plt.title('TM3-TM6 Distance over Time')
plt.grid(True)
plt.show()
