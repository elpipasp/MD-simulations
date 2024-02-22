import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

trajectory = md.load('t.dcd', top='t.pdb')

tm_indices = [
    (start_index_tm1, end_index_tm1),
    (start_index_tm2, end_index_tm2),
    (start_index_tm3, end_index_tm3),
    (start_index_tm4, end_index_tm4),
    (start_index_tm5, end_index_tm5),
    (start_index_tm6, end_index_tm6),
    (start_index_tm7, end_index_tm7),
]
#calculate the angle between two vectors
def calculate_angle(v1, v2):
    dot_product = np.dot(v1, v2)
    norm_v1 = np.linalg.norm(v1)
    norm_v2 = np.linalg.norm(v2)
    cos_theta = dot_product / (norm_v1 * norm_v2)
    angle_rad = np.arccos(np.clip(cos_theta, -1.0, 1.0))
    angle_deg = np.degrees(angle_rad)
    return angle_deg
#output file for writing
with open('vectors.dat', 'w') as output_file:
    output_file.write("TM_Helix_Angle_Degrees\n")
    # calculate angles for each TM helix
    for tm_start, tm_end in tm_indices:
        # extract coordinates of CA atoms for the TM helix
        tm_coords_start = trajectory.xyz[:, tm_start, :]
        tm_coords_end = trajectory.xyz[:, tm_end, :]
        #calculate vectors from CA atoms
        tm_vector = tm_coords_end - tm_coords_start
        z_axis_vector = np.array([0, 0, 1])  # Principal z-axis
        #calculate the angle between the TM vector and the z-axis
        angle_deg = calculate_angle(tm_vector, z_axis_vector)
        #write angle to the output file
        output_file.write(f"{angle_deg:.2f}\n")
#plot
angles = np.loadtxt('vectors.dat', skiprows=1)  # Skip header row
plt.plot(angles, marker='o', linestyle='-', label='TM Helix Angles')
plt.xlabel('TM Helix Index')
plt.ylabel('Angle (Degrees)')
plt.title('Angles of TM Helices with Principal Z-Axis')
plt.legend()
plt.show()
