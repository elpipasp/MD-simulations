import numpy as np
import matplotlib.pyplot as plt

# define the residue indices for the first and last TMs
tm_indices = {
    'TM1': (33, 65),
    'TM2': (72, 101),
    'TM3': (107, 141),
    'TM4': (149, 174),
    'TM5': (199, 236),
    'TM6': (242, 278),
    'TM7': (284, 310)
}

# read PDB file and extract Cα atom coordinates based on residue indices
def read_pdb_file(pdb_filename):
    ca_coordinates = {}
    with open(pdb_filename, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith('ATOM') and line[13:15] == 'CA':
                residue_index = int(line[23:26].strip())
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                if any(start <= residue_index <= end for start, end in tm_indices.values()):
                    ca_coordinates[(residue_index, 'CA')] = (x, y, z)
    return ca_coordinates


pdb_filename = '7mt8.pdb'


ca_coordinates = read_pdb_file(pdb_filename)


def calculate_tm_angles(tm1, tm2):
    indices_tm1 = tm_indices[tm1]
    indices_tm2 = tm_indices[tm2]

    ca_coords_tm1 = [np.array(ca_coordinates[(i, 'CA')]) for i in range(indices_tm1[0], indices_tm1[1] + 1)]
    ca_coords_tm2 = [np.array(ca_coordinates[(i, 'CA')]) for i in range(indices_tm2[0], indices_tm2[1] + 1)]

    vector_tm1 = np.mean(ca_coords_tm1, axis=0)
    vector_tm2 = np.mean(ca_coords_tm2, axis=0)

    dot_product = np.dot(vector_tm1, vector_tm2)
    magnitude_tm1 = np.linalg.norm(vector_tm1)
    magnitude_tm2 = np.linalg.norm(vector_tm2)

    cosine_angle = dot_product / (magnitude_tm1 * magnitude_tm2)
    angle_degrees = np.arccos(cosine_angle) * 180 / np.pi

    return angle_degrees


tm_list = list(tm_indices.keys())

angles_matrix = np.zeros((len(tm_list), len(tm_list)))


for i, tm1 in enumerate(tm_list):
    for j, tm2 in enumerate(tm_list):
        if i != j:  
            angle = calculate_tm_angles(tm1, tm2)
            angles_matrix[i, j] = angle


for i, tm1 in enumerate(tm_list):
    plt.figure(figsize=(8, 6))
    angles = np.zeros(len(tm_list) - 1)
    labels = []
    x = []

    for j, tm2 in enumerate(tm_list):
        if i != j: 
            angles[j - (j > i)] = angles_matrix[i, j]
            labels.append(tm2)
            x.append(j - (j > i))

    width = 0.2 
    x = np.array(x)  

    plt.bar(x, angles, width=width, color='royalblue')

    plt.ylabel('Angle (degrees)', fontsize=14, fontweight='bold')
    plt.title(f'Angles of {tm1}', fontsize=14, fontweight='bold')

    for j in range(len(labels)):
        angle = angles[j]
        plt.annotate(f'{angle:.2f}', xy=(x[j], angle), ha='center', va='bottom', rotation=45, fontweight='bold')


    plt.xticks(x, labels, fontsize=12, fontweight='bold')
    plt.ylim(0, 2, 4, 6)
    plt.yticks([0, 7], fontsize=12)
    plt.tight_layout()
    plt.show()
