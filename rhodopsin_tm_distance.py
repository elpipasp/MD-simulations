import numpy as np
import matplotlib.pyplot as plt

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


pdb_filename = '7mt8.pdb' #active rhodopsin

s
ca_coordinates = read_pdb_file(pdb_filename)

# function to calculate pairwise TM distances in Ångströms
def calculate_tm_distances(tm1, tm2):
    indices_tm1 = tm_indices[tm1]
    indices_tm2 = tm_indices[tm2]

    ca_coords_tm1 = [np.array(ca_coordinates[(i, 'CA')]) for i in range(indices_tm1[0], indices_tm1[1] + 1)]
    ca_coords_tm2 = [np.array(ca_coordinates[(i, 'CA')]) for i in range(indices_tm2[0], indices_tm2[1] + 1)]

    # Calculate distances to Angstroms
    distance_beginning = np.linalg.norm(ca_coords_tm1[0] - ca_coords_tm2[0])
    distance_middle = np.linalg.norm(np.mean(ca_coords_tm1, axis=0) - np.mean(ca_coords_tm2, axis=0))
    distance_end = np.linalg.norm(ca_coords_tm1[-1] - ca_coords_tm2[-1])

    return distance_beginning, distance_middle, distance_end

# create a list of TMs
tm_list = list(tm_indices.keys())

# create a separate plot for each TM comparison, excluding self-comparisons
for i, tm1 in enumerate(tm_list):
    plt.figure(figsize=(8, 6))
    distances = np.zeros((len(tm_list) - 1, 3))  # 3 distances for each pair
    labels = []
    x = []

    for j, tm2 in enumerate(tm_list):
        if i != j:  # avoid calculating distances for the same TM
            distances[j - (j > i)] = calculate_tm_distances(tm1, tm2)
            labels.append(tm2)
            x.append(j - (j > i))

    width = 0.2  
    x = np.array(x)  

    colors = ['grey', 'rosybrown', 'tan']

    for k in range(3):  # Beginning, Middle, End
        plt.bar(x + k * width, distances[:, k], width=width, label=f'{["Beginning", "Middle", "End"][k]}', color=colors[k])


    plt.ylabel('Distance (Å)', fontsize=14, fontweight='bold')
    plt.title(f'Distances of {tm1}', fontsize=14, fontweight='bold')
    
    for j in range(len(labels)):
        for k in range(3):
            distance = distances[j, k]
            plt.annotate(f'{distance:.2f}', xy=(x[j] + k * width, distance), ha='center', va='bottom', rotation=90)

    plt.xticks(x + width, labels, fontsize=12, fontweight='bold')
    plt.ylim(0, 100)
    plt.yticks([0, 10, 20, 30, 40, 50, 60], fontsize=12)
    plt.legend()
    plt.tight_layout()
    plt.show()
