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


pdb_filename = '7mt8.pdb' #rhodopsin
ca_coordinates = read_pdb_file(pdb_filename)

def calculate_tm_distances(tm1, tm2):
    indices_tm1 = tm_indices[tm1]
    indices_tm2 = tm_indices[tm2]
    ca_coords_tm1 = [np.array(ca_coordinates[(i, 'CA')]) for i in range(indices_tm1[0], indices_tm1[1] + 1)]
    ca_coords_tm2 = [np.array(ca_coordinates[(i, 'CA')]) for i in range(indices_tm2[0], indices_tm2[1] + 1)]
    distance_beginning = np.sqrt(np.sum((ca_coords_tm1[0] - ca_coords_tm2[0]) ** 2))
    distance_middle = np.sqrt(np.sum(((np.mean(ca_coords_tm1, axis=0) - np.mean(ca_coords_tm2, axis=0))) ** 2))
    distance_end = np.sqrt(np.sum((ca_coords_tm1[-1] - ca_coords_tm2[-1]) ** 2))

    return distance_beginning, distance_middle, distance_end

#calculate pairwise TM distances for all combinations
tm_list = list(tm_indices.keys())
for i in range(len(tm_list)):
    for j in range(i + 1, len(tm_list)):
        tm1 = tm_list[i]
        tm2 = tm_list[j]
        distances = calculate_tm_distances(tm1, tm2)
        print(f'Distances between {tm1} and {tm2}: Beginning: {distances[0]}, Middle: {distances[1]}, End: {distances[2]}')

# plot
fig, ax = plt.subplots()
for tm, indices in tm_indices.items():
    ca_coords = [np.array(ca_coordinates[(i, 'CA')]) for i in range(indices[0], indices[1] + 1)]
    ca_coords = np.array(ca_coords)
    ax.plot(ca_coords[:, 0], ca_coords[:, 1], label=tm)
ax.set_xlabel('X Coordinate')
ax.set_ylabel('Y Coordinate')
ax.legend()
plt.show()
