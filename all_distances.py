from Bio.PDB import PDBParser
import matplotlib.pyplot as plt
import numpy as np

def calculate_distance(pdb_file, chain_id, residue1, residue2):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    model = structure[0] 


    for chain in model:
        if chain.id == chain_id:
            ca_atom1 = None
            ca_atom2 = None

            for residue in chain.get_residues():
                if residue.id[1] == residue1:
                    for atom in residue:
                        if atom.name == 'CA':
                            ca_atom1 = atom
                            break
                elif residue.id[1] == residue2:
                    for atom in residue:
                        if atom.name == 'CA':
                            ca_atom2 = atom
                            break

                if ca_atom1 and ca_atom2:
                    break

            if ca_atom1 is None or ca_atom2 is None:
                return None

            # calculate the distance between the CA atoms
            distance = ca_atom1 - ca_atom2
            return distance

    return None

pdb_files = [
    {"file": "Inactive_GnRH1R.pdb", "chain_id": "P", "TM3": 139, "TM6": 262, "TM7": 323},
    {"file": "Inactive_940ns.pdb", "chain_id": "P", "TM3": 139, "TM6": 262, "TM7": 323},
    {"file": "1GNRH_1000ns.pdb", "chain_id": "P", "TM3": 139, "TM6": 262, "TM7": 323},
    {"file": "2GNRH_1000ns.pdb", "chain_id": "P", "TM3": 139, "TM6": 262, "TM7": 323},
    {"file": "1ROS1_600ns.pdb", "chain_id": "P", "TM3": 139, "TM6": 262, "TM7": 323},
    {"file": "2ROS1_600ns.pdb", "chain_id": "P", "TM3": 139, "TM6": 262, "TM7": 323},
    {"file": "1ROS2_600ns.pdb", "chain_id": "P", "TM3": 139, "TM6": 262, "TM7": 323},
    {"file": "2ROS2_600ns.pdb", "chain_id": "P", "TM3": 139, "TM6": 262, "TM7": 323},
    {"file": "1u19.pdb", "chain_id": "A", "TM3": 135, "TM6": 247, "TM7": 306},
    {"file": "7mt8.pdb", "chain_id": "R", "TM3": 135, "TM6": 247, "TM7": 306},
]

distances_TM3_TM6_list = []
distances_TM3_TM7_list = []
file_labels = []

for pdb_info in pdb_files:
    pdb_file = pdb_info["file"]
    chain_id = pdb_info["chain_id"]
    TM3 = pdb_info["TM3"]
    TM6 = pdb_info["TM6"]
    TM7 = pdb_info["TM7"]

    distance_TM3_TM6 = calculate_distance(pdb_file, chain_id, TM3, TM6)
    distance_TM3_TM7 = calculate_distance(pdb_file, chain_id, TM3, TM7)

    file_labels.append(pdb_file)

    if distance_TM3_TM6 is not None:
        distances_TM3_TM6_list.append(distance_TM3_TM6)
    else:
        distances_TM3_TM6_list.append(np.nan)

    if distance_TM3_TM7 is not None:
        distances_TM3_TM7_list.append(distance_TM3_TM7)
    else:
        distances_TM3_TM7_list.append(np.nan)


bar_positions = np.arange(len(file_labels))
bar_width = 0.35
fig, ax = plt.subplots(figsize=(12, 6))
bar1 = ax.bar(bar_positions - bar_width/2, distances_TM3_TM6_list, bar_width, label='TM3-TM6 Distances', color='orange')
bar2 = ax.bar(bar_positions + bar_width/2, distances_TM3_TM7_list, bar_width, label='TM3-TM7 Distances', color='blue')

ax.set_xticks(bar_positions)
ax.set_xticklabels(file_labels, rotation=45, fontsize=8, ha='right')

for bar1, bar2 in zip(bar1, bar2):
    height1 = bar1.get_height()
    height2 = bar2.get_height()
    ax.annotate(f'{height1:.2f} Å', xy=(bar1.get_x() + bar1.get_width() / 2, height1),
                xytext=(0, 3), textcoords='offset points', ha='center', va='bottom', fontsize=8)
    ax.annotate(f'{height2:.2f} Å', xy=(bar2.get_x() + bar2.get_width() / 2, height2),
                xytext=(0, 3), textcoords='offset points', ha='center', va='bottom', fontsize=8)

ax.set_ylabel('Distance (Å)')
ax.set_title('Distances Between CA Atoms for TM3-TM6 and TM3-TM7')
ax.legend()
plt.tight_layout()
plt.show()
