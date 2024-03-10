import numpy as np
import matplotlib.pyplot as plt

data = []
with open("pipistack.txt", "r") as file: #grep PIPISTACK etc first to create clean files
    for line in file:
        parts = line.strip().split("\t")
        residue1 = int(parts[0].split(":")[1])
        residue2 = int(parts[2].split(":")[1])
        distance = float(parts[3])
        data.append((residue1, residue2, distance))

residues = sorted(set(residue for pair in data for residue in pair[:2]))
contact_map = np.zeros((len(residues), len(residues)))

for pair in data:
    i = residues.index(pair[0])
    j = residues.index(pair[1])
    contact_map[i, j] = pair[2]
    #mirror the contact map along the diagonal
    contact_map[j, i] = pair[2]

#plot
plt.figure(figsize=(10, 8))
img = plt.imshow(contact_map, cmap='viridis', interpolation='nearest', vmin=0, vmax=np.max(contact_map), extent=[0, len(residues), 0, len(residues)], origin='lower', aspect='auto')
plt.colorbar()
tick_positions = np.arange(len(residues)) + 0.5
plt.xticks(ticks=tick_positions, labels=residues, rotation='vertical', fontsize=12)
plt.yticks(ticks=tick_positions, labels=residues, fontsize=12)

magenta_ticks = [2, 3, 5]
for tick in magenta_ticks:
    plt.xticks()[1][residues.index(tick)].set_color('magenta')
    plt.yticks()[1][residues.index(tick)].set_color('magenta')

plt.axhline(y=np.where(np.array(residues) == 5)[0][0] + 1, color='magenta', linestyle='--', linewidth=2)
plt.axvline(x=np.where(np.array(residues) == 5)[0][0] + 1, color='magenta', linestyle='--', linewidth=2)

#grid
for i in range(len(residues) + 1):
    plt.axhline(y=i, color='black', linestyle='-', linewidth=0.5)
    plt.axvline(x=i, color='black', linestyle='-', linewidth=0.5)

plt.title("ROS-1 π-π stacking", fontsize=18, fontweight='bold')
plt.xlabel("GnRH/GnRH1R residues", fontsize=16, fontweight='bold')
plt.ylabel("GnRH/GnRH1R residues", fontsize=16, fontweight='bold')
plt.gca().set_facecolor('white')
plt.show()
