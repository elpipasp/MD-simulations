import matplotlib.pyplot as plt

fig = plt.figure(figsize=(10, 10))

data_files = [
    {"file": "rmsf_inactive1_50ns.txt", "label": "Inactive GnRH1R", "color": "black"},
    {"file": "rmsf_1ROS2.txt", "label": "ROS-2 (1)", "color": "red"},
    {"file": "rmsf_2ROS2.txt", "label": "ROS-2 (2)", "color": "green"},
]


legend_handles = []
for data_info in data_files:
    with open(data_info["file"], "r") as file:
        lines = file.readlines()
    x_values = []
    y_values = []
    for line in lines:
        parts = line.split()
        if len(parts) == 2:
            x_values.append(int(parts[0]))
            y_values.append(float(parts[1]))
    plt.plot(x_values, y_values, label=data_info["label"], color=data_info["color"])
    legend_handles.append(plt.Line2D([0], [0], color=data_info["color"], label=data_info["label"]))
plt.xlabel('Residue number', fontsize=14, weight='bold')
plt.ylabel('RMSF ($\AA$)', fontsize=14, weight='bold')
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
plt.xlim(1, 328)
plt.xticks([33, 66, 74, 104, 110, 145, 153, 178, 204, 244, 257, 293, 301, 327], rotation=90)
plt.ylim(0, 10)
midpoints = [(33 + 66) / 2, (74 + 104) / 2, (110 + 145) / 2, (153 + 178) / 2, (204 + 244) / 2, (257 + 293) / 2, (301 + 327) / 2]
top_y_coordinate = 8.0 
for midpoint, label in zip(midpoints, ['TM1', 'TM2', 'TM3', 'TM4', 'TM5', 'TM6', 'TM7']):
    plt.text(midpoint, top_y_coordinate, label, ha='center', va='bottom', color='black', fontsize=14, weight='bold')

plt.axvspan(33, 66, zorder=0, alpha=0.2, color='grey', label='TM1')
plt.axvspan(74, 104, zorder=0, alpha=0.2, color='grey', label='TM2')
plt.axvspan(110, 145, zorder=0, alpha=0.2, color='grey', label='TM3')
plt.axvspan(153, 178, zorder=0, alpha=0.2, color='grey', label='TM4')
plt.axvspan(204, 244, zorder=0, alpha=0.2, color='grey', label='TM5')
plt.axvspan(257, 293, zorder=0, alpha=0.2, color='grey', label='TM6')
plt.axvspan(301, 327, zorder=0, alpha=0.2, color='grey', label='TM7')

plt.legend(handles=legend_handles, fontsize=12)
plt.show()
