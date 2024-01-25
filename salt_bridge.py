import matplotlib.pyplot as plt
import numpy as np
import re

data_files = [
    '',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
    '',


]

fig, ax1 = plt.subplots(figsize=(10, 8))
ax2 = ax1.twinx()

all_x_values = []
all_y_values = []
durations = []

color = 'blue'
markersize = 25

labels = []

for i, data_file in enumerate(data_files):
    with open(data_file, 'r') as file:
        data = [line.split() for line in file.readlines()]
        frames, distances = zip(*[(int(frame), float(distance)) for frame, distance in data])

    #filter data where distances are less than or equal to 3.5
    filtered_data = [(frame * 2 / 100, i * markersize + markersize / 2) for frame, distance in zip(frames, distances) if distance <= 3.5]

    #calculate duration
    duration_frames = len(filtered_data)

    #convert duration to nanoseconds
    duration_ns = duration_frames * 2 / 100

    all_x_values.extend([x for x, _ in filtered_data])
    all_y_values.extend([y for _, y in filtered_data])
    durations.append(duration_ns)

    #extract amino acid codes, chain information, and XXXnumber from datafile
    match = re.match(r'^saltbr-([A-Z]+)(\d+)_chain([A-Z]+)_segname[A-Z]+-([A-Z]+)(\d+)_chain([A-Z]+)_segname[A-Z]+.dat$', data_file)
    if match:
        amino_acid1 = match.group(1)
        number1 = match.group(2)
        chain1 = match.group(3)
        amino_acid2 = match.group(4)
        number2 = match.group(5)
        chain2 = match.group(6)

        amino_acid1 = amino_acid1.replace('GLU', 'E').replace('ARG', 'R').replace('ASP', 'D').replace('LYS', 'K')
        amino_acid2 = amino_acid2.replace('GLU', 'E').replace('ARG', 'R').replace('ASP', 'D').replace('LYS', 'K')

        label = f"{amino_acid1[0]}{number1}-{amino_acid2[0]}{number2}"
        labels.append(label)

        ax1.plot([x for x, _ in filtered_data], [y for _, y in filtered_data], '|', markersize=markersize, color=color, label=label)
    else:
        print(f"Skipping invalid filename: {data_file}")

ax1.set_xlabel('Time (ns)', fontsize=16, fontweight='bold')
ax1.set_ylabel('Salt bridge', color='black', fontsize=16, fontweight='bold')
ax1.set_title('Salt Bridges Over Time ROS-2 (1)', fontsize=16, fontweight='bold')
yticks = np.arange(markersize / 2, len(data_files) * markersize, markersize)
ax1.set_yticks(yticks)
ax1.set_yticklabels(labels)
ax1.set_ylim(min(all_y_values) - markersize / 2, max(all_y_values) + markersize / 2)
for ytick in yticks:
    ax1.axhline(y=ytick + markersize / 2, color='black', linestyle='-', linewidth=0.5, alpha=0.5)


x_max = max(all_x_values)
ax1.set_xlim(0, x_max)
yticks_right = yticks
yticklabels_right = [f'{duration:.2f}' for duration in durations]
ax2.set_yticks(yticks_right)
ax2.set_yticklabels(yticklabels_right)
ax2.set_ylabel('Duration (ns)', color='black', fontsize=16, fontweight='bold')
ax2.set_ylim(0, len(data_files) * markersize)
ax1.set_xticks(list(ax1.get_xticks()) + [x_max])
plt.tight_layout()
plt.show()
