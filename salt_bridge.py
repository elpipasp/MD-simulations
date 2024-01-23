import matplotlib.pyplot as plt
import numpy as np

data_files = [
    'saltbr-GLU111-LYS115.dat',
    'saltbr-GLU111-ARG179.dat',
    'saltbr-ASP185-LYS191.dat',
    'saltbr-ASP293-ARG299.dat',
    'saltbr-GLU295-ARG299.dat',
    'saltbr-GLU68-LYS77.dat',
    'saltbr-GLU68-LYS67.dat',
    'saltbr-GLU68-LYS69.dat',
    'saltbr-ASP98-LYS121.dat',
    'saltbr-GLU68-LYS71.dat',
    'saltbr-GLU68-LYS66.dat',
    'saltbr-ASP138-ARG139.dat',
    'saltbr-GLU248-ARG240.dat',
    'saltbr-GLU68-LYS62.dat',
    'saltbr-GLU68-LYS150.dat',
]


fig, ax1 = plt.subplots(figsize=(10, 8))
ax2 = ax1.twinx()


all_x_values = []
all_y_values = []
durations = []

color = 'blue'

# Define markersize
markersize = 30


labels = []

for i, data_file in enumerate(data_files):
    # Read data from the file
    with open(data_file, 'r') as file:
        data = [line.split() for line in file.readlines()]
        frames, distances = zip(*[(int(frame), float(distance)) for frame, distance in data])

    # filter data where distances are less than or equal to 3.5
    filtered_data = [(frame * 2 / 100, i * markersize + markersize / 2) for frame, distance in zip(frames, distances) if distance <= 3.5]

#calculate duration
    duration_frames = len(filtered_data)

    # convert duration to nanoseconds
    duration_ns = duration_frames * 2 / 100


    all_x_values.extend([x for x, _ in filtered_data])
    all_y_values.extend([y for _, y in filtered_data])
    durations.append(duration_ns)

    # extract amino acid codes from the file name and replace with 1-letter codes
    file_parts = data_file.split('-')
    amino_acid1 = file_parts[1].replace('GLU', 'E').replace('ARG', 'R').replace('ASP', 'D').replace('LYS', 'K')
    amino_acid2 = file_parts[2].split('.')[0].replace('GLU', 'E').replace('ARG', 'R').replace('ASP', 'D').replace('LYS', 'K')

    # create labels with 1-letter amino acid codes 
    label = amino_acid1[0] + file_parts[1][3:] + '-' + amino_acid2[0] + file_parts[2].split('.')[0][3:]
    labels.append(label)


    ax1.plot([x for x, _ in filtered_data], [y for _, y in filtered_data], '|', markersize=markersize, color=color, label=label)


ax1.set_xlabel('Time (ns)', fontsize=16, fontweight='bold')
ax1.set_ylabel('Salt bridge', color='black', fontsize=16, fontweight='bold')
ax1.set_title('Salt Bridges Over Time ROS-1 (1)', fontsize=16, fontweight='bold')


yticks = np.arange(markersize/2, len(data_files) * markersize, markersize)
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
ax2.set_ylabel('Duration (ns)', color='Black', fontsize=16, fontweight='bold')
ax2.set_ylim(0, len(data_files) * markersize)

ax1.set_xticks(list(ax1.get_xticks()) + [x_max])

plt.tight_layout()  
plt.show()
