import matplotlib.pyplot as plt
import numpy as np

# List of data file names
data_files = [
    'saltbr-GLU68-LYS67.dat',
    'saltbr-GLU68-LYS69.dat',
    'saltbr-ASP293-ARG299.dat',
    'saltbr-GLU295-ARG299.dat',
    'saltbr-GLU68-LYS66.dat',
    'saltbr-GLU111-LYS115.dat',
    'saltbr-ASP185-LYS191.dat',
    'saltbr-GLU111-ARG179.dat',
    'saltbr-ASP98-LYS121.dat',
    'saltbr-ASP98-ARG38.dat',
    'saltbr-ASP138-ARG139.dat',
    'saltbr-GLU68-LYS71.dat',
    'saltbr-ASP138-ARG75.dat',
    'saltbr-GLU90-LYS121.dat',
    'saltbr-GLU248-ARG240.dat',
    'saltbr-GLU68-LYS72.dat',
]

# Create a new figure with twin y-axis
fig, ax1 = plt.subplots(figsize=(10, 8))
ax2 = ax1.twinx()

# Concatenate and filter all data
all_x_values = []
all_y_values = []
durations = []

# Assign different colors to each data file
colors = plt.cm.tab20.colors

# Define markersize
markersize = 30

for i, data_file in enumerate(data_files):
    # Read data from the file
    with open(data_file, 'r') as file:
        data = [line.split() for line in file.readlines()]
        frames, distances = zip(*[(int(frame), float(distance)) for frame, distance in data])

    # Filter data where distances are less than or equal to 3.2
    filtered_data = [(frame * 2 / 100, i * markersize + markersize / 2) for frame, distance in zip(frames, distances) if distance <= 3.2]

    # Calculate duration in frames
    duration_frames = len(filtered_data)

    # Convert duration to nanoseconds
    duration_ns = duration_frames * 2 / 100

    # Append to the overall data and durations
    all_x_values.extend([x for x, _ in filtered_data])
    all_y_values.extend([y for _, y in filtered_data])
    durations.append(duration_ns)

    # Plot salt bridge data with a straightforward stacking and assign color
    ax1.plot([x for x, _ in filtered_data], [y for _, y in filtered_data], '|', markersize=markersize, color=colors[i])

# Customize plot
ax1.set_xlabel('Time (ns)', fontsize=16, fontweight='bold')
ax1.set_ylabel('Salt bridge', color='black', fontsize=16, fontweight='bold')
ax1.set_title('Salt Bridges Over Time', fontsize=16, fontweight='bold')

# Set y-axis ticks and labels for the left y-axis at the midpoint of the gridlines
yticks = np.arange(markersize/2, len(data_files) * markersize, markersize)
yticklabels = [file.split('saltbr-')[1].replace('.dat', '').replace('-', ' - ') for file in data_files]
ax1.set_yticks(yticks)
ax1.set_yticklabels(yticklabels)

# Set y-axis limits to cover the entire range of data points
ax1.set_ylim(min(all_y_values) - markersize / 2, max(all_y_values) + markersize / 2)

# Add grid lines at the midpoints of consecutive y-ticks
for ytick in yticks:
    ax1.axhline(y=ytick + markersize / 2, color='black', linestyle='-', linewidth=0.5, alpha=0.5)

# Set x-axis limits from 0 to the maximum value of the 'x_values' column
ax1.set_xlim(0, max(all_x_values))

# Set y-axis ticks and labels for the right y-axis at the midpoint of the markersize
yticks_right = yticks
yticklabels_right = [f'{duration:.2f}' for duration in durations]
ax2.set_yticks(yticks_right)
ax2.set_yticklabels(yticklabels_right)

# Plot durations on the right y-axis
ax2.set_ylabel('Duration (ns)', color='Black', fontsize=16, fontweight='bold')
ax2.set_ylim(0, len(data_files) * markersize)

plt.tight_layout()  # Adjust layout to prevent clipping of labels
plt.show()
