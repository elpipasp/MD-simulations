import os
import matplotlib.pyplot as plt

input_file_path = 'PROA_MEMB_hbonds-details_sorted.dat'
residue_lines = {}
with open(input_file_path, 'r') as file:
    for line in file:
        fields = line.strip().split()
        residue = fields[0].split('-')[0]  # Extract residue name
        if residue not in residue_lines:
            residue_lines[residue] = []
        residue_lines[residue].append(line)

output_dir = 'output_files'
os.makedirs(output_dir, exist_ok=True)

for residue, lines in residue_lines.items():
    output_file_path = os.path.join(output_dir, f'{residue}_outfit.txt')
    with open(output_file_path, 'w') as output_file:
        output_file.writelines(lines)

print("Output files created in the 'output_files' directory.")
occupancy_data = {}

for residue in residue_lines.keys():
    input_file_path = os.path.join(output_dir, f'{residue}_outfit.txt')
    with open(input_file_path, 'r') as input_file:
        total_occupancy = sum(float(line.split()[-1][:-1]) for line in input_file)
        occupancy_data[residue] = total_occupancy

labels = occupancy_data.keys()
sizes = occupancy_data.values()

plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90)
plt.axis('equal') 

plt.savefig('occupancy_pie_chart.png')
plt.show()
