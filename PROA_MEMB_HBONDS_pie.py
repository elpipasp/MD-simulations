import os
import matplotlib.pyplot as plt
from collections import defaultdict

data_directory = "/Users/nikipaspali/Desktop/PIPER"
data_file_name = "PROA_MEMB_hbonds-details_sorted.dat"
data_file = os.path.join(data_directory, data_file_name)
residue_data = defaultdict(list)

def custom_autopct(pct, residue):
    total = sum(entry['occupancy'] for entry in residue_data[residue])
    occupancy = round(pct * total / 100.0, 2)
    return f'{occupancy:.2f}%'

with open(data_file, 'r') as file:
    next(file)  #skip the header
    for line in file:
        line = line.strip().split('\t')

        if len(line) >= 3:
            donor = line[0]
            acceptor = line[1]
            occupancy = float(line[2].rstrip('%'))  
            residue_number = donor.split('-')[0]
            
            if occupancy > 5.0:
                label = f"{donor}-{acceptor}"
                residue_data[residue_number].append({'label': label, 'occupancy': occupancy})

#plot
fig, axes = plt.subplots(4, 7, figsize=(28, 16))
fig.subplots_adjust(hspace=0.5)
flat_axes = axes.flatten()
sorted_residue_data = sorted(residue_data.items(), key=lambda x: sum(entry['occupancy'] for entry in x[1]))

#plots with less/more than 4 legend entries
less_than_four = []
four_or_more = []

for residue, data in sorted_residue_data:
    if len(data) < 4:
        less_than_four.append((residue, data))
    else:
        four_or_more.append((residue, data))
final_plot_order = less_than_four + four_or_more

for (residue, data), ax in zip(final_plot_order, flat_axes):
    labels = [entry['label'] for entry in data]
    occupancies = [entry['occupancy'] for entry in data]

    colors = plt.cm.tab20c(range(len(occupancies)))
    patches, texts, autotexts = ax.pie(occupancies, startangle=140, colors=colors, autopct=lambda pct: custom_autopct(pct, residue))
    
    for text, color in zip(texts, colors):
        text.set_color(color)  
        text.set_fontsize(10)
        text.set_fontweight('bold') 

    for autotext in autotexts:
        autotext.set_fontsize(10)

    ax.set_title(f'{residue}', pad=1, fontsize=16)
    ax.legend(labels, loc='upper center', bbox_to_anchor=(0.5, -0.02), fontsize=10)
    plt.axis('equal')  

#save
output_file_path = os.path.join(data_directory, f"{os.path.splitext(data_file_name)[0]}_multiplot.pdf")
plt.savefig(output_file_path, bbox_inches='tight')
print(f"Multiplot PDF saved successfully: {output_file_path}")
plt.show()
