import matplotlib.pyplot as plt

data_file = "D20_D21_HBONDShbonds-details_sorted.dat"
labels = []
occupancies = []

def custom_autopct(pct):
    total = sum(occupancies)
    occupancy = round(pct * total / 100.0, 2)
    return f'{occupancy:.2f}%'

with open(data_file, 'r') as file:
    next(file)  #skip the header
    for line in file:
        line = line.strip().split('\t')
        donor = line[0]
        acceptor = line[1]
        occupancy = float(line[2].rstrip('%'))  #remove % and convert to float
        if occupancy > 3.0:
            label = f"{donor}-{acceptor}"
            labels.append(label)
            occupancies.append(occupancy)

fig, ax = plt.subplots(figsize=(8, 8))
colors = plt.cm.tab20c(range(len(occupancies)))
patches, texts, autotexts = ax.pie(occupancies, labels=labels, startangle=140, colors=colors, autopct=custom_autopct)

for text, color in zip(texts, colors):
    text.set_color(color) 
    text.set_fontsize(16)
    text.set_fontweight('bold') 
ax.set_title('Hydrogen Bond Occupancy (%)', pad=20, fontsize=18)
plt.axis('equal')  
plt.show()
