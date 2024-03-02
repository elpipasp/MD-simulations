import mdtraj as md
import matplotlib.pyplot as plt

trajectory_files = ["D1.dcd", "D2.dcd", "D3.dcd", "D4.dcd", "D5.dcd", "D6.dcd", "D7.dcd", "D8.dcd", "D9.dcd", "D10.dcd", "D11.dcd", "D12.dcd", "D13.dcd", "D14.dcd", "D15.dcd", "D16.dcd", "D17.dcd", "D18.dcd", "D19.dcd", "D20.dcd", "D21.dcd", "D22.dcd",]
topology_file = "ROS_1_SC.pdb"
combined_rg = []
combined_time = []

#loop through each trajectory
for traj_file in trajectory_files:
    trajectory = md.load(traj_file, top=topology_file)
    gpcr_atoms = trajectory.top.select("protein")

    #ensure dcd continuity
    if combined_time:
        last_time = combined_time[-1]
        time_shift = last_time + (trajectory.time[0] - trajectory.time[0])
        trajectory.time += time_shift

    #calculate the radius of gyration for the GPCR
    rg = md.compute_rg(trajectory.atom_slice(gpcr_atoms))
    combined_rg.extend(rg)
    combined_time.extend(trajectory.time)
combined_time = [time / 50 for time in combined_time] #time to ns

#save output
output_file = "rg_2ROS1_D1_D22.txt"
with open(output_file, "w") as file:
    file.write("Time (ns)\tRadius of Gyration (nm)\n")
    for time, rg in zip(combined_time, combined_rg):
        file.write(f"{time}\t{rg}\n")

#plot
plt.plot(combined_time, combined_rg, color='deepskyblue', linewidth=1)
plt.xlabel("Time (ns)", fontsize=16, fontweight='bold')
plt.ylabel("Radius of Gyration (nm)", fontsize=16, fontweight='bold')
plt.title("ROS-1 #1", fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.show()
