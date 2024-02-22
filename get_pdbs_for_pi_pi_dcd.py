import mdtraj as md
import os

def save_frames_as_pdb(input_dcd, topology_file, stride=50, output_dir='pdb_files'):
    # Create output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
      
    topology = md.load(topology_file).topology
    traj = md.load_dcd(input_dcd, top=topology)

    #iterate over frames with the specified stride
    for i, frame in enumerate(traj[::stride]):
        #select atoms with segname PROA or PROB
        selected_atoms = frame.atom_slice(frame.top.select("segname PROA or segname PROB"))

        #replace segname PROA with chain id P and PROB with chain id A
        for atom in selected_atoms.topology.atoms:
            if atom.segment_id == 'PROA':
                atom.chain_id = 'P'
            elif atom.segment_id == 'PROB':
                atom.chain_id = 'A'

        #replace residue HSD with HIS
        for residue in selected_atoms.topology.residues:
            if residue.name == 'HSD':
                residue.name = 'HIS'

        #generate PDB filename
        pdb_filename = os.path.join(output_dir, f"frame_{i*stride}.pdb")

        #save selected atoms as PDB
        selected_atoms.save_pdb(pdb_filename)

    print(f"{len(os.listdir(output_dir))} PDB files saved in {output_dir}")

if __name__ == "__main__":
    input_dcd_file = 'D20.dcd'
    topology_file = 'ROS_1_RING.pdb'
    save_frames_as_pdb(input_dcd_file, topology_file, stride=50, output_dir='pdb_files')
