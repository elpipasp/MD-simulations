from Bio.PDB import PDBParser
import numpy as np

#define function to calculate distance between two points
def calculate_distance(point1, point2):
    return np.linalg.norm(point1 - point2)

#define function to calculate the angle between two vectors
def calculate_angle(vector1, vector2):
    unit_vector1 = vector1 / np.linalg.norm(vector1)
    unit_vector2 = vector2 / np.linalg.norm(vector2)
    dot_product = np.dot(unit_vector1, unit_vector2)
    angle = np.arccos(dot_product) * (180 / np.pi)
    return angle

#define function to calculate the normal vector of a plane defined by three points
def calculate_normal_vector(point1, point2, point3):
    vector1 = point2 - point1
    vector2 = point3 - point1
    normal_vector = np.cross(vector1, vector2)
    return normal_vector

#define function to determine the type of pi-pi interaction based on orientation and distance
def determine_pi_interaction_type(distance, angle):
    if 4.5 <= distance <= 7.5 and 60.0 <= angle <= 95.0:
        return "T-shaped"
    elif 4.5 <= distance <= 7.5 and 0.0 <= angle <= 30.0:
        return "Parallel-displaced"

#define function to determine cation-pi interaction based on distance
def determine_cation_pi_interaction(distance):
    if distance <= 5.5: 
        return True
    return False

#load the protein structure using Biopython's PDBParser=proper pdb files only
parser = PDBParser()
structure = parser.get_structure("protein", "gnrh_receptor_combined_standard.pdb")

#iterate through the structure to identify and analyse pi-pi and cation-pi
for model in structure:
    for chain in model:
        for residue in chain:
            if residue.get_id()[0] == " " and residue.get_resname() in ["PHE", "TYR", "TRP", "ARG", "LYS"]:
                atoms = residue.get_atoms()
                coordinates = np.array([atom.get_coord() for atom in atoms])
                if residue.get_resname() in ["PHE", "TYR", "TRP"]:
                    centroid = np.mean(coordinates, axis=0)
                    normal_vector1 = calculate_normal_vector(coordinates[0], coordinates[1], coordinates[2])  #calculate normal vector for current residue
                elif residue.get_resname() in ["ARG", "LYS"]:
                    cation_position = np.mean(coordinates, axis=0)
                
                for other_chain in model:
                    for other_residue in other_chain:
                        if other_residue.get_id()[0] == " " and other_residue.get_resname() in ["PHE", "TYR", "TRP", "ARG", "LYS"]:
                            other_atoms = other_residue.get_atoms()
                            other_coordinates = np.array([atom.get_coord() for atom in other_atoms])
                            if other_residue.get_resname() in ["PHE", "TYR", "TRP"]:
                                other_centroid = np.mean(other_coordinates, axis=0)
                            elif other_residue.get_resname() in ["ARG", "LYS"]:
                                other_cation_position = np.mean(other_coordinates, axis=0)  #calculate cation position for other residue
                            
                            if residue.get_resname() in ["PHE", "TYR", "TRP"] and other_residue.get_resname() in ["PHE", "TYR", "TRP"]:
                                distance = calculate_distance(centroid, other_centroid)
                                centroid_vector = other_centroid - centroid
                                angle = calculate_angle(centroid_vector, normal_vector1)  #calculate angle using the normal vector of the current residue
                                interaction_type = determine_pi_interaction_type(distance, angle)
                                if interaction_type is not None:
                                    print(f"{interaction_type} pi-pi interaction found between {residue.get_id()[1]}{residue.get_id()[2]} and {other_residue.get_id()[1]}{other_residue.get_id()[2]}")
                            
                            elif residue.get_resname() in ["PHE", "TYR", "TRP"] and other_residue.get_resname() in ["ARG", "LYS"]:
                                distance = calculate_distance(centroid, other_cation_position)
                                if determine_cation_pi_interaction(distance):
                                    print(f"Cation-pi interaction found between {residue.get_id()[1]}{residue.get_id()[2]} and {other_residue.get_id()[1]}{other_residue.get_id()[2]}")
