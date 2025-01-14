### First steps of protein-ligand structure analysis 

# Libraries gemmi and numpy are required for manipulation of pdb files 
import gemmi
import numpy as np 

# Import PDB file 
oxa48_pdb = gemmi.read_structure("4s2n.pdb")

avibactam = "NXL"  # Define ligand name with the three letter code annotated in the structure

model = oxa48_pdb[0]  # First model in the PDB
ligand_residue = None
ligand_positions = [] # List to store ligand atom positions
ligands = [] # List to store ligand name and residue

# Search for the ligand across all chains, looping over chains of interest (here, chains A and B)
for chain in model:
    if chain.name in ['A', 'B']:
        for residue in chain:
            if residue.name == avibactam:
                ligand_residue = residue
                ligand_positions = [atom.pos for atom in ligand_residue]
                print(f"Ligand {avibactam} found in chain {chain.name} at positions: {ligand_positions}")
                ligands.append((chain.name, residue)) # Store the ligand name and residue 

# Print text confirming finding of ligands 
if ligand_residue:
    print(f"Ligand {avibactam} successfully located.")
else:
    print(f"Ligand {avibactam} not found in any chain.")
# The output will also contain the coordinates for the atoms in each ligand (inclusion of atom.pos)

# Determine the correct water annotation 
water_annotation = None

# Iterate through chains and residues to check for waters as either HOH or WAT 
for chain in model:
    if chain.name in ['A', 'B']:
        for residue in chain:
            if residue.name in ["HOH", "WAT"]:
                water_annotation = residue.name
                break  
        if water_annotation:
            break  

# Print text confirming water annotation as HOH or WAT
if water_annotation:
    print(f"Water molecules are annotated as {water_annotation}")


# We require the matplotlib.pyplot library to generate simple plots
import matplotlib.pyplot as plt
# We also require the counter within the collections library
from collections import Counter

# Initialize counters for individual amino acids and interaction types
amino_acid_counts_A = Counter()
amino_acid_counts_B = Counter()

interaction_counts_A = Counter({
    "Hydrophobic": 0,
    "Polar": 0,
    "Positively charged side chain": 0,
    "Negatively charged side chain": 0,
    "Water": 0
})

interaction_counts_B = Counter({
    "Hydrophobic": 0,
    "Polar": 0,
    "Positively charged side chain": 0,
    "Negatively charged side chain": 0,
    "Water": 0
})

# Define side-chain property groups
side_chain_groups = {
    "Hydrophobic": ["ALA", "VAL", "LEU", "ILE"],
    "Polar": ["THR", "SER", "TYR", "TRP"],
    "Positively charged side chain": ["ARG", "LYS", "HIS"],
    "Negatively charged side chain": ["ASP", "GLU"]
}

# Loop over each ligand and classify interactions
for ligand_index, (ligand_chain, ligand_residue) in enumerate(ligands):
    print(f"\nAnalysing all possible interactions for avibactam in chain {ligand_chain}:")
    ligand_positions = [atom.pos for atom in ligand_residue] 

    for chain in model:
        if chain.name in ['A', 'B']:  
            for residue in chain:
                for atom in residue:
                    for ligand_atom in ligand_residue:
                        distance = np.linalg.norm(
                            np.array([ligand_atom.pos.x, ligand_atom.pos.y, ligand_atom.pos.z]) -
                            np.array([atom.pos.x, atom.pos.y, atom.pos.z])
                        )
                        if distance <= 3.5:  # Define a cutoff distance for interactios in Angstroms
                            residue_number = residue.seqid.num
                            residue_name = residue.name
                            

                        # Update amino acid counts
                            if ligand_residue.name == "NXL":
                                if chain.name == 'A':
                                    amino_acid_counts_A[residue_name] += 1
                                elif chain.name == 'B':
                                    amino_acid_counts_B[residue_name] += 1

                            # Update interaction counts based on property groups 
                            if residue.name in ["ALA", "VAL", "LEU", "ILE"]:
                                print(f"Hydrophobic: {residue.name} {residue_number} - {avibactam} ({ligand_atom.name}) ({distance:.2f} Å)")
                                if chain.name == 'A':
                                    interaction_counts_A["Hydrophobic"] += 1
                                elif chain.name == 'B':
                                    interaction_counts_B["Hydrophobic"] += 1
                            elif residue.name in ["THR", "SER", "TYR", "TRP"]:
                                print(f"Polar: {residue.name} {residue_number} - {avibactam} ({ligand_atom.name}) ({distance:.2f} Å)")
                                if chain.name == 'A':
                                    interaction_counts_A["Polar"] += 1
                                elif chain.name == 'B':
                                    interaction_counts_B["Polar"] += 1
                            elif residue.name in ["ARG", "LYS", "HIS"]:
                                print(f"Positively charged side chain: {residue.name} {residue_number} - {avibactam} ({ligand_atom.name}) ({distance:.2f} Å)")
                                if chain.name == 'A':
                                    interaction_counts_A["Positively charged side chain"] += 1
                                elif chain.name == 'B':
                                    interaction_counts_B["Positively charged side chain"] += 1
                            elif residue.name in ["ASP", "GLU"]:
                                print(f"Negatively charged side chain: {residue.name} {residue_number} - {avibactam} ({ligand_atom.name}) ({distance:.2f} Å)")
                                if chain.name == 'A':
                                    interaction_counts_A["Negatively charged side chain"] += 1
                                elif chain.name == 'B':
                                    interaction_counts_B["Negatively charged side chain"] += 1
                            
                        
                            # Classify residues based on properties
                            for group, amino_acids in side_chain_groups.items():
                                if residue_name in amino_acids:
                                    if chain.name == 'A':
                                        interaction_counts_A[group] += 1
                                    elif chain.name == 'B':
                                       interaction_counts_B[group] += 1
                            if residue_name == "HOH":
                                if chain.name =='A' or 'B':
                                    interaction_counts_A["Water"] += 1
                                    interaction_counts_B["Water"] += 1

    # Classify surrounding waters
    for chain in model:
        for residue in chain:
            # Check if residue is a water molecule (e.g., "HOH" or "WAT")
            if residue.name in ["HOH"]:
                for atom in residue:
                    # For water molecules, check interactions with the ligand
                    for ligand_atom in ligand_residue:
                        distance = np.linalg.norm(
                            np.array([ligand_atom.pos.x, ligand_atom.pos.y, ligand_atom.pos.z]) -
                            np.array([atom.pos.x, atom.pos.y, atom.pos.z])
                        )
                        if distance <= 3.5:  # Define a cutoff distance
                            # Get residue number and name
                            residue_number = residue.seqid.num
                            residue_name = residue.name
                            # Classify interaction
                            print(f"Water interaction: {residue_name} {residue_number} - {avibactam} ({ligand_atom.name}) ({distance:.2f} Å)")
                           

# Function to plot individual amino acid counts
def plot_amino_acid_counts(amino_acid_counts, chain_name):
    # Filter out the ligand atoms from the counts to just include amino acids 
    filtered_counts = {key: value for key, value in amino_acid_counts.items() if key != "NXL"}
    labels = list(filtered_counts.keys())
    values = list(filtered_counts.values())

    colors = plt.cm.tab20(np.linspace(0, 1, len(labels)))

    plt.figure(figsize=(10, 6))
    plt.bar(labels, values, color=colors)
    plt.xlabel('Residues')
    plt.ylabel('Count')
    plt.title(f'Potential protein-ligand interactions in Chain {chain_name} (\u2264 3.5 \u00C5)') ##using unicode to include proper notations in title
    plt.xticks(rotation=45)
    plt.tight_layout()

# Function to plot side chain properties of interacting amino acids
def plot_interaction_counts(interaction_counts, chain_name):
    labels = list(interaction_counts.keys())
    values = list(interaction_counts.values())

    # Plot the bar chart
    plt.figure(figsize=(8, 6))
    plt.bar(labels, values, color=['blue', 'green', 'orange', 'red', 'purple'])
    plt.xlabel('Residue properties')
    plt.ylabel('Count')  
    plt.title(f'Properties of ligand interacting residues and waters in Chain {chain_name} (\u2264 3.5 \u00C5)')
    plt.xticks(rotation=45)
    plt.tight_layout()


# Plot amino acid counts
plot_amino_acid_counts(amino_acid_counts_A, 'A')
plot_amino_acid_counts(amino_acid_counts_B, 'B')

# Plot grouped interaction counts
plot_interaction_counts(interaction_counts_A, 'A')
plot_interaction_counts(interaction_counts_B, 'B')

# Show all plots
plt.show()