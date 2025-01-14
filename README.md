# Elementary-protein-ligand-structure-analysis

READ.ME file for first steps to analysis of protein-ligand structures
This is a pipeline written in Python to enable elementary analysis of any protein-ligand structures as long as a pdb file has been generated from your own macro-crystallography experiments and the three-letter code for the ligand is known. 
To run:
1.	Download the Python script as firststepstoanalysisofprotein-ligandstructures.py in this repository along with the example pdb file of OXA-48 β-lactamase bound to avibactam (PDB: 4S2N) into a single working directory. 
2.	The libraries gemmi and numpy need importing before any analysis of imported pdb files can be initiated.  
3.	The script can be edited with your pdb file of choice and the ligand name must be also edited in line. The number of protein chains to be analysed can be edited if intra-protein binding site differences are being investigated. If only a specific few amino acids are of interest, the amino acid count can be edit easily 
4.	The maximum distance governing significant interactions was deemed to be 3.5 Å, however this can easily be adjusted to define a larger or smaller cut-off distance for the counts. 
5.	The script will generate both a list of the interactions in truncated text form which can be viewed fully as a scrollable element or in a text editor, as well as graphs illustrating amino acid name and property counts. 

