# Code to simulate and analyze PopZ condensates
By Andrew P. Latham
If you use this code, please cite:

Scholl D, Boyd T, Latham AP, Salazar A, Khan A, Boeynaems S, Holehouse AS, Lander GC, Sali A, Park D, Deniz AA, Lasker K Cellular Function of a Biomolecular Condensate Is Determined by Its Ultrastructure. bioRxiv 2024, DOI:10.1101/2024.12.27.630454.

Liu S, Wang C, Latham AP, Ding X, Zhang B OpenABC enables flexible, simplified, and efficient GPU accelerated simulations of biomolecular condensates. PLoS Comput Biol 2023, 19, e1011442.

Regy RM, Thompson J, Kim YC, Mittal J Improved coarse-grained model for studying sequence dependent phase separation of disordered proteins. Protein Sci 2021, 30, 1371.

Here, we walk through an example of the simulations and analysis used in our manuscript.
All steps for performing simulations and analysis are demonstrated using the WT PopZ system, 
and are implemented in OpenABC (https://github.com/ZhangGroup-MITChemistry/OpenABC/), version 1.06.

## Running HPS+SMOG simulations -------------------------
Slab simulations of PopZ are run through a 4 step procedure. They begin from an alpha-fold structure of the PopZ mutant, in a set oligomerization state. The included scripts are examples for WT-PopZ trimer at 150 mM monovalent salt, WT-PopZ trimer at 600 mM monovalent salt, WT-PopZ trimer at 150 mM monovalent salt and pH~4, and OD-PopZ trimer at 150 mM monovalent salt. 

Dependencies: numpy, pandas, OpenMM, OpenABC (https://github.com/ZhangGroup-MITChemistry/OpenABC)

1. run_npt.py - initializes the system and runs an npt simulation at low temperature (150 K), producing a dense phase

(Note: for different PopZ mutants, the PopZ_parser variable needs to be appropriately changed to define oligomerization state and the ordered regions)

(Note 2: the __location__ variable must be set to the openabc/forcefields folder)

2. run_eq.py - starts by expanding the Z-axis to prepare for slab simulations. In this new box, it runs an annealing simulation, warming the condensate from 150 K to 300 K.

3. run_prod.py - initializes the production simulation. This simulation lasts for 5 micro-seconds at nvt in the slab configuration.

4. restart_prod.py - finishes the production simulation by looping over the number of steps until  5 micro-seconds is reached.

 
## Analyzing simulations of PopZ condensates -------------
Simulations are analyzed for a variety of quantities. We include the scripts for producing these quantities as well as a sample output for WT-PopZ trimer at 150 mM monovalent salt.

Dependencies: numpy, MDAnalysis, xml

Inputs / Parameters (same for all scripts, and are input manually in the first few lines):

DCD - dcd file used to run the simulation

PDB - pdb file used to set the topology

SYSTEM - xml file used to define the topology system. Input as sys.argv[1]

NPT -  xml file used to define box lengths

nchain - number of protein chains in the simulation box (this should be the number of independently moving oligomers)

eq - number of equilibration frames

cut - distance cutoff for interactions (in Angstroms)

1. cluster_PopZ_mindist.py - Function to calculate density and cluster size from slab simulations

Outputs:

cluster_size_mindist.txt - size of the largest 4 clusters, in number of molecules

protein_hist_mindist.txt - histogram of protein density, centered at the largest cluster

(Note: The protein density needs to be converted into mg/mL and normalized by the number of timesteps / box size. See Z_dist_cluster_size.m for an example of final data analysis. Note that the xy variable needs to be manually set to the length of the x/y axis for each simulation.)

2. contact_map_PopZ.py - Function to inter and intra chain contact maps from slab simulations

Outputs:

inter_map.txt - inter-molecular contact map from slab simulation. Ordered according to protein sequence

intra_map.txt - intra-molecular contact map from slab simulation. Ordered according to protein sequence

(Note: These contact maps need to be appropriately normalized before plotting. This is done by plot_contacts_log2.m, which normalizes and plots the contact map, with intra molecular contacts on the upper diagonal and intermolecular contacts on the bottom diagonal)
