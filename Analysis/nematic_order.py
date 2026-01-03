# Function to calculate the nematic ordrer parameter from slab simulations of PopZ trimers
# Written by Andrew Latham
# Note the inputs at top. These need to be customized for the system.
# Here, we assume 100 molecules, with 2000 equilibration frames, and a COM cutoff for the local nematic order parameter of 50 Ang

# Output files:
# nematic_order.txt - nematic order parameter of the entire system
# nematic_director.txt - vector of the nematic ordering
# local_order - nematic order parameter calculated only on chains within a COM cutoff
# local_chains - average number of neighboring chains

import sys
import os
import math
import numpy
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
from MDAnalysis import transformations
import matplotlib.pyplot as plt
import freud
import xml.etree.ElementTree as ET

# Inputs
# dcd file used to run the simulation
DCD='prod.dcd'
# pdb file used to set the topology
PDB='start.pdb'
# xml file used to define the system
SYSTEM=sys.argv[1]
# xml file used to define the topology, specifically the box length
NPT='npt.xml'
# number of protein chains in the simulation box
nchain=100
# number of equilibration frames
eq=2000
# distance cutoff
cut=10
# cutoff for local nematic order parameter
local_cut=50




# function that runs depth first search algorithm on a graph
def dfs_iterative(graph, start):
    stack, path = [start], []
    while stack:
        vertex = stack.pop()
        if vertex in path:
            continue
        path.append(vertex)
        if vertex in graph:
            for neighbor in graph[vertex]:
                stack.append(neighbor)
    return path


def mat_to_dict(mat):
    dict1={}
    for i in range(0,len(mat)):
        for j in range(0,len(mat[i])):
            if i in dict1 and mat[i][j]==1:
                dict1[i].append(j)
            elif i not in dict1 and mat[i][j]==1:
                dict1[i] = [j]
    return dict1

def wrap_coord(pos,side):
    while pos<0:
        pos=pos+side
    while pos>side/2:
        pos=pos-side
    return pos

def wrap_coord_Z(ts,box,nchain,pos_com,XYZ=2):
    N = ts.n_atoms
    atoms_chain = int(N / nchain)
    print('Number of atoms per chain:')
    print(atoms_chain)
    # translate1 - places center of mass of largest cluster at 0
    for i in range(0,nchain):
        ts.positions[i*atoms_chain:(i+1)*atoms_chain,XYZ]=ts.positions[i*atoms_chain:(i+1)*atoms_chain,XYZ]-pos_com
        COG_pos=numpy.mean(ts.positions[i*atoms_chain:(i+1)*atoms_chain,XYZ])
        while COG_pos<0:
            # update atom positions
            ts.positions[i*atoms_chain:(i+1)*atoms_chain, XYZ]=ts.positions[i*atoms_chain:(i+1)*atoms_chain,XYZ]+box[XYZ]
            # update COG positions
            COG_pos=COG_pos+box[XYZ]
        while COG_pos>box[XYZ]/2:
            # update atom positions
            ts.positions[i*atoms_chain:(i+1)*atoms_chain, XYZ]=ts.positions[i*atoms_chain:(i+1)*atoms_chain,XYZ]-box[XYZ]
            # update COG positions
            COG_pos=COG_pos-box[XYZ]
    # translate2 - shifts the entire simulation by box / 2 (thus, the largest cluster is now centered at box/2)
    ts.positions[:, XYZ]=ts.positions[:,XYZ]+box[XYZ]/2
    return ts

# function to find the mass of particles from the system XML file
def find_mass(xml_file):
    # read in system XML file
    tree=ET.parse(xml_file)
    root=tree.getroot()
    # list for masses of each particle
    mass_list=[]
    # go over all branches of the XML
    for branch in root:
        # the branch that stores mass is named 'Particles'
        if branch.tag=='Particles':
            for leaf in branch:
                # append these masses to the mass_list
                mass_list.append(float(leaf.attrib['mass']))
    return mass_list

# function to find the X/Y box edges from the topology XML file
def find_box(xml_file,find_z=False,default_Z=5000):
    # read in system XML file
    tree=ET.parse(xml_file)
    root=tree.getroot()
    box=numpy.zeros(6)
    # go over all branches of the XML
    for branch in root:
        # the branch that stores box length is named 'PeriodicBoxVectors'
        if branch.tag=='PeriodicBoxVectors':
            # multiply by 10 to convert nm to Ang.
            box[0]=float(branch[0].attrib['x'])*10
            box[1] = float(branch[1].attrib['y'])*10
            # Z has not yet been resized. Save only if asked
            if find_z:
                box[2] = branch[2].attrib['z']
            else:
                box[2]=default_Z
    # set angles for box edges. This assumes a rectangular box
    box[3]=90
    box[4]=90
    box[5]=90
    return box


def nematic_order_parameter(dcd_file,pdb_file,mass_xml,npt_xml,nchain,start,cutoff,local_cutoff):
    u1 = mda.Universe(pdb_file,dcd_file)

    # caclulate number of atoms per protein chain and number of timesteps
    protien=u1.select_atoms("all")
    print(protien.atoms)
    N1=int(len(protien)/nchain)
    print('Number of protein atoms per chain: '+str(N1))
    timesteps=len(u1.trajectory)-(start)
    print('Number of timesteps: '+str(timesteps))
    N = len(u1.atoms)

    # Atoms to select. Assumes 3 chains. From the start of H3, resid 136 to the end of H4, resid 171
    sel1=135
    sel2=170

    # set mass of atoms
    atom_masses=find_mass(mass_xml)
    for i in range(N):
        u1.atoms[i].mass=atom_masses[i]

    # read in box dimensions
    box=find_box(npt_xml)

    # counter for number of timesteps, initial empty pas
    count=0
    path_tot_timestep = []
    Z_mat2 = numpy.zeros((timesteps, 4))
    # Empty output for nematic order parameter
    nematic_order=numpy.zeros((timesteps, 1))
    nematic_director=numpy.zeros((timesteps, 3))
    local_order=numpy.zeros((timesteps, nchain))
    local_chains=numpy.zeros((timesteps, 1))


    for ts in u1.trajectory:
        print('frame: ' + str(ts.frame))
        if ts.frame < start:
            pass
        else:
            # Calculate COM of the largest cluster -----------------------------------------------------------------------------------------
            com = numpy.zeros((nchain, 3))
            mass = numpy.zeros((nchain, 1))
            # convert to adjacency matrix
            mat = numpy.zeros((nchain, nchain))
            # calculate com of each protein chain
            for i in range(0, nchain-1):
                # select atoms in each chain
                index1 = i * N1
                index2 = ((i + 1) * N1)
                atoms1 = u1.atoms[index1:index2]
                com[i, :] = atoms1.center_of_mass()
                mass[i] = atoms1.total_mass()
                for j in range(i+1,nchain):
                    # find all pairs of chain
                    index3 = j * N1
                    index4 = ((j + 1) * N1)
                    atoms2 = u1.atoms[index3:index4]
                    # caclulate total distance matrix
                    dist2 = distance_array(atoms1.positions, atoms2.positions, box)
                    # use minimum of distance matrix to calculate contact matrix
                    dist = numpy.amin(dist2[:])
                    if dist < cutoff:
                        mat[i, j] = mat[i, j] + 1
                        mat[j, i] = mat[i, j]
            # need com/mass for all chains. Add for the last chain here
            chain_index=nchain-1
            index1 = chain_index * N1
            index2 = ((chain_index + 1) * N1)
            atoms1 = u1.atoms[index1:index2]
            com[chain_index, :] = atoms1.center_of_mass()
            mass[chain_index] = atoms1.total_mass()

            # Calculate dfs at each timestep from matrix
            graph2 = mat_to_dict(mat)
            visited = []
            path_tot = []
            for i in range(0, len(mat)):
                flag = 0
                for k in range(0, len(visited)):
                    if visited[k] == i:
                        flag = 1

                if flag == 0:
                    path1 = []
                    path1 = dfs_iterative(graph2, i)
                    for j in range(0, len(path1)):
                        visited.append(path1[j])
                    path_tot.append(path1)
                path_tot_timestep.append(path_tot)
                #print(path_tot_timestep)

            # find atoms in largest cluster
            n = len(path_tot)
            l2 = -1
            index = -1
            #print(path_tot)
            l_list=[]
            for i in range(0, n):
                l = len(path_tot[i])
                l_list.append(l)
                if l > l2:
                    # reset old variables
                    l2 = l
                    index = i
            atoms_in_cluser = path_tot[index]
            mass_tot = 0
            for i in range(0, len(atoms_in_cluser)):
                index = atoms_in_cluser[i]
                mass_tot = mass_tot + mass[index]
            com2 = (mass[:, 0] / mass_tot) * numpy.cos((com[:, 2] / box[2]) * 2 * numpy.pi)
            com3 = (mass[:, 0] / mass_tot) * numpy.sin((com[:, 2] / box[2]) * 2 * numpy.pi)

            # Calculate COM of cluster. Use angles to avoid issues with periodicity
            Z1 = 0
            Z2 = 0
            # l_list: list of cluster sizes
            l_list=numpy.asarray(l_list)
            l_list[::-1].sort()
            l_list.resize((4))
            Z_mat2[count,:]=l_list
            print(Z_mat2[count,:])
            for i in range(0, l2):
                index = atoms_in_cluser[i]
                Z1 = Z1 + com2[index]
                Z2 = Z2 + com3[index]
            Z1 = Z1 / l2
            Z2 = Z2 / l2
            theta = numpy.arctan2(-1 * Z2, -1 * Z1) + numpy.pi
            Z_com = (box[2] / (2 * numpy.pi)) * theta


            # Add transformation to set new postions
            if count==0:
                # Add box dimensions to the universe before writing
                transform = transformations.boxdimensions.set_dimensions(box)
                u1.trajectory.add_transformations(transform)
            # Wrap Z-coordinate according to the COM of the box
            ts = wrap_coord_Z(ts, box, nchain, Z_com)

            # Calculate the average nematic order parameter -----------------------------------------------------------------------------------------
            
            # Get positions and orientations of each PopZ trimer
            pos = numpy.zeros((nchain, 3))
            orientations = numpy.zeros((nchain, 3))
            for i in range(0, nchain):
                # select atoms in each chain
                index1 = i * N1
                index2 = ((i + 1) * N1)
                atoms1 = u1.atoms[index1:index2]
                pos[i, :] = atoms1.center_of_mass()
                # start of H3, resid 136
                atoms_start=u1.select_atoms(f'index {index1+sel1} or index {index1+sel1+int(N1/3)} or index {index1+sel1+int(2*N1/3)}')
                # end of H4, resid 171
                atoms_end=u1.select_atoms(f'index {index1+sel2} or index {index1+sel2+int(N1/3)} or index {index1+sel2+int(2*N1/3)}')
                # draw unit vector through the center of masses and save as the orientation
                vector=atoms_start.center_of_mass()-atoms_end.center_of_mass()
                orientations[i,:]=vector/numpy.linalg.norm(vector)
            # Compute and save nematic order prameter
            nematic = freud.order.Nematic()
            nematic.compute(orientations)
            nematic_order[count]=nematic.order
            nematic_director[count,:]=nematic.director

            if count==timesteps-1:
                # make a 3D plot of the system
                fig = plt.figure()
                ax = fig.add_subplot(111, projection="3d")
                # Plot as quiver3D. Put Z-axis (long axis) on X, X on Y, and Y on Z.
                # Normalize Z so that the middle of the box is at 0, and convert Ang to nm
                ax.quiver3D(
                    (pos[:, 2]-box[2]/2)/10,
                    (pos[:, 0])/10,
                    (pos[:, 1])/10,
                    orientations[:, 0],
                    orientations[:, 1],
                    orientations[:, 2],
                    length=4,
                    pivot='middle',
                    normalize=True,
                    arrow_length_ratio=0.3,
                    color="k"
                )
                ax.set_xlabel("Box length [nm]", fontname="Helvetica", fontsize=16)
                ax.set_ylabel("X Axis [nm]", fontname="Helvetica", fontsize=16)
                ax.set_zlabel("Y Axis [nm]", fontname="Helvetica", fontsize=16)
                fig.savefig('PopZ_orientation.png',transparent=True,format='png',dpi=1200)
                fig.savefig('PopZ_orientation.eps',transparent=True,format='eps',dpi=1200)
            
            # Calculate the local nematic order parameter -----------------------------------------------------------------------------------------
            com_dist_mat = distance_array(pos, pos, box)
            N_neighbors=[]
            # For each chain
            for i in range(nchain):
                orientations2 = []
                for j in range(nchain):
                    # find other chains within the cutoff, note the orientation of those chains
                    if com_dist_mat[i,j]<local_cutoff:
                        orientations2.append(orientations[j,:])
                orientations2=numpy.asarray(orientations2)
                # Calculate the order parameter if there are other chains nearby
                if len(orientations2)>1:
                    local = freud.order.Nematic()
                    local.compute(orientations2)
                    local_order[count,i]=local.order
                    N_neighbors.append(len(orientations2))
                # place nan as place holder
                else:
                    local_order[count,i]=numpy.nan
            local_chains[count]=numpy.mean(N_neighbors)
            
            count = count + 1
    
    # Save nematic order parameter and director to a new file
    numpy.savetxt('nematic_order.txt',nematic_order)
    numpy.savetxt('nematic_director.txt',nematic_director)
    # Save the local nematic order and the average number of chains used to calculate the local nematic order to files
    numpy.savetxt('local_order.txt',local_order)
    numpy.savetxt('local_chains.txt',local_chains)
    return



nematic_order_parameter(DCD,PDB,SYSTEM,NPT,nchain,eq,cut,local_cut)

