# Function to inter and intra chain contact maps from slab simulations
# Written by Andrew Latham
# Note the inputs at top. These need to be customized for the system.
# Here, we assume 100 rigidbodies, with 2000 equilibration frames, and a bead-wise cutoff of 10 Ang

# Output files:
# inter_map.txt - inter-molecular contact map from slab simulation. Ordered according to protein sequence
# intra_map.txt - intra-molecular contact map from slab simulation. Ordered according to protein sequence

import sys
import os
import math
import numpy
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
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

# wraps coordinate back into the simulation box
def wrap_coord(pos,side):
    while pos<0:
        pos=pos+side
    while pos>side/2:
        pos=pos-side
    return pos

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


def contact_mat(dcd_file,pdb_file,mass_xml,npt_xml,nchain,start,cutoff):
    u1 = mda.Universe(pdb_file,dcd_file)

    # caclulate number of atoms per protein chain and number of timesteps
    protien=u1.select_atoms("all")
    print(protien.atoms)
    N1=int(len(protien)/nchain)
    print('Number of protein atoms per chain: '+str(N1))
    timesteps=len(u1.trajectory)-(start)
    print('Number of timesteps: '+str(timesteps))
    N = len(u1.atoms)

    # set mass of atoms
    atom_masses=find_mass(mass_xml)
    for i in range(N):
        u1.atoms[i].mass=atom_masses[i]

    # read in box dimensions
    box=find_box(npt_xml)

    # counter for number of timesteps
    count=0
    inter_map=numpy.zeros((N1,N1))
    intra_map=numpy.zeros((N1,N1))


    for ts in u1.trajectory:
        print('frame: ' + str(ts.frame))
        if ts.frame < start:
            pass
        else:
            # loop over all proteins
            for i in range(0, nchain):
                index1 = i * N1
                index2 = ((i + 1) * N1)
                atoms1 = u1.atoms[index1:index2]
                # loop over proteins again
                for j in range(i,nchain):
                    index3 = j * N1
                    index4 = ((j + 1) * N1)
                    atoms2 = u1.atoms[index3:index4]
                    # calculate distance between atoms in first protein and atoms in 2nd protein
                    dist = distance_array(atoms1.positions, atoms2.positions, box)
                    # convert distance to contact matrix
                    contacts=numpy.where(dist<cutoff,1.0,0.0)
                    # if i==j, contacts are within the same chain (intra)
                    if i==j:
                        intra_map+=contacts
                    # otherwise, they are between chains
                    else:
                        inter_map+=contacts
            # keep track of number of timesteps
            count=count+1

    # normalize contact maps by number of proteins
    intra_map=intra_map/(count*nchain)
    inter_map = inter_map / (count * nchain)

    return inter_map,intra_map


# run main function
inter_map,intra_map=contact_mat(DCD,PDB,SYSTEM,NPT,nchain,eq,cut)

# save outputs
numpy.savetxt('inter_map.txt',inter_map)
numpy.savetxt('intra_map.txt',intra_map)

