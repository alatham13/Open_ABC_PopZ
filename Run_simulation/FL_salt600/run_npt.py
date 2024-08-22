# run simulation of a protein.
import numpy as np
import pandas as pd
import os
import sys
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit

sys.path.insert(0, "/wynton/home/sali/aplatham/Programs/OpenABC/")
from openabc3.forcefields.parsers import MOFFParser
from openabc3.forcefields.moff_mrg_model import MOFFMRGModel
from openabc3.forcefields import functional_terms
from openabc3.lib import _amino_acids, _kcal_to_kj
from openabc3.utils.insert import insert_molecules

"""
Define a new force field that use nonbonded terms from HPS while other terms from SMOG.

"""

__location__ = '/wynton/home/sali/aplatham/Programs/OpenABC/openabc3/forcefields' # path of openabc/forcefields

class SMOGHPSModel(MOFFMRGModel):
    # rewrite add_contacts, deprecate add_elec_switch, replace with add_elec
    # for add_contacts, only protein is supported, since RNA is not supported by SMOG
    def add_contacts(self, hydropathy_scale='Urry', epsilon=0.2 * _kcal_to_kj, mu=1, delta=0.08, force_group=2):
        """
        Add nonbonded contacts.

        The raw hydropathy scale is scaled and shifted by: mu*lambda - delta

        Parameters
        ----------
        hydropathy_scale : str
            Hydropathy scale, can be KR or Urry.

        epsilon : float or int
            Contact strength.

        mu : float or int
            Hydropathy scale factor.

        delta : float or int
            Hydropathy shift factor.

        force_group : int
            Force group.

        """
        print('Add nonbonded contacts.')
        resname_list = self.atoms['resname'].tolist()
        atom_types = [_amino_acids.index(x) for x in resname_list]
        if hydropathy_scale == 'KR':
            print('Use KR hydropathy scale.')
            df_contact_parameters = pd.read_csv(f'{__location__}/parameters/HPS_KR_parameters.csv')
        elif hydropathy_scale == 'Urry':
            print('Use Urry hydropathy scale.')
            df_contact_parameters = pd.read_csv(f'{__location__}/parameters/HPS_Urry_parameters.csv')
        else:
            sys.exit(f'Error: hydropathy scale {hydropathy_scale} cannot be recognized!')
        sigma_ah_map, lambda_ah_map = np.zeros((20, 20)), np.zeros((20, 20))
        for i, row in df_contact_parameters.iterrows():
            atom_type1 = _amino_acids.index(row['atom_type1'])
            atom_type2 = _amino_acids.index(row['atom_type2'])
            sigma_ah_map[atom_type1, atom_type2] = row['sigma']
            sigma_ah_map[atom_type2, atom_type1] = row['sigma']
            lambda_ah_map[atom_type1, atom_type2] = row['lambda']
            lambda_ah_map[atom_type2, atom_type1] = row['lambda']
        print(f'Scale factor mu = {mu} and shift delta = {delta}.')
        lambda_ah_map = mu * lambda_ah_map - delta
        force = functional_terms.ashbaugh_hatch_term(atom_types, self.exclusions, self.use_pbc, epsilon,
                                                     sigma_ah_map, lambda_ah_map, force_group)
        self.system.addForce(force)

    def add_dh_elec(self, ldby=(1/1.26)*unit.nanometer, dielectric_water=80.0, cutoff=3.5*unit.nanometer,
                    force_group=3):
        """
        Add Debye-Huckel electrostatic interactions.

        Parameters
        ----------
        ldby : Quantity
            Debye length.

        dielectric_water : float or int
            Dielectric constant of water.

        cutoff : Quantity
            Cutoff distance.

        force_group : int
            Force group.

        """
        print('Add Debye-Huckel electrostatic interactions.')
        print(f'Set Debye length as {ldby.value_in_unit(unit.nanometer)} nm.')
        print(f'Set water dielectric as {dielectric_water}.')
        charges = self.atoms['charge'].tolist()
        # convert charges from MOFF to HPS
        for i in range(len(charges)):
            if charges[i]==0.25:
                charges[i]=0.5
        force = functional_terms.dh_elec_term(charges, self.exclusions, self.use_pbc, ldby, dielectric_water,
                                              cutoff, force_group)
        self.system.addForce(force)

    def add_elec_switch(self):
        # deprecate
        return None

tot_t=10000000
dt=int(tot_t/100)

# set simulation platform
platform_name = 'CUDA'

print('Successful import')

# start from predicted PDB
PopZ_parser = MOFFParser.from_atomistic_pdb('PopZ_FL_trimer_AF.pdb', 'PopZ_FL_trimer_CA.pdb')

# PopZ_FL_trimer: we want to model the complex by rigidizing H3-H4 as a complex, and keeping H1 and H2 as a flexible helix.
# Note: indexing goes across all residues, regardless of chain and starts at 0
old_native_pairs = PopZ_parser.native_pairs.copy()
new_native_pairs = pd.DataFrame(columns=old_native_pairs.columns)
# full length of an individual monomer
Nres=177
# keep interactions between H3-H4 across chains
H3H4_1=np.arange(135,171)
print(H3H4_1)
H3H4_2=H3H4_1+Nres
H3H4_3=H3H4_1+Nres*2
H3H4_tot=np.concatenate((H3H4_1,H3H4_2,H3H4_3))
# keep interactions within H2 on the same chain
H2_1=np.arange(101,128)
H2_2=H2_1+Nres
H2_3=H2_1+Nres*2
# keep interactions within H1 on the same chain
H1_1=np.arange(9,24)
print(H1_1)
H1_2=H1_1+Nres
print(H1_2)
H1_3=H1_1+Nres*2
print(H1_3)
for i, row in old_native_pairs.iterrows():
    a1, a2 = int(row['a1']), int(row['a2'])
    if a1 > a2:
        a1, a2 = a2, a1
    flag1=( (a1 in H2_1) and (a2 in H2_1) )
    flag2=( (a1 in H2_2) and (a2 in H2_2) )
    flag3=( (a1 in H2_3) and (a2 in H2_3) )
    flag4=( (a1 in H3H4_tot) and (a2 in H3H4_tot) )
    flag5=( (a1 in H1_1) and (a2 in H1_1) )
    flag6=( (a1 in H1_2) and (a2 in H1_2) )
    flag7=( (a1 in H1_3) and (a2 in H1_3) )
    if flag1 or flag2 or flag3 or flag4 or flag5 or flag6 or flag7:
        new_native_pairs.loc[len(new_native_pairs.index)] = row
PopZ_parser.native_pairs = new_native_pairs
PopZ_parser.parse_exclusions() # update exclusions based on the new native pairs

# setup initial topology
n_mol = 100
box_a=100
box_b=100
box_c=100
insert_molecules('PopZ_FL_trimer_CA.pdb', 'start.pdb', n_mol, box=[box_a, box_b, box_c])


# set up protein model
condensate = SMOGHPSModel()
# add protein to model
for i in range(n_mol):
    # append multiple FG parser instances
    condensate.append_mol(PopZ_parser)
# setup initial topology
top = app.PDBFile('start.pdb').getTopology()
condensate.create_system(top, box_a=box_a, box_b=box_b, box_c=box_c, remove_cmmotion=True)
# setup variables necessary for simulations
T=150
salt_concentration = 600*0.001
eps_w=80
eps_0=8.8541878128*(10**-12)
R=8.31446261815324
F=9.64853321233100184*(10**4)
# Caclulate and convert from SI units to nm
kappa=np.sqrt((eps_w*eps_0*R*300)/(2*(10**3)*(F**2)*salt_concentration) )*(10**9)
Debye_L=kappa*unit.nanometer
temperature = T*unit.kelvin
# add force terms to simulation
condensate.add_protein_bonds(force_group=1)
condensate.add_protein_angles(force_group=2,clip=True)
condensate.add_protein_dihedrals(force_group=3)
condensate.add_native_pairs(force_group=4)
condensate.add_contacts(force_group=5)
condensate.add_dh_elec(Debye_L, force_group=6)
pressure = 1*unit.bar
condensate.system.addForce(mm.MonteCarloBarostat(pressure, temperature))
# setup integrator
friction_coeff = 1/unit.picosecond
timestep = 10*unit.femtosecond
integrator = mm.LangevinMiddleIntegrator(temperature, friction_coeff, timestep)
# setup initial coordinates
init_coord = app.PDBFile('start.pdb').getPositions()
# setup simulation
condensate.set_simulation(integrator, platform_name, init_coord=init_coord, properties={'Precision': 'mixed'})
# perform energy minimization
condensate.simulation.minimizeEnergy()
output_interval = dt
output_dcd = 'npt.dcd'
# Run short simulation
condensate.add_reporters(output_interval, output_dcd)
condensate.simulation.reporters.append(app.checkpointreporter.CheckpointReporter('npt.cpt',output_interval))
condensate.simulation.reporters.append(app.StateDataReporter('npt.csv', output_interval, step=True, potentialEnergy=True, temperature=True, volume=True, density=True, speed=True, time=True))
condensate.simulation.context.setVelocitiesToTemperature(temperature)
condensate.simulation.step(tot_t)

# print final box vectors
state = condensate.simulation.context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, getParameters=True, enforcePeriodicBox=True)
with open('npt.xml', 'w') as f:
    f.write(mm.XmlSerializer.serialize(state))
# save final system
condensate.save_system('PopZ_FL_trimer_system.xml')
