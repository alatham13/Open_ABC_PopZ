# run simulation of a protein.
import numpy as np
import pandas as pd
import os
import sys
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit

tot_t=500000000
nframes=5000
dt=int(tot_t/nframes)

Temperature=300

# set simulation platform
platform_name = 'CUDA'
properties={'CudaPrecision':'mixed'}

print('Successful import')

system_xml = 'PopZ_FL_trimer.xml'
with open(system_xml, 'r') as f:
    system = mm.XmlSerializer.deserialize(f.read())

top = app.PDBFile('start.pdb').getTopology()

# initiate simulation
# setup integrator
friction_coeff = 0.01/unit.picosecond
timestep = 10*unit.femtosecond
temperature = Temperature*unit.kelvin
integrator = mm.LangevinMiddleIntegrator(temperature, friction_coeff, timestep)
platform = mm.Platform.getPlatformByName(platform_name)
# begin simulation
simulation = app.Simulation(top, system, integrator, platform,properties)
# Load checkpoint file
simulation.loadCheckpoint('eq1.cpt')
# add dcd. Append if available
output_dcd = 'prod.dcd'
simulation.reporters.append(app.DCDReporter(output_dcd, dt))
# Append state reporter
simulation.reporters.append(app.StateDataReporter('prod.csv', dt, step=True, time=True, potentialEnergy=True, temperature=True, volume=True, density=True, speed=True))
# Append checkpoint file
simulation.reporters.append(app.CheckpointReporter('prod.cpt', dt))
print('Running Production...')
simulation.step(tot_t)
