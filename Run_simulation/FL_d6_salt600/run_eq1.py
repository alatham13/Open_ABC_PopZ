# run simulation of a protein.
import numpy as np
import pandas as pd
import os
import sys
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit

tot_t=10000000
nsteps=100
dt=int(tot_t/nsteps)
box_Z=np.asarray([0.0,0.0,500.0])

Tstart=150
Tend=300
deltaT=(Tend-Tstart)/nsteps


# set simulation platform
platform_name = 'CUDA'
properties={'CudaPrecision':'mixed'}

print('Successful import')

# Load system
system_xml = 'PopZ_FL_trimer_d6_system.xml'
with open(system_xml, 'r') as f:
    system = mm.XmlSerializer.deserialize(f.read())
# remove barostat to run nvt simulation. Assumes the barostat is the last force
system.removeForce(system.getNumForces()-1)

top = app.PDBFile('start.pdb').getTopology()

# initiate simulation
# setup integrator
friction_coeff = 0.01/unit.picosecond
timestep = 10*unit.femtosecond
temperature = Tstart*unit.kelvin
integrator = mm.LangevinMiddleIntegrator(temperature, friction_coeff, timestep)
platform = mm.Platform.getPlatformByName(platform_name)
# begin simulation
simulation = app.Simulation(top, system, integrator, platform, properties)
# Load checkpoint file
simulation.loadCheckpoint('npt.cpt')
state=simulation.context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, getParameters=True, enforcePeriodicBox=True)
# change box vectors
box_vec = state.getPeriodicBoxVectors(asNumpy=True)
simulation.context.setPeriodicBoxVectors(box_vec[0], box_vec[1], box_Z)
system=simulation.context.getSystem()
# save new system to xml
with open('PopZ_FL_trimer_d6.xml', 'w') as output:
    output.write(mm.XmlSerializer.serialize(system))
# add dcd. Append if available
output_dcd = 'eq1.dcd'
simulation.reporters.append(app.DCDReporter(output_dcd, dt))
# Append state reporter
simulation.reporters.append(app.StateDataReporter('eq1.csv', dt, step=True, time=True, potentialEnergy=True, temperature=True, volume=True, density=True, speed=True))
# Append checkpoint file
simulation.reporters.append(app.CheckpointReporter('eq1.cpt', dt))
print('Running Production...')
for i in range(nsteps):
    integrator.setTemperature((Tstart + deltaT*(i+1))*unit.kelvin)
    simulation.step(dt)
    state=simulation.context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, getParameters=True, enforcePeriodicBox=True)
    box_vec = state.getPeriodicBoxVectors(asNumpy=True)
    print(box_vec)
