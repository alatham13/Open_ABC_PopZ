# run simulation of a protein.
import numpy as np
import pandas as pd
import os
import sys
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit

# total length the simulation in number of steps (including previous 20,000,000 equibration steps)
tot_t=520000000
nframes=5000
dt=int(100000)
print('dt is: '+str(dt))

Temperature=300

out_dcd='prod.dcd'


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
simulation = app.Simulation(top, system, integrator, platform,properties)

restart = False
time_so_far = 0
# check if checkpoint exists
if os.path.isfile('prod.cpt'):
    restart = True
    print('Loading a restart...')
    simulation.loadCheckpoint('prod.cpt')
    #simulation.loadState('checkpnt.xml')
    # find out how many steps to run left
    sim_log = [ i for i in open('prod.csv').readlines() ]
    last_line = sim_log[-1]
    time_so_far = float(last_line.split(',')[1]) # how much time run so far
    tot_t -= int(round(time_so_far) / (10/1000))
    print('number of remaining steps: '+str(tot_t))
    if tot_t==0:
        print('Simulation complete. Closing...')
        exit()
else:
    print('Checkpoint file not found!!! Closing...')
    exit()
# add dcd. Append if available
simulation.reporters.append(app.DCDReporter(out_dcd, dt, append=restart))
# Append state reporter
simulation.reporters.append(app.StateDataReporter('prod.csv', dt, step=True, time=True, potentialEnergy=True, temperature=True, volume=True, density=True, speed=True))
# Append checkpoint file
simulation.reporters.append(app.CheckpointReporter('prod.cpt', dt))
print('Running Production...')
simulation.step(int(tot_t))
