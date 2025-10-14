from netpyne import specs, sim
from neuron import h
import matplotlib.pyplot as plt
import numpy as np

###############################################################################
# 1. NetPyNE network and cell definition (just dummy “cells” to host voxels)
###############################################################################
netParams = specs.NetParams()

# simple spherical cell with a soma section
cellRule = {'conds': {'cellType': 'voxelhost', 'cellModel': 'HH'},
            'secs': {'soma': {'geom': {'L':10, 'diam':10}, 'mechs': {}}}}
netParams.cellParams['VoxelHost'] = cellRule

# population of 3 “cells” (each will host one voxel)
netParams.popParams['VoxelPop'] = {'cellType': 'voxelhost','cellModel': 'HH' , 'numCells': 3}

###############################################################################
# 2. Simulation configuration
###############################################################################
simConfig = specs.SimConfig()
simConfig.duration = 5000   # ms
simConfig.dt = 0.1
simConfig.verbose = False
simConfig.recordStep = 0.1

###############################################################################
# 3. Create network
###############################################################################
sim.initialize(simConfig=simConfig, netParams=netParams)  # create network object and set cfg and net params
sim.net.createPops()  # instantiate network populations
sim.net.createCells()  # instantiate network cells based on defined populations


cells = sim.net.cells

# This function makes each pointer for neighboring concentrations point
# to its own soma as a placeholder before being overwritten later.
# This prevents this voxels at end points and corners from pointing outside of the voxel lattice
def init_voxel(vox):
    h.setpointer(vox._ref_conc, 'conc_xp', vox)
    h.setpointer(vox._ref_conc, 'conc_xn', vox)
    h.setpointer(vox._ref_conc, 'conc_yp', vox)
    h.setpointer(vox._ref_conc, 'conc_yn', vox)
    h.setpointer(vox._ref_conc, 'conc_zp', vox)
    h.setpointer(vox._ref_conc, 'conc_zn', vox)

# Insert no_voxel point process at soma(0.5) for each
voxels = []
for cell in cells:
    vox = h.no_voxel(cell.secs['soma']['hObj'](0.5))
    init_voxel(vox)
    voxels.append(vox)

###############################################################################
# 4. Set voxel parameters
###############################################################################

# define time points (ms)
tvec = h.Vector([0, 100, 200, 5000])     # simulation goes to 5000 ms
# corresponding values for F (nM/ms)
fvec = h.Vector([0,   5,   0,    0])

# play into voxel[0].F
fvec.play(voxels[0]._ref_F, tvec, 1)   # the last arg = interpolation flag

# Source in voxel 0
# voxels[0].F = 2     # production term [nM/ms]
voxels[0].lam = 0.01  # decay [1/ms]
# middle voxel
voxels[1].F = 0
voxels[1].lam = 0.01
# last voxel
voxels[2].F = 0
voxels[2].lam = 0.01

# symmetric diffusion coupling (1 <-> 2 <-> 3)
D = 0.05  # /ms
voxels[0].dx_pos = D
voxels[1].dx_neg = D; voxels[1].dx_pos = D
voxels[2].dx_neg = D

# set POINTER links

h.setpointer(voxels[1]._ref_conc, 'conc_xp', voxels[0])  # voxel0 → voxel1
h.setpointer(voxels[0]._ref_conc, 'conc_xn', voxels[1])  # voxel1 ← voxel0
h.setpointer(voxels[2]._ref_conc, 'conc_xp', voxels[1])  # voxel1 → voxel2
h.setpointer(voxels[1]._ref_conc, 'conc_xn', voxels[2])  # voxel2 ← voxel1

###############################################################################
# 5. Record concentrations
###############################################################################
t = h.Vector().record(h._ref_t)
c0 = h.Vector().record(voxels[0]._ref_conc)
c1 = h.Vector().record(voxels[1]._ref_conc)
c2 = h.Vector().record(voxels[2]._ref_conc)

###############################################################################
# 6. Run
###############################################################################
sim.runSim()

###############################################################################
# 7. Plot
###############################################################################

plt.plot(t, c0, label='Voxel 0')
plt.plot(t, c1, label='Voxel 1 (source)')
plt.plot(t, c2, label='Voxel 2')
plt.xlabel('Time (ms)')
plt.ylabel('NO conc (nM)')
plt.legend()
plt.savefig('test_conc_lines.png')


