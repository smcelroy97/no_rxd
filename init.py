from no_utils import plotting
from netpyne import sim
from netParams import netParams, cfg
from neuron import h
import numpy as np
pc = h.ParallelContext()


sim.initialize(simConfig=cfg, netParams=netParams)  # create network object and set cfg and net params
sim.net.createPops()  # instantiate network populations
sim.net.createCells()  # instantiate network cells based on defined populations

def get_cell_coords(cell):
    tags = getattr(cell, 'tags', {}) or (cell.get('tags', {}) if isinstance(cell, dict) else {})
    x = tags['x']
    y = tags['y']
    z = tags['z']
    # snap to nearest grid node to kill tiny float noise
    q = lambda v: int(round(v / cfg.cube_side_len)) * cfg.cube_side_len
    return (q(x), q(y), q(z))


voxels = {}
for cell in sim.net.cells: 
    if cell.tags['pop'] == 'VoxelPop':                     # sim.net.cells are local on this rank
        vox = h.no_voxel(cell.secs['soma']['hObj'](0.5))
        voxels[get_cell_coords(cell)] = vox


###############################################################################
# Set voxel parameters
###############################################################################

offsets = {
    (11, 0, 0): 'conc_xp',
    (-11, 0, 0): 'conc_xn',
    (0, 11, 0): 'conc_yp',
    (0, -11, 0): 'conc_yn',
    (0, 0, 11): 'conc_zp',
    (0, 0, -11): 'conc_zn',
}

for vox in voxels:
    for (dx, dy, dz), pname in offsets.items():
        neighbor_key = (vox[0] + dx, vox[1] + dy, vox[2] + dz)
        if neighbor_key in voxels:
            h.setpointer(voxels[neighbor_key]._ref_conc, pname, voxels[vox])
        else:
            h.setpointer(voxels[vox]._ref_conc, pname, voxels[vox])  # boundary = self

lam = 1000  # Decay constant in ms

D_phys = 3.3          # µm²/ms
dx = 11.0             # µm
D = D_phys / (dx**2)  # 1/ms  -> 0.0273
lam_actual = np.log(2)/lam  # decay constant (/ms)

for vox in voxels.values():
    vox.dx_pos = D
    vox.dx_neg = D
    vox.dy_pos = D
    vox.dy_neg = D
    vox.dz_pos = D
    vox.dz_neg = D
    vox.lam = lam_actual

xs = sorted({x for x, _, _ in voxels})
ys = sorted({y for _, y, _ in voxels})
zs = sorted({z for _, _, z in voxels})

# explicit grid constants — keep these consistent with your NetPyNE sizes
GRID = cfg.cube_side_len  # 11.0 µm
CENTER_KEY = (55, 55, 55)  # known center for 110-µm box (or compute it explicitly)

# ... after you've created voxels from local cells ...
if CENTER_KEY in voxels:
    print(f"[host {int(pc.id())}] driving center voxel {CENTER_KEY}")
if CENTER_KEY in voxels:
    tvec = h.Vector([0, 420, 470, 570, 620, cfg.duration])
    fvec = h.Vector([0,   0, 250, 250,   0,           0])  # nM/ms

    # KEEP STRONG REFERENCES so GC won’t free them during the run
    if not hasattr(sim, '_Fplays'):
        sim._Fplays = []
    sim._Fplays.append((tvec, fvec, voxels[CENTER_KEY]))

    # fvec.play(voxels[CENTER_KEY]._ref_F, tvec, 1)

def snap(v):  # nearest voxel center
    return int(round(v / GRID)) * int(GRID)

###############################################################################
# Record and store [NO]
###############################################################################

# sim.simData['NO_conc'] = {}
# for key, vox in voxels.items():
#     vec = h.Vector().record(vox._ref_conc)
#     sim.simData['NO_conc'][f"conc_{key}"] = vec

###############################################################################
#  Run
###############################################################################

sim.net.connectCells()


# cell_to_voxel = {}
# for c in sim.net.cells:
#     if c.tags['pop'] == 'VoxelPop':
#         continue
#     else:
#         x, y, z = c.tags['x'], c.tags['y'], c.tags['z']
#         vkey = (snap(x), snap(y), snap(z))
#         if vkey not in voxels:
#             raise ValueError(f'No voxel at {vkey} for cell gid {c.gid}')
#         cell_to_voxel[c.gid] = vkey

# # 2) loop over postsynaptic cells and their conns; set POINTER on GABAA_NO synapses
# for c in sim.net.cells:
#     if c.tags['pop'] == 'VoxelPop':
#         continue
#     else:
#         vkey = cell_to_voxel[c.gid]  # NO field at postsynaptic site
#         for conn in c.conns:
#             # conn['synMech'] can be a string or list; normalize to list
#             mechs = conn['synMech'] if isinstance(conn['synMech'], list) else [conn['synMech']]
#             if 'GABAA_NO' in mechs:
#                 syn = conn['hSyn']  # hoc POINT_PROCESS MyExp2SynBB_NO
#                 # wire voxel -> synapse
#                 h.setpointer(voxels[vkey]._ref_conc, 'no_voxel', syn)


sim.net.addStims()
sim.setupRecording()  # setup variables to record for each cell (spikes, V traces, etc)


sim.runSim()
sim.gatherData()
sim.saveData()
sim.analysis.plotData()  # plot spike raster etc
