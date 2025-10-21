from netpyne import specs, sim
from neuron import h
from no_utils import plotting
import pickle
import numpy as np

lam_vals = [500, 1000, 2000, 5000]
diff_results = {lam: [] for lam in lam_vals}
for lam in lam_vals:

    ###############################################################################
    # Simulation configuration
    ###############################################################################

    cfg = specs.SimConfig()
    cfg.duration = 2000   # ms
    cfg.dt = 0.05   # Internal Integration Time Step
    cfg.verbose = False
    cfg.scaleDensity = 1.0  # Should be 1.0 unless need lower cell density for test simulation or visualization

    # cfg.recordTraces = {# 'V_soma': {'sec': 'soma', 'loc': 0.5, 'var': 'v'},
    #                     'NO_conc': {'sec': 'soma', 'loc': 0.5, 'var': 'conc'}
    #                     }

    cfg.recordStim = False  # Seen in M1 cfg.py
    cfg.recordTime = True  # SEen in M1 cfg.py
    cfg.recordStep = 0.05  # St ep size (in ms) to save data -- value from M1 cfg.py

    cfg.simLabel = '3d_no_demo'
    cfg.saveFolder = 'simOutput/' + cfg.simLabel  # Set file output name
    cfg.savePickle = False  # Save pkl file
    cfg.saveJson = True # Save json file
    cfg.saveDataInclude = ['simData', 'simConfig', 'netParams', 'net']
    cfg.backupCfgFile = None
    cfg.gatherOnlySimData = False
    cfg.saveCellSecs = False
    cfg.saveCellConns = False
    # cfg.recordCells = 'VoxelPop'

    ###############################################################################
    # NetPyNE network and cell definition
    ###############################################################################

    netParams = specs.NetParams()

    netParams.sizeX = 110
    netParams.sizeY = 110
    netParams.sizeZ = 110

    cube_side_len = 11

    # simple spherical cell with a soma section
    cellRule = {'conds': {'cellType': 'voxelhost', 'cellModel': 'HH'},
                'secs': {'soma': {'geom': {'L': 10, 'diam': 10}, 'mechs': {}}}}
    netParams.cellParams['VoxelHost'] = cellRule

    # population of 3 “cells” (each will host one voxel)
    netParams.popParams['VoxelPop'] = {
        'cellType': 'voxelhost',
        'cellModel': 'HH',
        'gridSpacing': 11,
        'gridDim': [cube_side_len, cube_side_len, cube_side_len]
    }


    ###############################################################################
    # Create network
    ###############################################################################

    sim.initialize(simConfig=cfg, netParams=netParams)  # create network object and set cfg and net params
    sim.net.createPops()  # instantiate network populations
    sim.net.createCells()  # instantiate network cells based on defined populations

    def get_cell_coords(cell):
        tags = getattr(cell, 'tags', {}) or (cell.get('tags', {}) if isinstance(cell, dict) else {})
        x = tags['x']
        y = tags['y']
        z = tags['z']
        # snap to nearest grid node to kill tiny float noise
        q = lambda v: int(round(v / cube_side_len)) * cube_side_len
        return (q(x), q(y), q(z))

    # Insert no_voxel point process at soma(0.5) for each
    voxels = {}
    for cell in sim.net.cells:
        vox = h.no_voxel(cell.secs['soma']['hObj'](0.5))
        voxels[get_cell_coords(cell)] = vox

    # for cell in sim.net.p:
    #     x, y, z = cell.tags['x'], cell.tags['y'], cell.tags['z']
    #
    #     # find nearest voxel center
    #     vx = round(x / 11) * 11
    #     vy = round(y / 11) * 11
    #     vz = round(z / 11) * 11
    #     voxel_key = (vx, vy, vz)
    #
    #     if voxel_key in voxels:
    #         gaba_mech = h.no_gaba(cell.secs['soma']['hObj'](0.5))
    #         h.setpointer(voxels[voxel_key]._ref_conc, 'conc_NO', gaba_mech)

    ###############################################################################
    # 4. Set voxel parameters
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

    D_phys = 3.3          # µm²/ms
    dx = 11.0             # µm
    D = D_phys / (dx**2)  # 1/ms  -> 0.0273
    lam_actual = np.log(2)/lam  # decay constant (/ms)
    # lam_actual = 1/lam
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

    center = (xs[len(xs)//2], ys[len(ys)//2], zs[len(zs)//2])
    right = (110, 55, 55)
    left = (00, 55, 55)
    top = (55, 77, 55)
    bottom = (55, 00, 55)

    voxels[center].conc0 = 240  # nM

    assert center in voxels, f"center key {center} not in voxels!"
    # optional: verify it maps to the middle index
    idx = (xs.index(center[0]), ys.index(center[1]), zs.index(center[2]))
    assert idx == (len(xs)//2, len(ys)//2, len(zs)//2)

    # Example case to drop in a pulse of NO to this cube
    tvec = h.Vector([0, 420, 470, 570, 620, cfg.duration])
    fvec = h.Vector([0, 0, 250, 250, 0, 0])
    fvec.play(voxels[center]._ref_F, tvec, 1)

    # fvec_top = h.Vector([0, 0, 0, 0, 0, 0])          # values (nM/ms)
    # fvec_top.play(voxels[top]._ref_F, tvec, 1)

    ###############################################################################
    # Record concentrations
    ###############################################################################
    # t = h.Vector().record(h._ref_t)
    # sim.simData['t'] = t

    for key, vox in voxels.items():
        vec = h.Vector().record(vox._ref_conc)
        sim.simData[f"conc_{key}"] = vec

    ###############################################################################
    # 6. Run
    ###############################################################################
    sim.runSim()
    # sim.gatherData()
    # sim.saveData()

    x_voxels = range(55, 121, 11)  # select voxels from source to edge in any direction

    # grab the actual concentration vectors for all of these locations
    x_vox_section = {}
    for vox in x_voxels:
        x_vox_section[(vox, 55, 55)] = sim.simData[f"conc_({vox}, 55, 55)"]

    for vox in x_vox_section:
        idx_to_max = np.argmax(x_vox_section[vox])
        diff_results[lam].append(((idx_to_max * cfg.dt)/1000))

    dist_vec = []
    for vox in x_vox_section:
        dist_vec.append(vox[0] - 55)

###############################################################################
# 7. Plot
###############################################################################

plot_time2max = True
plot_grid = False
plot_conc_heatmap = False
plot_max_conc_by_dist = False

if plot_time2max:
    plotting.time_to_max_conc(dist_vec, diff_results)

if plot_grid:
    plotting.voxel_net(sim)

if plot_conc_heatmap:
    plotting.conc_heat_map(sim, voxels)


if plot_max_conc_by_dist:
    plotting.max_conc_by_dist(sim, voxels)
