from netpyne import specs, sim
from neuron import h
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import LogNorm
import numpy as np

lam_vals = [500, 1000, 2000, 5000]
diff_results = {lam : [] for lam in lam_vals}
for lam in lam_vals:
    ###############################################################################
    # Simulation configuration
    ###############################################################################
    cfg = specs.SimConfig()
    cfg.duration = 2000   # ms
    cfg.dt = 0.1
    cfg.verbose = False
    cfg.recordStep = 0.1


    ###############################################################################
    # NetPyNE network and cell definition (just dummy “cells” to host voxels)
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
        # 'numCells': pow(cube_side_len, 3),
        'gridSpacing': 11,
        'gridDim': [cube_side_len, cube_side_len, cube_side_len]
    }


    ###############################################################################
    # Create network
    ###############################################################################
    sim.initialize(simConfig=cfg, netParams=netParams)  # create network object and set cfg and net params
    sim.net.createPops()  # instantiate network populations
    sim.net.createCells()  # instantiate network cells based on defined populations


    cells = sim.net.cells

    def get_cell_coords(cell):
        tags = getattr(cell, 'tags', {}) or (cell.get('tags', {}) if isinstance(cell, dict) else {})
        x = tags['x']; y = tags['y']; z = tags['z']
        # snap to nearest grid node to kill tiny float noise
        q = lambda v: int(round(v / cube_side_len)) * cube_side_len
        return (q(x), q(y), q(z))

    # Insert no_voxel point process at soma(0.5) for each
    voxels = {}
    for cell in cells:
        vox = h.no_voxel(cell.secs['soma']['hObj'](0.5))
        voxels[get_cell_coords(cell)] = vox

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
    lam_actual = np.log(2)/lam # decay constant (/ms)

    for vox in voxels.values():
        vox.dx_pos = D
        vox.dx_neg = D
        vox.dy_pos = D
        vox.dy_neg = D
        vox.dz_pos = D
        vox.dz_neg = D
        vox.lam = lam_actual


    xs = sorted({x for x,_,_ in voxels})
    ys = sorted({y for _,y,_ in voxels})
    zs = sorted({z for _,_,z in voxels})

    center = (xs[len(xs)//2], ys[len(ys)//2], zs[len(zs)//2])
    right = (110, 55, 55)
    left = (00, 55, 55)
    top = (55, 77, 55)
    bottom = (55, 00, 55)

    assert center in voxels, f"center key {center} not in voxels!"
    # optional: verify it maps to the middle index
    idx = (xs.index(center[0]), ys.index(center[1]), zs.index(center[2]))
    assert idx == (len(xs)//2, len(ys)//2, len(zs)//2)

    # Example case to drop in a pulse of NO to this cube
    tvec = h.Vector([0, 500, 550, 650, 700, cfg.duration])
    fvec = h.Vector([0, 0, 2.5, 2.5, 0, 0])
    fvec.play(voxels[center]._ref_F, tvec, 1)

    # fvec_top = h.Vector([0, 0, 0, 0, 0, 0])          # values (nM/ms)
    # fvec_top.play(voxels[top]._ref_F, tvec, 1)

    ###############################################################################
    # Record concentrations
    ###############################################################################
    t = h.Vector().record(h._ref_t)
    sim.simData['t'] = t

    for key, vox in voxels.items():
        vec = h.Vector().record(vox._ref_conc)
        sim.simData[f"conc_{key}"] = vec

    ###############################################################################
    # 6. Run
    ###############################################################################
    sim.runSim()

    x_voxels = range(55, 121, 11)  # select voxels from source to edge in any direction

    # grab the actual concentration vectors for all of these locations
    x_vox_section = {}
    for vox in x_voxels:
        x_vox_section[(vox, 55, 55)] = sim.simData[f"conc_({vox}, 55, 55)"]

    for vox in x_vox_section:
        idx_to_max = np.argmax(x_vox_section[vox])
        diff_results[lam].append(((idx_to_max * cfg.dt)/1000))


###############################################################################
# 7. Plot
###############################################################################
plot_grid = False
plot_conc_trace = False
plot_diff_over_space = False
plot_max_no_by_dist = False

if plot_grid:
    def get_cell_info(cell):
        # Works with both Cell objects (cell.tags) and dicts (cell['tags'])
        tags = getattr(cell, 'tags', None)
        if tags is None and isinstance(cell, dict):
            tags = cell.get('tags', {})
        if not tags:
            return None
        if 'pos' in tags and isinstance(tags['pos'], (list, tuple)) and len(tags['pos']) == 3:
            x, y, z = tags['pos']
        else:
            x = tags.get('x', 0.0)
            y = tags.get('y', 0.0)
            z = tags.get('z', 0.0)
        pop = tags.get('pop', 'unknown')
        return x, y, z, pop
    # If running in parallel, you may have sim.net.allCells available with all ranks' cells
    cells = getattr(sim.net, 'allCells', None) or sim.net.cells
    pts = [get_cell_info(c) for c in cells]
    pts = [p for p in pts if p is not None]
    xs, ys, zs, pops = zip(*pts) if pts else ([], [], [], [])
    # Color by population
    uniq_pops = sorted(set(pops))
    pop_to_idx = {p: i for i, p in enumerate(uniq_pops)}
    colors = [pop_to_idx[p] for p in pops]
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection='3d')
    sc = ax.scatter(xs, ys, zs, c=colors, cmap='tab20', s=10, alpha=0.85)
    ax.set_xlabel('X (µm)')
    ax.set_ylabel('Y (µm)')
    ax.set_zlabel('Z (µm)')
    ax.set_title('Cell soma positions')
    # Make axes scale equally
    def set_axes_equal(ax):
        xr = np.array(ax.get_xlim3d())
        yr = np.array(ax.get_ylim3d())
        zr = np.array(ax.get_zlim3d())
        r = np.array([xr, yr, zr]); spans = r[:,1]-r[:,0]
        center = r.mean(axis=1)
        radius = spans.max()/2
        ax.set_xlim(center[0]-radius, center[0]+radius)
        ax.set_ylim(center[1]-radius, center[1]+radius)
        ax.set_zlim(center[2]-radius, center[2]+radius)
    set_axes_equal(ax)
    plt.tight_layout()
    plt.savefig('3Dnetfig.png')

if plot_conc_trace:

    t = np.array(sim.simData['t'])
    c_center = np.array(sim.simData[f'conc_{center}'])
    c_right = np.array(sim.simData[f'conc_{right}'])
    c_left = np.array(sim.simData[f'conc_{left}'])
    c_top = np.array(sim.simData[f'conc_{top}'])
    c_bottom = np.array(sim.simData[f'conc_{bottom}'])

    plt.figure()
    plt.plot(t, c_center, label='center')
    # plt.plot(t, c_right, label='right edge')
    # plt.plot(t, c_left, label='left edge')
    # plt.plot(t, c_top, label='top voxel')
    # plt.plot(t, c_bottom, label='bottom voxel')
    plt.xlabel('Time (ms)')
    plt.ylabel('NO concentration (nM)')
    plt.legend()
    plt.savefig('diff_cube_NOpulse_traces.png')
    plt.close()



if plot_diff_over_space:

    # get all voxel coordinates and sort them for consistency
    coords = list(voxels.keys())

    # determine grid size and shape
    xs = sorted(set([x for x, _, _ in coords]))
    ys = sorted(set([y for _, y, _ in coords]))
    zs = sorted(set([z for _, _, z in coords]))

    netParams.sizeX, netParams.sizeY, netParams.sizeZ = len(xs), len(ys), len(zs)

    # map coords to array indices
    coord_to_idx = {(x, y, z): (xs.index(x), ys.index(y), zs.index(z)) for (x, y, z) in coords}

    # Compute global min and max across all voxels and all timepoints
    all_conc = []
    for key in voxels:
        conc_vec = np.array(sim.simData[f"conc_{key}"])
        all_conc.append(conc_vec)
    all_conc = np.array(all_conc)

    timepoints = [500, 550, 1000, 1500, 2000]

    for time in timepoints:
        time_ms = time   # or 1000, etc.
        t_idx = (np.abs(t - time_ms)).argmin()  # nearest recorded index

        # fill 3D array with final concentration values
        conc_grid = np.zeros((netParams.sizeX, netParams.sizeY, netParams.sizeZ))
        for (x, y, z), vox in voxels.items():
            idx = coord_to_idx[(x, y, z)]
            conc_vec = np.array(sim.simData[f"conc_{(x, y, z)}"])
            conc_grid[idx] = conc_vec[t_idx]

        plt.figure()

        mid_z = len(zs)//2
        img = conc_grid[:, :, mid_z].T   # transpose so x→cols, y→rows
        vals = conc_grid.flatten()
        # vmin, vmax = np.percentile(vals, [0, 10])  # ignore outliers
        plt.imshow(img, origin='lower', cmap='plasma',
                   # vmin=0, vmax=0.2,
                   extent=[xs[0], xs[-1], ys[0], ys[-1]])

        plt.xlabel('X (µm)')
        plt.ylabel('Y (µm)')
        plt.title(f'NO concentration (z={zs[mid_z]} µm)')
        plt.colorbar(label='NO (nM)')
        plt.savefig(f'diff_cube_NOpulse_{time_ms}.png')
        plt.close()


if plot_max_no_by_dist:

    x_voxels = range(55, 121, 11) # select voxels from source to edge in any direction

    # grab the actual concentration vectors for all of these locations
    x_vox_section = {}
    for vox in x_voxels:
        x_vox_section[(vox, 55, 55)] = sim.simData[f"conc_({vox}, 55, 55)"]

    # Grab max value from source
    max = x_vox_section[(55, 55, 55)].max()

    # normalize values to the max
    norm_x_vox_section = {}
    for vox in x_vox_section:
        norm_x_vox_section[vox] = [x / max for x in x_vox_section[vox]]

    # create a vector of distances from the center voxel for plotting
    dist_vec = []
    for vox in norm_x_vox_section:
        dist_vec.append(vox[0] - 55)

    # create a vector of the max normalized concentration value in each voxel
    conc_vec = []
    for max in norm_x_vox_section:
        conc_vec.append(np.max(norm_x_vox_section[max]))

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(dist_vec, conc_vec)
    ax.set_xlabel('distance (µm)')
    ax.set_xticks(np.arange(0, 56, 5))
    ax.set_ylabel('[NOmax]/[NOmax Global]')
    ax.set_yticks(np.arrange(0, 1.1, 0.1))
    ax.set_title('NO concentration over distance')
    fig.savefig('NO_conc_by_dist.png')
