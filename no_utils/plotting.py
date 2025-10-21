from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np


def set_axes_equal(ax):
    xr = np.array(ax.get_xlim3d())
    yr = np.array(ax.get_ylim3d())
    zr = np.array(ax.get_zlim3d())
    r = np.array([xr, yr, zr])
    spans = r[:, 1] - r[:, 0]
    center = r.mean(axis=1)
    radius = spans.max() / 2
    ax.set_xlim(center[0] - radius, center[0] + radius)
    ax.set_ylim(center[1] - radius, center[1] + radius)
    ax.set_zlim(center[2] - radius, center[2] + radius)


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


def time_to_max_conc(
        dist_vec,  # A vector of distances from the source
        diff_results  # List of times in seconds that it took to reach max concentration at distances in dist_vec
):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    for trace in diff_results:
        ax.plot(dist_vec, diff_results[trace], label=f'{trace}')

    plt.legend()
    plt.ylabel('Time (s)')
    plt.xlabel('Distance (µm)')
    plt.show()
    fig.savefig('figs/time_to_max_conc.png')


def voxel_net(sim):

    # If running in parallel, you may have sim.net.allCells available with all ranks' cells
    cells = getattr(sim.net, 'allCells', None) or sim.net.cells
    pts = [get_cell_info(c) for c in cells if c.tags['pop'] == 'VoxelPop']
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
    set_axes_equal(ax)
    plt.tight_layout()
    plt.savefig('figs/3Dnetfig.png')


def conc_heat_map(sim, voxels):

    # get all voxel coordinates and sort them for consistency
    coords = list(voxels.keys())

    # determine grid size and shape
    xs = sorted(set([x for x, _, _ in coords]))
    ys = sorted(set([y for _, y, _ in coords]))
    zs = sorted(set([z for _, _, z in coords]))

    sim.net.params.sizeX, sim.net.params.sizeY, sim.net.params.sizeZ = len(xs), len(ys), len(zs)

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
        t_idx = (np.abs(sim.simData['t'] - time_ms)).argmin()  # nearest recorded index

        # fill 3D array with final concentration values
        conc_grid = np.zeros((sim.net.params.sizeX, sim.net.params.sizeY, sim.net.params.sizeZ))
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
        plt.savefig(f'figs/diff_cube_NOpulse_{time_ms}.png')
        plt.close()


def max_conc_by_dist(sim, voxels):

    x_voxels = range(55, 121, 11)  # select voxels from source to edge in any direction

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
    fig.savefig('figs/NO_conc_by_dist.png')
