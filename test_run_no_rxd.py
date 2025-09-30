from netpyne import specs, sim
from neuron import h, rxd
import numpy as np

"""
NetPyNE + NEURON RxD template for a 3-D cube (voxel) nitric-oxide model
that exactly follows the paper's linear multi-compartment diffusion:
    dC_i/dt = sum_neighbors D_i,neighbor (C_neighbor - C_i) - lambda_i C_i + F_i

Key knobs to sweep/optimize:
- DX_UM: cube edge (µm)
- D_PHYS_CM2_S: physical diffusion (cm^2/s)
- TORT: tortuosity (dimensionless) → D_eff = D_phys / TORT^2
- LAMBDA0_S: baseline decay (s^-1)

This scaffold:
- Builds a small population and places cells in 3D.
- Creates a 3D ECS lattice with dx=DX_UM.
- Creates NO species with d=0 (we supply operator via `rxd.Rate`).
- Computes D coefficients (s^-1) from D_eff/(dx_cm^2).
- Uses a 6-neighbor stencil with open boundaries.
- Adds spike-driven NO sources (alpha pulses) into somatic voxels.
- Feeds local NO back to cells (example: scales g_pas; swap to NMDA/HCN as needed).
"""

# =====================
# 0) Hyperparameters
# =====================
DX_UM        = 10.0           # cube size (µm)
D_PHYS_CM2_S = 1.5e-5         # free diffusion (cm^2/s)
TORT         = 1.7            # tortuosity
LAMBDA0_S    = 5.0            # baseline decay (s^-1)

VOL_X_UM, VOL_Y_UM, VOL_Z_UM = 120.0, 120.0, 80.0  # volume dimensions

# Spike source kernel (per spike into its soma voxel)
TAU_NO_MS = 15.0
A_NO      = 5e-5              # amplitude scaling (tune to land in ~10–50 nM peaks)

# Feedback strength (example: leak scaling)
ALPHA_LEAK = 5e-2

# =====================
# 1) NetPyNE config
# =====================
netParams = specs.NetParams()
simConfig = specs.SimConfig()

simConfig.duration = 300  # ms
simConfig.dt = 0.025
simConfig.hParams = {'celsius': 34, 'v_init': -65, 'cvode_active': True}
simConfig.verbose = False
simConfig.recordStep = 0.1

cellModels = ['HH_reduced', 'HH_full']

cellParamLabels = ['TC_reduced']

for ruleLabel in cellParamLabels:
    netParams.loadCellParamsRule(label=ruleLabel, fileName='cells/' + ruleLabel + '_cellParams.json')  # Load cellParams for each of the above cell subtype


# Population
N_CELLS = 16
netParams.popParams['TC'] = {'cellType': 'TC', 'cellModel': 'HH_reduced', 'size': N_CELLS, 'ynormRange': [1.2, 1.4]}

# Two IClamps to make a subset spike
netParams.stimSourceParams['icl_short'] = {'type': 'IClamp', 'del': 100, 'dur': 50,  'amp': 0.2}
netParams.stimSourceParams['icl_long']  = {'type': 'IClamp', 'del': 100, 'dur': 90,  'amp': 0.2}

# Build net structures (no run yet)
# sim.createSimulateAnalyze(netParams=netParams, simConfig=simConfig, output=False,
#                           create=True, simulate=False, analyze=False)

sim.initialize(simConfig=simConfig, netParams=netParams)  # create network object and set cfg and net params
sim.net.createPops()  # instantiate network populations
sim.net.createCells()  # instantiate network cells based on defined populations

# Place cells in 3D grid inside the volume
nx_cells = int(np.ceil(N_CELLS ** (1/3)))
coords = []
for i, cell in enumerate(sim.net.cells):
    ix = i % nx_cells
    iy = (i // nx_cells) % nx_cells
    iz = i // (nx_cells * nx_cells)
    # spread within volume
    x = (ix + 0.5) * (VOL_X_UM / nx_cells) - VOL_X_UM/2
    y = (iy + 0.5) * (VOL_Y_UM / nx_cells) - VOL_Y_UM/2
    z = (iz + 0.5) * (VOL_Z_UM / nx_cells) - VOL_Z_UM/2
    cell.tags.update({'x': x, 'y': y, 'z': z})
    sec = cell.secs['soma']['hObj']
    h.pt3dclear(sec=sec)
    h.pt3dadd(x, y, z, 20, sec=sec)
    coords.append((x, y, z))

# Attach IClamps to half the cells (alternating)
odd = [c.gid for c in sim.net.cells if (c.gid % 2) == 1]
even = [c.gid for c in sim.net.cells if (c.gid % 2) == 0]
netParams.stimTargetParams['short_to_even'] = {'source': 'icl_short', 'sec': 'soma', 'loc': 0.5, 'conds': {'cellList': even}}
netParams.stimTargetParams['long_to_odd']   = {'source': 'icl_long',  'sec': 'soma', 'loc': 0.5, 'conds': {'cellList': odd}}

# Re-realize new stim targets
sim.net.createCells(); sim.net.connectCells(); sim.net.addStims()

# =====================
# 2) RxD 3-D ECS + NO species (d=0)
# =====================
# ECS bounds centered on (0,0,0)
XLO = -VOL_X_UM/2; XHI = VOL_X_UM/2
YLO = -VOL_Y_UM/2; YHI = VOL_Y_UM/2
ZLO = -VOL_Z_UM/2; ZHI = VOL_Z_UM/2

ecs = rxd.Extracellular(XLO, YLO, ZLO, XHI, YHI, ZHI, dx=DX_UM)
no  = rxd.Species(ecs, name='no', charge=0, d=0.0)  # disable built-in diffusion
nodes = list(no.nodes)

# Grid shape
nx = int(ecs._nx) + 1
ny = int(ecs._ny) + 1
nz = int(ecs._nz) + 1
M  = nx*ny*nz

# Map each soma to its voxel index

def voxel_index_of(x, y, z):
    i = int(round((x - ecs._xlo)/ecs._dx[0]))
    j = int(round((y - ecs._ylo)/ecs._dx[1]))
    k = int(round((z - ecs._zlo)/ecs._dx[2]))
    return i + j*nx + k*nx*ny

voxel_index = [voxel_index_of(*p) for p in coords]

# =====================
# 3) Diffusion coefficients from (D_phys, dx, tort)
# =====================
DX_CM = DX_UM * 1e-4
D_eff = D_PHYS_CM2_S / (TORT**2)
Dcoef = D_eff / (DX_CM**2)  # s^-1

# 6-neighbor isotropic coefficients
Dx_pos = np.full(M, Dcoef); Dx_neg = np.full(M, Dcoef)
Dy_pos = np.full(M, Dcoef); Dy_neg = np.full(M, Dcoef)
Dz_pos = np.full(M, Dcoef); Dz_neg = np.full(M, Dcoef)

# Decay (can be spatially varying)
lam = np.full(M, LAMBDA0_S)

# Source vector (rebuilt each step)
F = np.zeros(M)

# Neighbor helper (open boundaries)

def nindex(ix, iy, iz, jx, jy, jz):
    if not (0 <= jx < nx and 0 <= jy < ny and 0 <= jz < nz):
        return None
    return jx + jy*nx + jz*nx*ny

# =====================
# 4) Operator application: work = H*C + F
# =====================
work  = np.zeros(M)
last_t = -1.0


def update_work():
    C = np.fromiter((nd.concentration for nd in nodes), dtype=float, count=M)
    w = - (Dx_pos + Dx_neg + Dy_pos + Dy_neg + Dz_pos + Dz_neg + lam) * C + F
    for iz in range(nz):
        for iy in range(ny):
            row = iy*nx + iz*nx*ny
            for ix in range(nx):
                i = ix + row
                # +x / -x
                j = nindex(ix, iy, iz, ix+1, iy, iz)
                if j is not None: w[i] += Dx_pos[i] * C[j]
                j = nindex(ix, iy, iz, ix-1, iy, iz)
                if j is not None: w[i] += Dx_neg[i] * C[j]
                # +y / -y
                j = nindex(ix, iy, iz, ix, iy+1, iz)
                if j is not None: w[i] += Dy_pos[i] * C[j]
                j = nindex(ix, iy, iz, ix, iy-1, iz)
                if j is not None: w[i] += Dy_neg[i] * C[j]
                # +z / -z
                j = nindex(ix, iy, iz, ix, iy, iz+1)
                if j is not None: w[i] += Dz_pos[i] * C[j]
                j = nindex(ix, iy, iz, ix, iy, iz-1)
                if j is not None: w[i] += Dz_neg[i] * C[j]
    work[:] = w

# One tiny Rate per node that reads from `work`

def rhs_factory(i):
    def rhs(node):
        global last_t
        if h.t != last_t:
            update_work()
            last_t = h.t
        return work[i]
    return rhs

rates = [rxd.Rate(no, rhs_factory(i), regions=[nd.region]) for i, nd in enumerate(nodes)]

# =====================
# 5) Spike → source pulses; NO → feedback
# =====================
# # Record spikes per cell
# spike_vecs = [h.Vector() for _ in sim.net.cells]
# for cell, v in zip(sim.net.cells, spike_vecs):
#     nc = h.NetCon(cell.secs['soma']['hObj'](0.5)._ref_v, None, sec=cell.secs['soma']['hObj'])
#     nc.threshold = 0.0
#     nc.record(v)
#
# active_spikes = []  # (voxel_idx, t_sp)
#
#
# def gather_new_spikes():
#     for i, vec in enumerate(spike_vecs):
#         while int(vec.size()) > 0:
#             t_sp = vec.x[0]
#             vec.remove(0)
#             active_spikes.append((voxel_index[i], t_sp))
#
#
# def rebuild_F():
#     t = h.t
#     F.fill(0.0)
#     keep = []
#     for vidx, ts in active_spikes:
#         a = (t - ts) / TAU_NO_MS
#         if 0 <= a < 10:
#             F[vidx] += A_NO * a * np.exp(1 - a)
#             keep.append((vidx, ts))
#     active_spikes[:] = keep
#
# # Example feedback: scale g_pas at soma by local NO (swap to NMDA/HCN as desired)
#
# def apply_feedback():
#     for i, cell in enumerate(sim.net.cells):
#         NOi = nodes[voxel_index[i]].concentration
#         cell.secs['soma']['hObj'].g_pas = 1e-4 * (1.0 + ALPHA_LEAK * NOi)
#
# def before_step():
#     gather_new_spikes()
#     rebuild_F()
#     apply_feedback()
#
# cv = h.CVode()
# h.CVode().extra_scatter_gather(before_step)

# =====================
# 6) Run
# =====================
sim.runSim()

# Optional: plot two example traces
try:
    from netpyne.analysis import plotTraces
    gids = [sim.net.cells[0].gid, sim.net.cells[-1].gid]
    plotTraces(include=gids, oneFig=True)
except Exception as e:
    print('Plotting skipped:', e)
