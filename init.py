from netpyne import sim
from neuron import h
from netParams import netParams, cfg
import numpy as np

pc = h.ParallelContext()
rank = int(pc.id())
nhost = int(pc.nhost())

# ----------------------------
# 1) Build the regular network
# ----------------------------
sim.initialize(simConfig=cfg, netParams=netParams)
sim.net.createPops()
sim.net.createCells()

# --- gid -> pop map (used to read prePop from conns) ---
gid2pop = {}
for pop_name, pop in sim.net.pops.items():
    gids = getattr(pop, 'cellGids', []) or getattr(pop, 'gids', []) or []
    for gid in gids:
        gid2pop[gid] = pop_name

def pre_pop_of_conn(conn):
    pg = conn.get('preGid', None)
    if pg is not None and pg in gid2pop:
        return gid2pop[pg]
    return conn.get('prePop', None)


sim.net.connectCells()
sim.net.addStims()
sim.setupRecording()   # spikes, V traces, etc. (doesn't include NO; we'll handle NO below)

# ---------- helpers ----------
GRID = cfg.cube_side_len      # e.g. 11.0 µm
BOX  = (cfg.sizeX, cfg.sizeY, cfg.sizeZ)  # e.g. (110,110,110)

def snap(v): return int(round(v/GRID))*int(GRID)

# ----------------------------------------------------------
# 2) Rank 0: build the entire NO lattice and link neighbors
# ----------------------------------------------------------
voxels_rank0 = {}   # {(x,y,z): hoc no_voxel}
offsets = {
    ( GRID, 0, 0): 'conc_xp', (-GRID, 0, 0): 'conc_xn',
    (0,  GRID, 0): 'conc_yp', (0, -GRID, 0): 'conc_yn',
    (0, 0,  GRID): 'conc_zp', (0, 0, -GRID): 'conc_zn',
}

# --- who to target ---
PRE_OK   = {'IRE', 'IREM'}
POST_OK  = {'TC', 'TCM', 'HTC'}
GABA_MECHS = {'RETCSynMech'}   # edit to your actual names

# --- lattice geometry (same as rank-0 lattice) ---
GRID = cfg.cube_side_len                  # e.g., 11.0
xs = list(range(0, cfg.sizeX + 1, int(GRID)))
ys = list(range(0, cfg.sizeY + 1, int(GRID)))
zs = list(range(0, cfg.sizeZ + 1, int(GRID)))
nx, ny, nz = len(xs), len(ys), len(zs)

def clamp_idx(i, n): return max(0, min(n-1, i))
def pos_to_idx(x, y, z):
    ix = clamp_idx(int(round(x/GRID)), nx-1)
    iy = clamp_idx(int(round(y/GRID)), ny-1)
    iz = clamp_idx(int(round(z/GRID)), nz-1)
    return ix, iy, iz
def idx3(ix, iy, iz): return (iz*ny + iy)*nx + ix

# --- collect local postsynaptic synapses to be frequency-driven ---
freq_targets = []   # (syn_hoc, ix, iy, iz)

for cell in sim.net.cells:  # local postsynaptic cells on this rank
    post_pop = cell.tags.get('pop')
    if post_pop not in POST_OK:
        continue
    x, y, z = cell.tags['x'], cell.tags['y'], cell.tags['z']
    ix, iy, iz = pos_to_idx(x, y, z)   # anchor at soma; OK for first pass

    for conn in cell.conns:
        mechs = conn['synMech'] if isinstance(conn['synMech'], list) else [conn['synMech']]
        if not any(m in GABA_MECHS for m in mechs):
            continue
        if pre_pop_of_conn(conn) not in PRE_OK:
            continue
        syn = conn.get('hSyn', None)
        if syn is None:
            continue
        # If your GABAA mech has NO amplitude modulation, neutralize it:
        if hasattr(syn, 'alpha'):  # MyExp2SynBB_NO style
            syn.alpha = 0.0
        freq_targets.append((syn, ix, iy, iz))

print(f"[rank {rank}] NO-freq targets on this rank: {len(freq_targets)}")

# --- NO -> rate parameters (tune) ---
R0   = 1.0    # Hz baseline minis
RMAX = 10.0   # Hz additional at high NO
KNO  = 100.0  # nM half-saturation
W    = 1.0    # NetCon weight (keeps amplitude constant)

def rate_from_NO(NO_nM):
    return max(0.0, R0 + RMAX * (NO_nM / (KNO + NO_nM)))

freq_drivers = []  # (nc, rng, ix, iy, iz)
seed_base = 54321 + rank*100000

for k, (syn, ix, iy, iz) in enumerate(freq_targets):
    nc = h.NetCon(None, syn)                     # programmatic event source
    nc.weight[0] = W
    nc.delay     = max(0.1, float(sim.cfg.dt))   # >0 keeps parallel mindelay happy
    rng = np.random.default_rng(seed_base + k)   # deterministic per rank
    freq_drivers.append((nc, rng, ix, iy, iz))

if rank == 0:
    # create a dummy Section host for all voxels so they live on rank 0
    host_sec = h.Section(name='no_host_rank0')

    # instantiate all voxels on the host section at 0.5
    xs = list(range(0, BOX[0]+1, int(GRID)))
    ys = list(range(0, BOX[1]+1, int(GRID)))
    zs = list(range(0, BOX[2]+1, int(GRID)))
    for x in xs:
        for y in ys:
            for z in zs:
                pp = h.no_voxel(host_sec(0.5))
                voxels_rank0[(x,y,z)] = pp

    # link neighbors (all local on rank 0 → safe)
    for (x,y,z), v in voxels_rank0.items():
        for (dx,dy,dz), pname in offsets.items():
            nk = (x+dx, y+dy, z+dz)
            target = voxels_rank0.get(nk, v)  # self at boundary
            h.setpointer(target._ref_conc, pname, v)

    # set global physical params on all voxels
    D_phys = 3.3              # µm^2/ms
    lam_ms = cfg.no_t_half_ms # or any source you use; lam = ln(2)/t_half
    lam = np.log(2)/lam_ms
    D = D_phys / (GRID**2)    # 1/ms
    for v in voxels_rank0.values():
        v.dx_pos = v.dx_neg = v.dy_pos = v.dy_neg = v.dz_pos = v.dz_neg = D
        v.lam = lam

    # optional initial condition / F schedule (center voxel)
    cx, cy, cz = xs[len(xs)//2], ys[len(ys)//2], zs[len(zs)//2]
    center_key = (cx, cy, cz)
    voxels_rank0[center_key].conc0 = 240.0  # nM, if you want that initial condition

    # Example: play any F schedule you want (only on rank 0)
    tvec = h.Vector([0, 420, 470, 570, 620, cfg.duration])
    fvec = h.Vector([0,   0, 250, 250,   0,           0])  # nM/ms
    fvec.play(voxels_rank0[center_key]._ref_F, tvec, 1)

pc.barrier()

# ----------------------------------------------------------------
# 3) Precompute synapse → lattice mapping (local synapses per rank)
#    (we’ll set syn.no_local from the broadcasted grid each sync)
# ----------------------------------------------------------------
# Gather all local synapses that use the NO-aware mechanism
syn_list = []  # list of (syn_hoc, ix, iy, iz, weights[8]) for trilinear OR (syn_hoc, ix, iy, iz) for nearest
use_trilinear = False  # start with nearest-neighbor (fast & simple)

# Build lattice coordinate arrays (same as rank 0):
xs = list(range(0, BOX[0]+1, int(GRID)))
ys = list(range(0, BOX[1]+1, int(GRID)))
zs = list(range(0, BOX[2]+1, int(GRID)))
nx, ny, nz = len(xs), len(ys), len(zs)

def clamp_idx(i, n): return max(0, min(n-1, i))

for cell in sim.net.cells:  # local postsynaptic cells on this rank
    x, y, z = cell.tags['x'], cell.tags['y'], cell.tags['z']
    # snap (or compute trilinear) indices
    ix = clamp_idx(int(round(x/GRID)), nx-1)
    iy = clamp_idx(int(round(y/GRID)), ny-1)
    iz = clamp_idx(int(round(z/GRID)), nz-1)

    for conn in cell.conns:
        # Normalize synMech to list and check for our NO-enabled syn
        mechs = conn['synMech'] if isinstance(conn['synMech'], list) else [conn['synMech']]
        if 'GABAA_NO' in mechs:
            syn = conn['hSyn']  # POINT_PROCESS MyExp2SynBB_NO with RANGE no_local
            if use_trilinear:
                # compute 8-corner weights here (omitted for brevity)
                # syn_list.append((syn, ix0,iy0,iz0, w8tuple))
                pass
            else:
                syn_list.append((syn, ix, iy, iz))

pc.barrier()

# ---------------------------------------------------------
# 4) Run in chunks; broadcast the grid; update syn.no_local
# ---------------------------------------------------------
tstop = cfg.duration
sync_dt = cfg.recordStep  # how often to sync NO field (ms)
t = 0.0

# Helper to pack/unpack the grid
def pack_grid_rank0():
    # Return a flat float64 array of shape (nx*ny*nz,)
    arr = np.empty(nx*ny*nz, dtype=np.float64)
    idx = 0
    for iz, z in enumerate(zs):
        for iy, y in enumerate(ys):
            for ix, x in enumerate(xs):
                arr[idx] = voxels_rank0[(x,y,z)].conc
                idx += 1
    return arr

def update_syns_from_grid(flat):
    # flat is 1D np.array length nx*ny*nz
    # nearest-neighbor assignment:
    for syn, ix, iy, iz in syn_list:
        idx = (iz*ny + iy)*nx + ix
        syn.no_local = float(flat[idx])

# Main loop
while t < tstop - 1e-9:
    tnext = min(t + sync_dt, tstop)

    # advance simulation to tnext
    pc.psolve(tnext)   # NEURON parallel solve to tnext
    t = tnext

    # rank 0 packs grid and broadcasts
    if rank == 0:
        grid_flat = pack_grid_rank0()
    else:
        grid_flat = None

    # Broadcast to all ranks (NEURON supports Python object broadcast in modern versions)
    try:
        grid_flat = pc.py_broadcast(grid_flat, 0)  # broadcast from rank 0
    except:
        # Fallback: use bytes (older NEURON). You can implement your own broadcast via pc.broadcast if needed.
        raise RuntimeError("Your NEURON lacks pc.py_broadcast; use bytes-based broadcast here.")

    window_ms = tnext - t
    window_s = window_ms * 1e-3

    events_this_step = 0
    for (nc, rng, ix, iy, iz) in freq_drivers:
        NO = float(grid_flat[idx3(ix, iy, iz)])
        lam_hz = rate_from_NO(NO)
        n = rng.poisson(lam_hz * window_s)
        if n == 0: continue
        events_this_step += int(n)
        if n == 1:
            nc.event(float(t + rng.random() * window_ms))
        else:
            U = rng.random(n)
            for te in (t + U * window_ms):
                nc.event(float(te))

    # all ranks update their local synapses’ no_local from the received grid
    update_syns_from_grid(grid_flat)

pc.barrier()

# -----------------------------------
# 5) Finish and save/plot with NetPyNE
# -----------------------------------
sim.gatherData()
sim.saveData()
sim.analysis.plotData()  # optional
