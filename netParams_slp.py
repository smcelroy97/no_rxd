from netpyne import specs
from cellRules import build_cell_rules

netParams = specs.NetParams()
for label, rule in build_cell_rules().items():
    netParams.cellParams[label] = rule

# ---------------- Scene size (µm) ----------------
# Make room for cortex above thalamus cube
netParams.sizeX = 1200
netParams.sizeY = 1200
netParams.sizeZ = 1200

# ---------------- Thalamus cube placement ----------------
# Align this with your NO voxel lattice. If your voxel grid’s origin is (0,0,0),
# this puts the thalamus cube at x,y,z ∈ [0,110] µm. Change TH_X0/TH_Y0/TH_Z0 to shift it.
TH_SIDE = 110.0
TH_X0, TH_Y0, TH_Z0 = 0.0, 0.0, 0.0          # <-- adjust if your voxel cube starts elsewhere
TH_X1, TH_Y1, TH_Z1 = TH_X0+TH_SIDE, TH_Y0+TH_SIDE, TH_Z0+TH_SIDE

# ---------------- Cortex slab (example) ----------------
# Put cortex above the thalamus (choose any ranges you like)
CX_Y0, CX_Y1 = 400.0, 1000.0   # cortical band in Y
CX_X0, CX_X1 = 0.0, 1200.0
CX_Z0, CX_Z1 = 0.0, 1200.0

# ---------------- Populations (keep your counts) ----------------
netParams.popParams['Pyr'] = {
    'cellType': 'Pyr', 'cellModel': 'HH', 'numCells': 500,
    'xRange': [CX_X0, CX_X1], 'yRange': [CX_Y0, CX_Y1], 'zRange': [CX_Z0, CX_Z1],
}
netParams.popParams['Inh'] = {
    'cellType': 'Inh', 'cellModel': 'HH', 'numCells': 100,
    'xRange': [CX_X0, CX_X1], 'yRange': [CX_Y0, CX_Y1], 'zRange': [CX_Z0, CX_Z1],
}

# Thalamus inside the 110 µm cube (aligned to your NO voxels)
netParams.popParams['RE'] = {
    'cellType': 'RE', 'cellModel': 'HH', 'numCells': 100,
    'xRange': [TH_X0, TH_X1], 'yRange': [TH_Y0, TH_Y1], 'zRange': [TH_Z0, TH_Z1],
}
netParams.popParams['TC'] = {
    'cellType': 'TC', 'cellModel': 'HH', 'numCells': 100,
    'xRange': [TH_X0, TH_X1], 'yRange': [TH_Y0, TH_Y1], 'zRange': [TH_Z0, TH_Z1],
}

# ---------------- “Neighbor” connectivity in 3D ----------------
# Use local radii; keep your original weights/delays.
same_slab_tol = 25.0  # you can remove this if not needed

R_pyr_pyr = 150.0
R_pyr_inh = 80.0
R_inh_pyr = 150.0
R_re_tc   = 60.0
R_tc_re   = 60.0
R_tc_pyr  = 180.0

netParams.synMechParams['AMPA'] = {'mod': 'AMPA'}
netParams.synMechParams['GABA_A'] = {'mod': 'GABA_A'}
netParams.synMechParams['GABA_B'] = {'mod': 'GABA_B'}
netParams.synMechParams['AMPA_D2'] = {'mod': 'AMPA_D2', 'seed1': 123, 'seed2': 0}
netParams.synMechParams['NMDA_D1'] = {'mod': 'NMDA_D1'}
netParams.synMechParams['GABA_A_D2'] = {'mod': 'GABA_A_D2', 'seed1': 456, 'seed2': 0}

# Pyr -> Pyr (AMPA + NMDA in one rule)
netParams.connParams['Pyr->Pyr'] = {
    'preConds':  {'pop': 'Pyr'},
    'postConds': {'pop': 'Pyr'},
    'synMech':   ['AMPA_D2', 'NMDA_D1'],
    'weight':    [0.03,      0.0075],   # example values you had
    'delay':     0.1,
    'probability': f'(dist_3D <= {R_pyr_pyr})'
}

# Pyr -> Inh
netParams.connParams['Pyr->Inh'] = {
    'preConds':  {'pop': 'Pyr'},
    'postConds': {'pop': 'Inh'},
    'synMech':   ['AMPA_D2', 'NMDA_D1'],
    'weight':    [0.75*0.12/5.0, 0.0],  # keep NMDA if you want; else 0
    'delay':     0.1,
    'probability': f'(dist_3D <= {R_pyr_inh})'
}

# Inh -> Pyr
netParams.connParams['Inh->Pyr'] = {
    'preConds':  {'pop': 'Inh'},
    'postConds': {'pop': 'Pyr'},
    'synMech':   ['GABA_A_D2'],         # add GABA_B here if you want
    'weight':    [0.75*0.24],
    'delay':     0.1,
    'probability': f'(dist_3D <= {R_inh_pyr})'
}

# RE -> TC (GABAA + GABAB)
netParams.connParams['RE->TC'] = {
    'preConds':  {'pop': 'RE'},
    'postConds': {'pop': 'TC'},
    'synMech':   ['GABA_A', 'GABA_B'],
    'weight':    [0.05, 0.002],
    'delay':     0.1,
    'probability': f'(dist_3D <= {R_re_tc})'
}

# TC -> RE (AMPA)
netParams.connParams['TC->RE'] = {
    'preConds':  {'pop': 'TC'},
    'postConds': {'pop': 'RE'},
    'synMech':   ['AMPA'],
    'weight':    [0.5*0.05],
    'delay':     0.1,
    'probability': f'(dist_3D <= {R_tc_re})'
}

# TC -> Pyr (feedforward thalamo-cortical)
netParams.connParams['TC->Pyr'] = {
    'preConds':  {'pop': 'TC'},
    'postConds': {'pop': 'Pyr'},
    'synMech':   ['AMPA'],
    'weight':    [0.75*0.2/5.0],
    'delay':     0.1,
    'probability': f'(dist_3D <= {R_tc_pyr})'
}

# Thalamus↔Cortex projections (if you keep them): pick a modest radius to
# target the patch of cortex above the thalamic cube.
R_tc_pyr = 180.0
netParams.connParams['TC->Pyr_AMPA'] = {
    'preConds': {'pop': 'TC'}, 'postConds': {'pop': 'Pyr'},
    'synMech': 'AMPA', 'weight': 0.75*0.2/5.0, 'delay': 0.1,
    'probability': f'(dist_3D <= {R_tc_pyr})'
}
