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
    'cellType': 'Pyr', 'cellModel': 'HH', 'numCells': 50,
    'xRange': [CX_X0, CX_X1], 'yRange': [CX_Y0, CX_Y1], 'zRange': [CX_Z0, CX_Z1],
}
netParams.popParams['Inh'] = {
    'cellType': 'Inh', 'cellModel': 'HH', 'numCells': 10,
    'xRange': [CX_X0, CX_X1], 'yRange': [CX_Y0, CX_Y1], 'zRange': [CX_Z0, CX_Z1],
}

# Thalamus inside the 110 µm cube (aligned to your NO voxels)
netParams.popParams['RE'] = {
    'cellType': 'RE', 'cellModel': 'HH', 'numCells': 10,
    'xRange': [TH_X0, TH_X1], 'yRange': [TH_Y0, TH_Y1], 'zRange': [TH_Z0, TH_Z1],
}
netParams.popParams['TC'] = {
    'cellType': 'TC', 'cellModel': 'HH', 'numCells': 10,
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
R_ctx_tc = 120.0
R_ctx_re = 120.0

# --- Fast glutamate (AMPA): ~2 ms decay
netParams.synMechParams['AMPA'] = {
    'mod': 'Exp2Syn',
    'tau1': 0.2,      # ms  (rise)
    'tau2': 2.0,      # ms  (decay)
    'e': 0.0          # mV
}

# --- “NMDA-like” slow glutamate (no Mg block; just long decay)
# If you want closer-to-NMDA timecourse, lengthen tau2 (e.g., 50–100 ms).
netParams.synMechParams['NMDA'] = {
    'mod': 'Exp2Syn',
    'tau1': 2.0,      # ms
    'tau2': 80.0,     # ms
    'e': 0.0          # mV
}

# --- Fast GABA_A: 6–10 ms decay typical in cortex/thalamus
netParams.synMechParams['GABA_A'] = {
    'mod': 'Exp2Syn',
    'tau1': 0.3,      # ms
    'tau2': 6.0,      # ms
    'e': -70.0        # mV  (set to your E_Cl)
}

# --- Slow GABA_B: approximate with biphasic Exp2Syn (slow decay)
# If you want even slower, increase tau2 to 200–300+ ms.
netParams.synMechParams['GABA_B'] = {
    'mod': 'Exp2Syn',
    'tau1': 30.0,     # ms (slow rise)
    'tau2': 200.0,    # ms (slow decay)
    'e': -95.0        # mV (K+ dominated)
}

# Pyr -> Pyr (AMPA + NMDA in one rule)
netParams.connParams['Pyr->Pyr'] = {
    'preConds':  {'pop': 'Pyr'},
    'postConds': {'pop': 'Pyr'},
    'synMech':   ['AMPA', 'NMDA'],
    'weight':    [0.03, 0.0075],
    'delay':     0.1,
    'probability': f'(abs(dist_y) <= {same_slab_tol}) * (dist_3D <= {R_pyr_pyr})',
    'sec': 'dend',
    'loc': 0.5
}

# Pyr -> Inh
netParams.connParams['Pyr->Inh'] = {
    'preConds':  {'pop': 'Pyr'},
    'postConds': {'pop': 'Inh'},
    'synMech':   ['AMPA', 'NMDA'],
    'weight':    [0.75*0.12/5.0, 0.0],  # keep NMDA if you want; else 0
    'delay':     0.1,
    'sec': 'dend',
    'loc': 0.5,
    'probability': f'(dist_3D <= {R_pyr_inh})'
}

# Inh -> Pyr
netParams.connParams['Inh->Pyr'] = {
    'preConds':  {'pop': 'Inh'},
    'postConds': {'pop': 'Pyr'},
    'synMech':   ['GABA_A'],         # add GABA_B here if you want
    'weight':    [0.75*0.24],
    'delay':     0.1,
    'sec': 'dend',
    'loc': 0.5,
    'probability': f'(dist_3D <= {R_inh_pyr})'
}

# RE -> TC (GABAA + GABAB)
netParams.connParams['RE->TC'] = {
    'preConds':  {'pop': 'RE'},
    'postConds': {'pop': 'TC'},
    'synMech':   ['GABA_A', 'GABA_B'],
    'weight':    [0.05, 0.002],
    'delay':     0.1,
    'sec': 'soma',
    'loc': 0.5,
    'probability': f'(dist_3D <= {R_re_tc})'
}

# TC -> RE (AMPA)
netParams.connParams['TC->RE'] = {
    'preConds':  {'pop': 'TC'},
    'postConds': {'pop': 'RE'},
    'synMech':   ['AMPA'],
    'weight':    [0.5*0.05],
    'delay':     0.1,
    'sec': 'soma',
    'loc': 0.5,
    'probability': f'(dist_3D <= {R_tc_re})'
}

# TC -> Pyr (feedforward thalamo-cortical)
netParams.connParams['TC->Pyr'] = {
    'preConds':  {'pop': 'TC'},
    'postConds': {'pop': 'Pyr'},
    'synMech':   ['AMPA'],
    'weight':    [0.75*0.2/5.0],
    'delay':     0.1,
    'sec': 'dend',
    'loc': 0.5,
    'probability': 0.2
}

# --- Corticothalamic: Pyr -> TC (AMPA on TC soma) ---
netParams.connParams['Pyr->TC'] = {
    'preConds'   : {'pop': 'Pyr'},
    'postConds'  : {'pop': 'TC'},
    'synMech'    : ['AMPA'],           # same AMPA mech you already defined
    'weight'     : [0.20],             # <-- tune to match your g_AMPA_CX_TC scaling
    'delay'      : 0.5,                # ms
    'probability': 0.2, # <-- use your radius or variable if you have one
    'sec'        : 'soma',
    'loc'        : 0.5
}

# --- Corticothalamic: Pyr -> RE (AMPA on RE soma) ---
netParams.connParams['Pyr->RE'] = {
    'preConds'   : {'pop': 'Pyr'},
    'postConds'  : {'pop': 'RE'},
    'synMech'    : ['AMPA'],
    'weight'     : [0.15],             # <-- tune to match your g_AMPA_CX_RE scaling
    'delay'      : 0.5,
    'probability': 0.2, # <-- use your radius or variable if you have one
    'sec'        : 'soma',
    'loc'        : 0.5
}



# ----- Tonic DC clamps (match published model) -----

# Sources
netParams.stimSourceParams['tonicPyr'] = {'type': 'IClamp', 'del': 0.0, 'dur': 1e9, 'amp': -0.00068805}
netParams.stimSourceParams['tonicInh'] = {'type': 'IClamp', 'del': 0.0, 'dur': 1e9, 'amp': 0.005455}
netParams.stimTargetParams['tonicPyr->Pyr'] = {'source': 'tonicPyr', 'conds': {'pop': 'Pyr'}, 'sec': 'dend', 'loc': 0.5}
netParams.stimTargetParams['tonicInh->Inh'] = {'source': 'tonicInh', 'conds': {'pop': 'Inh'}, 'sec': 'dend', 'loc': 0.5}
