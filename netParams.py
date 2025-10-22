from cfg import cfg
from netpyne.batchtools import specs
import pickle
import json

netParams = specs.NetParams()  # object of class NetParams to store the network parameters

# ------------------------------------------------------------------------------
# VERSION
# ------------------------------------------------------------------------------

netParams.version = 1

# ------------------------------------------------------------------------------
#
# NETWORK PARAMETERS
#
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# General network parameters
# ------------------------------------------------------------------------------

###############################################################################
# NetPyNE network and cell definition
###############################################################################

netParams = specs.NetParams()

# Load cellParams
cellParamLabels = ['RE_reduced', 'TC_reduced', 'HTC_reduced', 'TI_reduced']

for ruleLabel in cellParamLabels:
    netParams.loadCellParamsRule(label=ruleLabel, fileName='cells/' + ruleLabel + '_cellParams.json')  # Load cellParams for each of the above cell subtype

# ------------------------------------------------------------------------------
# General connectivity parameters
# ------------------------------------------------------------------------------
netParams.scaleConnWeight = 1.0  # Connection weight scale factor (default if no model specified)
netParams.scaleConnWeightModels = {'HH_reduced': 1.0, 'HH_full': 1.0}  # scale conn weight factor for each cell model
netParams.scaleConnWeightNetStims = 1.0  # 0.5  # scale conn weight factor for NetStims
netParams.defaultThreshold = 0.0  # spike threshold, 10 mV is NetCon default, lower it for all cells
netParams.defaultDelay = 2.0  # default conn delay (ms)
netParams.propVelocity = 500.0  # propagation velocity (um/ms)
netParams.probLambda = 100.0  # length constant (lambda) for connection probability decay (um)
ThalamicCoreLambda = 50.0
netParams.sizeX = cfg.sizeX  # x-dimension (horizontal length) size in um
netParams.sizeY = cfg.sizeY  # y-dimension (vertical height or cortical depth) size in um
netParams.sizeZ = cfg.sizeZ  # z-dimension (horizontal depth) size in um

# ------------------------------------------------------------------------------
# Population parameters
# ------------------------------------------------------------------------------

with open('cells/cellDensity.pkl', 'rb') as fileObj:
    density = pickle.load(fileObj)['density']
density = {k: [x * cfg.scaleDensity for x in v] for k, v in density.items()}  # Scale densities

# simple spherical cell with a soma section
cellRule = {'conds': {'cellType': 'voxelhost', 'cellModel': 'HH'},
            'secs': {'soma': {'geom': {'L': 10, 'diam': 10}, 'mechs': {}}}}
netParams.cellParams['VoxelHost'] = cellRule

netParams.popParams['VoxelPop'] = {
    'cellType': 'voxelhost',
    'cellModel': 'HH',
    'gridSpacing': 11,
    'gridDim': [cfg.cube_side_len, cfg.cube_side_len, cfg.cube_side_len]
}

### THALAMIC POPULATIONS (from prev model)
thalDensity = density[('A1', 'PV')][2] * 1.25  # temporary estimate (from prev model)

netParams.popParams['TC'] = {'cellType': 'TC', 'cellModel': 'HH_reduced', 'yRange': [55, 110], 'density': 0.75 * thalDensity}
netParams.popParams['TCM'] = {'cellType': 'TC', 'cellModel': 'HH_reduced', 'yRange': [55, 110], 'density': thalDensity}
netParams.popParams['HTC'] = {'cellType': 'HTC', 'cellModel': 'HH_reduced', 'yRange': [55, 110], 'density': 0.25 * thalDensity}
netParams.popParams['IRE'] = {'cellType': 'RE', 'cellModel': 'HH_reduced', 'yRange': [0, 50], 'density': thalDensity}
netParams.popParams['IREM'] = {'cellType': 'RE', 'cellModel': 'HH_reduced', 'yRange': [0, 50], 'density': thalDensity}
netParams.popParams['TI'] = {'cellType': 'TI', 'cellModel': 'HH_reduced', 'yRange': [55, 110], 'density': 0.33 * thalDensity}  # Winer & Larue 1996; Huang et al 1999
netParams.popParams['TIM'] = {'cellType': 'TI', 'cellModel': 'HH_reduced', 'yRange': [55, 110], 'density': 0.33 * thalDensity}  # Winer & Larue 1996; Huang et al 1999

# ------------------------------------------------------------------------------
# Synaptic mechanism parameters
# ------------------------------------------------------------------------------

netParams.synMechParams['NMDA'] = {
    'mod': 'MyExp2SynNMDABB',
    'tau1NMDA': 15,
    'tau2NMDA': 150,
    'e': 0
}

netParams.synMechParams['AMPA'] = {
    'mod': 'MyExp2SynBB',
    'tau1': 0.05,
    'tau2': 5.3,
    'e': 0
}

netParams.synMechParams['GABAA'] = {
    'mod': 'MyExp2SynBB',
    'tau1': 0.07,
    'tau2': 18.2,
    'e': -80
}

netParams.synMechParams['GABAB'] = {
    "mod": "MyExp2SynBB",
    "tau1": 41,
    "tau2": 642,
    "e": -105
}

netParams.synMechParams['GABAASlow'] = {
    'mod': 'MyExp2SynBB',
    'tau1': 2,
    'tau2': 100,
    'e': -80
}

netParams.synMechParams['GABAA_NO'] = {
    'mod': 'MyExp2SynBB_NO',
    'tau1': 0.07,
    'tau2': 18.2,
    'e': -80,
    'gmax_base': 1e-3,  # or whatever your MyExp2SynBB expects for weight=1
    'alpha': 2e-3,  # tune
    'K': 100  # nM, tune
}

ESynMech = ['AMPA', 'NMDA']
ThalIESynMech = ['GABAA_NO', 'GABAB']
ThalIISynMech = ['GABAASlow']

with open('conn/conn.pkl', 'rb') as fileObj:
    connData = pickle.load(fileObj)
pmat = connData['pmat']
lmat = connData['lmat']
wmat = connData['wmat']
bins = connData['bins']
connDataSource = connData['connDataSource']

# ------------------------------------------------------------------------------
# Thalamic connectivity parameters
# ------------------------------------------------------------------------------

TEpops = ['TC', 'TCM', 'HTC']
TIpops = ['IRE', 'IREM', 'TI', 'TIM']

for pre in TEpops + TIpops:
    for post in TEpops + TIpops:
        gain = cfg.intraThalamicGain
        if post in pmat[pre]:
            # for syns use ESynMech, ThalIESynMech and ThalIISynMech
            if pre in TEpops:  # E->E/I
                syn = ESynMech
                synWeightFactor = cfg.synWeightFractionEE
                if post in TEpops:
                    gain *= cfg.intraThalamicEEGain
                else:
                    gain = cfg.intraThalamicEIGain
            elif post in TEpops:  # I->E
                syn = ThalIESynMech
                synWeightFactor = cfg.synWeightFractionThalIE
                gain *= cfg.intraThalamicIEGain
            else:  # I->I
                syn = ThalIISynMech
                synWeightFactor = cfg.synWeightFractionThalII
                gain *= cfg.intraThalamicIIGain
            # use spatially dependent wiring between thalamic core excitatory neurons
            # if (pre == 'TC' and (post == 'TC' or post == 'HTC')) or (pre == 'HTC' and (post == 'TC' or post == 'HTC')):
            #     prob = '%f * exp(-dist_x/%f)' % (pmat[pre][post], ThalamicCoreLambda)
            # else:
            prob = pmat[pre][post]
            netParams.connParams['ITh_' + pre + '_' + post] = {
                'preConds': {'pop': pre},
                'postConds': {'pop': post},
                'synMech': syn,
                'probability': 1,
                'weight': wmat[pre][post] * gain,
                'synMechWeightFactor': synWeightFactor,
                'delay': 'defaultDelay+dist_3D/propVelocity',
                'synsPerConn': 1,
                'sec': 'soma'}

# ------------------------------------------------------------------------------
# Background inputs
# ------------------------------------------------------------------------------
if cfg.addBkgConn:
    # add bkg sources for E and I cells
    netParams.stimSourceParams['excBkg'] = {
        'type': 'NetStim',
        'start': cfg.startBkg,
        'rate': cfg.rateBkg['exc'],
        'noise': cfg.noiseBkg,
        'number': 1e9}
    netParams.stimSourceParams['inhBkg'] = {
        'type': 'NetStim',
        'start': cfg.startBkg,
        'rate': cfg.rateBkg['inh'],
        'noise': cfg.noiseBkg,
        'number': 1e9}

    # excBkg/I -> thalamus + cortex
    with open('cells/bkgWeightPops.json', 'r') as f:
        weightBkg = json.load(f)
    pops = list(cfg.allpops)

    for pop in ['TC', 'TCM', 'HTC']:
        weightBkg[pop] *= cfg.EbkgThalamicGain

    for pop in ['IRE', 'IREM', 'TI', 'TIM']:
        weightBkg[pop] *= cfg.IbkgThalamicGain

    for pop in pops:
        netParams.stimTargetParams['excBkg->'+pop] = {
            'source': 'excBkg',
            'conds': {'pop': pop},
            'sec': 'soma',
            'loc': 0.5,
            'synMech': ESynMech,
            'weight': weightBkg[pop],
            'synMechWeightFactor': cfg.synWeightFractionEE,
            'delay': cfg.delayBkg}

        netParams.stimTargetParams['inhBkg->'+pop] = {
            'source': 'inhBkg',
            'conds': {'pop': pop},
            'sec': 'soma',
            'loc': 0.5,
            'synMech': 'GABAA',
            'weight': weightBkg[pop],
            'delay': cfg.delayBkg}
