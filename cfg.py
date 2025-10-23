from netpyne import specs
import pickle
import json
import numpy as np


cfg = specs.SimConfig()

# Insert params from previous tuning
with open('data/initCfg.json', 'rb') as f:
    cfgLoad = json.load(f)

for key, value in cfgLoad.items():
    setattr(cfg, key, value)


# ------------------------------------------------------------------------------
# Run parameters
# ------------------------------------------------------------------------------

cfg.duration = 1000  # Duration of the sim, in ms
cfg.dt = 0.05   # 0.025  # Internal Integration Time Step
cfg.verbose = 0  # Show detailed messages
cfg.progressBar = 0  # even more detailed message
cfg.hParams['celsius'] = 37
cfg.createNEURONObj = 1
cfg.createPyStruct = 1
cfg.printRunTime = 0.1

cfg.connRandomSecFromList = False  # set to false for reproducibility
cfg.cvode_active = False
cfg.cvode_atol = 1e-6
cfg.cache_efficient = True
cfg.printRunTime = 0.1  			# specified above
cfg.oneSynPerNetcon = False
cfg.includeParamsLabel = False
cfg.validateNetParams = True
cfg.use_coreNEURON = False

# ------------------------------------------------------------------------------
# Recording
# ------------------------------------------------------------------------------

cfg.TEpops = ['TC', 'TCM', 'HTC']
cfg.TIpops = ['IRE', 'IREM', 'TI', 'TIM']

cfg.allpops = cfg.TEpops + cfg.TIpops

# cfg.recordCells = ['all']

cfg.recordTraces = {'V_soma': {'sec': 'soma', 'loc': 0.5, 'var': 'v'}
                    # 'g_GABAA_NO': {'sec': 'soma', 'loc': 0.5, 'synMech': 'GABAA_NO', 'var': 'g'}
                    }

# ------------------------------------------------------------------------------
# Saving
# ------------------------------------------------------------------------------

cfg.simLabel = 'no_rxd_thal_v1'
cfg.saveFolder = 'simOutput/' + cfg.simLabel  # Set file output name
cfg.savePickle = True  # Save pkl file
cfg.saveJson = False  # Save json file
cfg.saveDataInclude = ['simData', 'simConfig', 'netParams', 'net']
cfg.backupCfgFile = None
cfg.gatherOnlySimData = False
cfg.saveCellSecs = False
cfg.saveCellConns = True

# ------------------------------------------------------------------------------
# Network
# ------------------------------------------------------------------------------

# These values taken from M1 cfg (https://github.com/Neurosim-lab/netpyne/bflob/development/examples/M1detailed/cfg.py)
cfg.singleCellPops = False
cfg.reducedPop = False  # insert number to declare specific number of populations, if going for full model set to False
cfg.singlePop = ''
cfg.removeWeightNorm = False
cfg.cube_side_len = 11
cfg.sizeX = cfg.cube_side_len * 10
cfg.sizeY = cfg.cube_side_len * 10
cfg.sizeZ = cfg.cube_side_len * 10
cfg.scaleDensity = 1.0  # Should be 1.0 unless need lower cell density for test simulation or visualization
cfg.connRandomSecFromList = False  # set to false for reproducibility
cfg.cvode_active = False
cfg.addBkgConn = 1.0
# ------------------------------------------------------------------------------
# Analysis and plotting
# ------------------------------------------------------------------------------

# cfg.analysis['plotTraces'] = {
#     'include': cfg.allpops,
#     'timeRange': [0, cfg.duration],
#     'oneFigPer': 'trace',
#     'overlay': True,
#     'saveFig': True,
#     'showFig': False,
#     'figSize': (12, 8)
# }

# ------------------------------------------------------------------------------
# Synapses
# ------------------------------------------------------------------------------

# General Synaptic Parameters
cfg.synWeightFractionEE = [0.5, 0.5]  # E->E AMPA to NMDA ratio
cfg.synWeightFractionEI = [0.5, 0.5]  # E->I AMPA to NMDA ratio
cfg.synWeightFractionIE = [0.9, 0.1]
cfg.synWeightFractionII = [1.0]

# Thalamic Synaptic Parameters
cfg.synWeightFractionThalIE = [0.9, 0.2]
cfg.synWeightFractionThalII = [1.0, 0.0]

cfg.intraThalamicGain = 1.0
cfg.intraThalamicEEGain = 1.0
cfg.intraThalamicEIGain = 0.3
cfg.intraThalamicIEGain = 0.1
cfg.intraThalamicIIGain = 1.0

# ------------------------------------------------------------------------------
# Connectivity
# ------------------------------------------------------------------------------

# full weight conn matrix
with open('conn/conn.pkl', 'rb') as fileObj:
    connData = pickle.load(fileObj)
cfg.wmat = connData['wmat']


# ------------------------------------------------------------------------------
# Background inputs
# ------------------------------------------------------------------------------
cfg.addBkgConn = 1.0
cfg.noiseBkg = 1.0  # firing rate random noise
cfg.delayBkg = 5.0  # (ms)
cfg.startBkg = 0  # start at 0 ms
cfg.rateBkg = {'exc': 40, 'inh': 40}

cfg.EbkgThalamicGain = 1.0  # 0.392
cfg.IbkgThalamicGain = 1.0  # 1.96
