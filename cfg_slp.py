from netpyne import specs

cfg = specs.SimConfig()

cfg.saveDataInclude = ['simData', 'simConfig', 'netParams', 'net']
cfg.duration = 2000  # ms
cfg.dt = 0.025
cfg.verbose = False
cfg.cvode_active = True
cfg.validateNetParams = True
cfg.cvode_atol = 1e-6
cfg.cache_efficient = True
cfg.printRunTime = 0.1
cfg.recordTime = True

# Recording
cfg.recordCells = [0]
cfg.recordTraces = {'V_soma': {'sec': 'soma', 'loc': 0.5, 'var': 'v'}}
cfg.recordStep = 0.1

cfg.savePickle = True

cfg.saveCellConns = True

cfg.analysis = {
    'plotRaster': {'saveFig': True},
    'plotConn':   {'feature': 'weight', 'saveFig': True},
}

cfg.analysis['plotTraces'] = {
    'include': [0],
    'timeRange': [0, cfg.duration],
    'oneFigPer': 'trace',
    'overlay': True,
    'saveFig': True,
    'showFig': False,
    'figSize': (12, 8)
}
