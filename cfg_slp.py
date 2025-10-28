from netpyne import specs

cfg = specs.SimConfig()
cfg.duration = 2000  # ms
cfg.dt = 0.05
cfg.verbose = False
cfg.cvode_active = False
cfg.cvode_atol = 1e-6
cfg.cache_efficient = True
cfg.printRunTime = 0.1

# Recording
cfg.recordCells = []
cfg.recordTraces = {}
cfg.recordStep = 0.1

cfg.savePickle = True

cfg.saveCellConns = False

cfg.analysis = {
    'plotRaster': {'saveFig': True},
    'plotConn':   {'feature': 'weight', 'saveFig': True},
}