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

netParams.scale = cfg.scale  # Scale factor for number of cells # NOT DEFINED YET! 3/11/19 # How is this different than scaleDensity?
netParams.sizeX = cfg.sizeX  # x-dimension (horizontal length) size in um
netParams.sizeY = cfg.sizeY  # y-dimension (vertical height or cortical depth) size in um
netParams.sizeZ = cfg.sizeZ  # z-dimension (horizontal depth) size in um
netParams.shape = 'cylinder'  # cylindrical (column-like) volume

# ------------------------------------------------------------------------------
# General connectivity parameters
# -----------------------------------------------------------------------------
netParams.scaleConnWeight = 1.0  # Connection weight scale factor (default if no model specified)
netParams.scaleConnWeightModels = {'HH_reduced': 1.0, 'HH_full': 1.0}  # scale conn weight factor for each cell model
netParams.scaleConnWeightNetStims = 1.0  # 0.5  # scale conn weight factor for NetStims
netParams.defaultThreshold = 0.0  # spike threshold, 10 mV is NetCon default, lower it for all cells
netParams.defaultDelay = 2.0  # default conn delay (ms)
netParams.propVelocity = 500.0  # propagation velocity (um/ms)
netParams.probLambda = 100.0  # length constant (lambda) for connection probability decay (um)
ThalamicCoreLambda = 50.0
# ------------------------------------------------------------------------------
# Cell parameters
# ------------------------------------------------------------------------------

cellModels = ['HH_reduced', 'HH_full']  # List of cell models

# II: 100-950, IV: 950-1250, V: 1250-1550, VI: 1550-2000
layer = {'thal': [1.2, 1.4],
         'cochlear': [1.6, 1.601]
         }  # normalized layer boundaries

cellParamLabels = ['RE_reduced', 'TC_reduced', 'HTC_reduced', 'TI_reduced']

for ruleLabel in cellParamLabels:
    netParams.loadCellParamsRule(label=ruleLabel, fileName='cells/' + ruleLabel + '_cellParams.json')  # Load cellParams for each of the above cell subtype


# ------------------------------------------------------------------------------
# Population parameters
# ------------------------------------------------------------------------------

## load densities
with open('cells/cellDensity.pkl', 'rb') as fileObj:
    density = pickle.load(fileObj)['density']
density = {k: [x * cfg.scaleDensity for x in v] for k, v in density.items()}  # Scale densities

### THALAMIC POPULATIONS (from prev model)
thalDensity = density[('A1', 'PV')][2] * 1.25  # temporary estimate (from prev model)

netParams.popParams['TC'] = {'cellType': 'TC', 'cellModel': 'HH_reduced', 'ynormRange': layer['thal'], 'density': 0.75 * thalDensity}
netParams.popParams['TCM'] = {'cellType': 'TC', 'cellModel': 'HH_reduced', 'ynormRange': layer['thal'], 'density': thalDensity}
netParams.popParams['HTC'] = {'cellType': 'HTC', 'cellModel': 'HH_reduced', 'ynormRange': layer['thal'], 'density': 0.25 * thalDensity}
netParams.popParams['IRE'] = {'cellType': 'RE', 'cellModel': 'HH_reduced', 'ynormRange': layer['thal'], 'density': thalDensity}
netParams.popParams['IREM'] = {'cellType': 'RE', 'cellModel': 'HH_reduced', 'ynormRange': layer['thal'], 'density': thalDensity}
netParams.popParams['TI'] = {'cellType': 'TI', 'cellModel': 'HH_reduced', 'ynormRange': layer['thal'], 'density': 0.33 * thalDensity}  # Winer & Larue 1996; Huang et al 1999
netParams.popParams['TIM'] = {'cellType': 'TI', 'cellModel': 'HH_reduced', 'ynormRange': layer['thal'], 'density': 0.33 * thalDensity}  # Winer & Larue 1996; Huang et al 1999

if cfg.singleCellPops:
    for pop in netParams.popParams.values():
        pop['numCells'] = 1

if cfg.reducedPop:
    for pop in netParams.popParams.values():
        pop['numCells'] = cfg.reducedPop

if hasattr(cfg, 'pops_active') and cfg.pops_active:
    pop_params_new = {}
    for pop in cfg.pops_active:
        if pop in netParams.popParams:
            pop_params_new[pop] = netParams.popParams[pop]
        else:
            print(f"Warning: pop '{pop}' not found in netParams.popParams")
    netParams.popParams = pop_params_new

# ------------------------------------------------------------------------------
# Synaptic mechanism parameters
# ------------------------------------------------------------------------------

### From M1 detailed netParams.py
netParams.synMechParams['NMDA'] = {'mod': 'MyExp2SynNMDABB', 'tau1NMDA': 15, 'tau2NMDA': 150, 'e': 0}
netParams.synMechParams['AMPA'] = {'mod': 'MyExp2SynBB', 'tau1': 0.05, 'tau2': 5.3, 'e': 0}
netParams.synMechParams['GABAASlow'] = {'mod': 'MyExp2SynBB', 'tau1': 2, 'tau2': 100, 'e': -80}
netParams.synMechParams['GABAB'] = {"mod": "MyExp2SynBB", "tau1": 41, "tau2": 642, "e": -105}


ESynMech = ['AMPA', 'NMDA']
ThalIESynMech = ['GABAASlow', 'GABAB']
ThalIISynMech = ['GABAASlow']

# ------------------------------------------------------------------------------
# Local connectivity parameters
# ------------------------------------------------------------------------------

## load data from conn pre-processing file
with open('conn/conn.pkl', 'rb') as fileObj:
    connData = pickle.load(fileObj)
pmat = connData['pmat']
lmat = connData['lmat']
wmat = connData['wmat']
bins = connData['bins']
connDataSource = connData['connDataSource']

wmat = cfg.wmat

# ------------------------------------------------------------------------------
# Thalamic connectivity parameters
# ------------------------------------------------------------------------------

TEpops = ['TC', 'TCM', 'HTC']
TIpops = ['IRE', 'IREM', 'TI', 'TIM']


def IsThalamicCore(ct):
    return ct == 'TC' or ct == 'HTC' or ct == 'IRE' or ct == 'TI'


def wireThal():
    # set intrathalamic connections
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
                if (pre == 'TC' and (post == 'TC' or post == 'HTC')) or (pre == 'HTC' and (post == 'TC' or post == 'HTC')):
                    prob = '%f * exp(-dist_x/%f)' % (pmat[pre][post], ThalamicCoreLambda)
                else:
                    prob = pmat[pre][post]
                netParams.connParams['ITh_' + pre + '_' + post] = {
                    'preConds': {'pop': pre},
                    'postConds': {'pop': post},
                    'synMech': syn,
                    'probability': prob,
                    'weight': wmat[pre][post] * gain,
                    'synMechWeightFactor': synWeightFactor,
                    'delay': 'defaultDelay+dist_3D/propVelocity',
                    'synsPerConn': 1,
                    'sec': 'soma'}


if cfg.addConn and cfg.addIntraThalamicConn:
    # print('adding intrathal')
    wireThal()


def prob2conv(prob, npre):
    # probability to convergence; prob is connection probability, npre is number of presynaptic neurons
    return int(0.5 + prob * npre)


# cochlea -> thal
def connectCochleaToThal():
    prob = '%f * exp(-dist_x/%f)' % (cfg.cochlearThalInput['probECore'], ThalamicCoreLambda)
    netParams.connParams['cochlea->ThalECore'] = {
        'preConds': {'pop': 'cochlea'},
        'postConds': {'pop': ['TC', 'HTC']},
        'sec': 'soma',
        'loc': 0.5,
        'synMech': ESynMech,
        'probability': prob,
        'weight': cfg.cochlearThalInput['weightECore'],
        'synMechWeightFactor': cfg.synWeightFractionEE,
        'delay': cfg.delayBkg}
    prob = '%f * exp(-dist_x/%f)' % (cfg.cochlearThalInput['probICore'], ThalamicCoreLambda)
    netParams.connParams['cochlea->ThalICore'] = {
        'preConds': {'pop': ['cochlea']},
        'postConds': {'pop': ['TI']},  # 'IRE',
        'sec': 'soma',
        'loc': 0.5,
        'synMech': ESynMech,
        'probability': prob,
        'weight': cfg.cochlearThalInput['weightICore'],
        'synMechWeightFactor': cfg.synWeightFractionEI,
        'delay': cfg.delayBkg}
    # cochlea -> Thal Matrix
    netParams.connParams['cochlea->ThalEMatrix'] = {
        'preConds': {'pop': ['cochlea']},
        'postConds': {'pop': ['TCM']},
        'sec': 'soma',
        'loc': 0.5,
        'synMech': ESynMech,
        'convergence': prob2conv(cfg.cochlearThalInput['probEMatrix'], numCochlearCells),
        'weight': cfg.cochlearThalInput['weightEMatrix'],
        'synMechWeightFactor': cfg.synWeightFractionEE,
        'delay': cfg.delayBkg}
    netParams.connParams['cochlea->ThalIMatrix'] = {
        'preConds': {'pop': ['cochlea']},
        'postConds': {'pop': ['TIM']},  # 'IREM',
        'sec': 'soma',
        'loc': 0.5,
        'synMech': ESynMech,
        'convergence': prob2conv(cfg.cochlearThalInput['probIMatrix'], numCochlearCells),
        'weight': cfg.cochlearThalInput['weightIMatrix'],
        'synMechWeightFactor': cfg.synWeightFractionEI,
        'delay': cfg.delayBkg}



# ------------------------------------------------------------------------------
# NetStim inputs (to simulate short external stimuli; not bkg)
# ------------------------------------------------------------------------------
if cfg.addNetStim:
    for key in [k for k in dir(cfg) if k.startswith('NetStim')]:
        params = getattr(cfg, key, None)
        [pop, ynorm, sec, loc, synMech, synMechWeightFactor, start, interval, noise, number, weight, delay] = \
            [params[s] for s in ['pop', 'ynorm', 'sec', 'loc', 'synMech', 'synMechWeightFactor', 'start', 'interval', 'noise', 'number', 'weight', 'delay']]
        netParams.stimSourceParams[key] = {'type': 'NetStim', 'start': start, 'interval': interval, 'noise': noise, 'number': number}

        if not isinstance(pop, list):
            pop = [pop]

        for eachPop in pop:
            # connect stim source to target

            netParams.stimTargetParams[key + '_' + eachPop] = {
                'source': key,
                'conds': {'pop': eachPop, 'ynorm': ynorm},
                'sec': sec,
                'loc': loc,
                'synMech': synMech,
                'weight': weight,
                'synMechWeightFactor': synMechWeightFactor,
                'delay': delay}
