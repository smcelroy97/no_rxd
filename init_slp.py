from netpyne import sim
from netParams_slp import netParams
from cfg_slp import cfg
from neuron import h

def seed_d2_syns():
    count = 0
    for cell in sim.net.cells:
        gid = int(getattr(cell, 'gid', 0))
        per_cell_index = 0
        for syn in cell.synMechs:
            mod = syn.get('mod') or syn.get('label')
            if mod in ('AMPA_D2', 'GABA_A_D2'):
                hsyn = syn['hObj']
                if hasattr(hsyn, 'setrand'):
                    # Give each syn its own deterministic stream using (gid, per_cell_index)
                    hsyn.setrand(gid, per_cell_index)
                    per_cell_index += 1
                    count += 1
    print(f"Seeded {count} D2 synapses")

sim.initialize(netParams=netParams, simConfig=cfg)
sim.net.createPops()
sim.net.createCells()
sim.net.connectCells()
sim.net.addStims()
sim.setupRecording()
sim.runSim()
sim.gatherData()
sim.saveData()
sim.analysis.plotData()
sim.close()
