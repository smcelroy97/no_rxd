from neuron import h, crxd as rxd
import matplotlib.pyplot as plt
import numpy as np

h.load_file('stdrun.hoc')

Tfinal = 200.0
Dt = 0.1

NO_Diff = 3.3E3  # Diffusion constant in um^2/s
t_init = 50  # onset time (ms)

rhoF = 0.575  # in 1/ms
tF = 40  # duration of the pulse (ms)

Nseg = 23  # number of segments (odd number), select between [3,7,11,15,19,23, ...]
centerNode = int((Nseg - 1) / 2)
quarterNode = int(((Nseg - 1) / 2 - 1) / 2)


class Cylinder:

    def __init__(self, gid):
        self._gid = gid
        self._setup_morphology()
        self._setup_biophysics()

    def _setup_morphology(self):
        self.soma = h.Section(name='soma', cell=self)
        self.all = [self.soma]
        self.soma.L = 2.0
        self.soma.diam = 1
        self.soma.nseg = Nseg

    def _setup_biophysics(self):
        self.soma.Ra = 100  # Axial resistance in Ohm * cm
        self.soma.cm = 1  # Membrane capacitance in micro Farads / cm^2
        self.soma.insert('pas')
        for seg in self.soma:
            seg.pas.g = 0.001  # Passive conductance in S/cm2
            seg.pas.e = -65  # Leak reversal potential mV

    def __repr__(self):
        return 'Cylinder[{}]'.format(self._gid)


my_cell = Cylinder(0)

###############################################################################
## RxD ##
###############################################################################

## REGIONS
cyt = rxd.Region(my_cell.all, nrn_region='i', geometry=rxd.FractionalVolume(volume_fraction=1, surface_fraction=1))

## SPECIES & STATES
NO = rxd.Species([cyt], d=NO_Diff, name='NO', charge=0, initial=0, atolscale=1e-15)
Fvar = rxd.State([cyt], initial=0, atolscale=1e-15)

# PRODUCTION & DECAY
NO_prod = rxd.Rate(NO[cyt], rhoF * Fvar)

###############################################################################

## SIMULATION
h.finitialize(-65)
h.dt = Dt

times = []
cNO_cyt = []
cFvar_cyt = []

Flag1 = Flag2 = 1
Nt = round(Tfinal / Dt)

for nt in range(0, Nt, 10):

    t = nt * Dt
    if t > t_init and Flag1:
        Fvar[cyt].nodes[centerNode].concentration = 1
        h.CVode().re_init()
        Flag1 = 0
    if t > t_init + tF and Flag2:
        Fvar[cyt].nodes[centerNode].concentration = 0
        h.CVode().re_init()
        Flag2 = 0

    h.continuerun(t)
    print(t)
    times.append(t)
    cNO_cyt.append(NO[cyt].concentration)
    cFvar_cyt.append(Fvar[cyt].concentration)

plt.plot(times, [cNO_cyt[n][centerNode] for n in range(np.shape(cNO_cyt)[0])])
plt.plot(times, [cNO_cyt[n][quarterNode] for n in range(np.shape(cNO_cyt)[0])])