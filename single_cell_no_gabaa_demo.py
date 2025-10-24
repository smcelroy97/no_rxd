from neuron import h
import numpy as np
import matplotlib.pyplot as plt

h.load_file('stdrun.hoc')

# Cell
soma = h.Section(name='soma')
soma.L = soma.diam = 20
soma.insert('pas'); soma.g_pas = 1e-4; soma.e_pas = -65

# Synapse (use your GABA kernel; amplitude fixed)
syn = h.MyExp2SynBB(soma(0.5))  # or Exp2Syn / MyExp2SynBB_NO with syn.alpha=0
syn.e = -80
syn.tau1 = 0.07
syn.tau2 = 18.2

# Injector
nc = h.NetCon(None, syn)
nc.weight[0] = 1.0    # amplitude per IPSC (constant)
nc.delay     = 0.1    # positive

# NO -> rate (Hz)
R0, RMAX, KNO = 1.0, 10.0, 1.0
def rate_from_NO(no):
    return max(0.0, R0 + RMAX * (no / (KNO + no)))
def NO_t(ms):
    return 0.0 if ms < 1000.0 else 2.5

# Run and schedule Poisson events
h.dt = 0.025
tstop  = 2000.0
win_ms = 10.0
t = 0.0

tvec = h.Vector().record(h._ref_t)
ivec = h.Vector().record(syn._ref_i)

rng = np.random.default_rng(1234)
h.finitialize(-65)

while t < tstop - 1e-9:
    tnext = min(t + win_ms, tstop)
    NO = NO_t(0.5*(t+tnext))
    lam = rate_from_NO(NO)                 # Hz
    n = rng.poisson(lam * (tnext - t) * 1e-3)
    if n > 0:
        U = rng.random(n)
        for te in (t + U*(tnext - t)):
            nc.event(float(te))
    h.continuerun(tnext)
    t = tnext

# Plot
import matplotlib.pyplot as plt
t = np.array(tvec); i = np.array(ivec)
plt.figure(figsize=(6,4))
plt.plot(t, i, lw=1)
plt.axvspan(0, 1000,   color='gray',   alpha=0.15, label='NO=0 nM')
plt.axvspan(1000, 2000, color='orange', alpha=0.15, label='NO=2.5 nM')
plt.xlabel('Time (ms)'); plt.ylabel('IPSC (nA)')
plt.title('NO increases IPSC frequency (amplitude fixed)')
plt.legend(); plt.tight_layout()
plt.savefig('demo_no_freq_single.png', dpi=150)
print('Saved demo_no_freq_single.png')

