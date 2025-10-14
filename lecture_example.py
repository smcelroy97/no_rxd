from neuron import h, gui

# Create the neuron
soma = h.Section(name='soma')
soma.L = 12.615  # Length in microns
soma.diam = 12.615  # Diameter in microns
soma.Ra = 100  # Axial resistance

# Insert Hodgkin-Huxley mechanisms and passive properties
soma.insert('hh')  # Insert the Hodgkin-Huxley mechanism
# soma.insert('pas')  # Insert passive properties

# Create a current clamp
stim = h.IClamp(soma(0.5))  # Insert at the center of the soma
stim.delay = 5  # Delay before current onset (ms)
stim.dur = 0.1  # Duration of the current pulse (ms)
stim.amp = 0.9   # Amplitude of the current pulse (nA)

# Create a second stim for 2 APs
stim2 = h.IClamp(soma(0.5))
stim2.delay = 25
stim2.dur = 0.1
stim2.amp = 0.9

# Create a vector to record the membrane potential
v = h.Vector().record(soma(0.5)._ref_v)  # Record voltage
t = h.Vector().record(h._ref_t)  # Record time

h.load_file("stdrun.hoc")
h.finitialize(-65)

# Run the simulation
h.tstop = 40  # Simulation time (ms)
h.run()

# Plot results
from matplotlib import pyplot as plt
plt.plot(t, v)
plt.xlabel('Time (ms)')
plt.ylabel('Membrane potential (mV)')
plt.title('Hodgkin-Huxley Neuron Simulation')
plt.savefig('lectureex.png')
