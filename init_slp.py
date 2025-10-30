from netpyne import sim
from netParams_slp import netParams
from cfg_slp import cfg
from neuron import h

def seed_d2_syns():
    total = 0
    for cell in sim.net.cells:
        gid = int(getattr(cell, 'gid', 0))
        for conn in cell.conns:
            names = conn.get('synMech')
            synobjs = conn.get('synMechObj')

            # normalize to lists
            if names is None:
                continue
            if not isinstance(names, (list, tuple)):
                names = [names]
            if synobjs is None:
                # fallback via NetCon target (may return one obj only)
                synobjs = []
                hobj = conn.get('hObj', None)
                if hobj is not None and hasattr(hobj, 'syn'):
                    try:
                        synobjs = [hobj.syn()]
                    except Exception:
                        synobjs = []
            if not isinstance(synobjs, (list, tuple)):
                synobjs = [synobjs]

            # zip over mech names and objects (some templates create same-length lists)
            for name, sobj in zip(names, synobjs):
                if name in ('AMPA_D2', 'GABA_A_D2') and sobj is not None:
                    # published model uses a small integer syn_index category
                    if hasattr(sobj, 'syn_index'):
                        sobj.syn_index = 0 if name == 'AMPA_D2' else 1
                    # Random123 seeding
                    if hasattr(sobj, 'setrand'):
                        sobj.setrand(gid, int(getattr(sobj, 'syn_index', 0)))
                        total += 1
    print(f"Seeded {total} D2 synapses")

sim.initialize(netParams=netParams, simConfig=cfg)
print("[synMechParams names]", list(netParams.synMechParams.keys()))
sim.net.createPops()
sim.net.createCells()
sim.net.connectCells()
# seed_d2_syns()

# ===== Auto-discover and wire POINTER fields on syn targets =====
from neuron import h

def _syn_from_nc(nc):
    try:
        return nc.syn()
    except Exception:
        return None

def autodiscover_and_wire_ptr(presyn_loc=('soma',0.5), max_test_per_syn=200):
    """
    For each NetCon target PP, find which attribute accepts setpointer(..., attr, syn),
    remember per-class which attr works, and wire it to presyn V.
    Prints a summary mapping: class -> pointer_field and counts.
    """
    # per-class cache: which attr worked as POINTER
    class2field = {}
    class_counts = {}
    wired = skipped = missing_pre = 0
    tested_attrs = 0

    def _attrs(obj):
        try:
            return [a for a in dir(obj) if not a.startswith('__')]
        except Exception:
            return []

    # First pass: discover
    for post in sim.net.cells:
        for conn in post.conns:
            nc = conn.get('hObj', None)
            if nc is None:
                continue
            syn = _syn_from_nc(nc)
            if syn is None:
                continue
            cname = syn.hname() if hasattr(syn,'hname') else type(syn).__name__
            class_counts[cname] = class_counts.get(cname, 0) + 1
            if cname in class2field:
                continue  # already discovered for this class

            # Try to find any attr that accepts setpointer
            attrs = _attrs(syn)[:max_test_per_syn]
            # Use current post as a temporary source; we only need to see if attr is a POINTER
            # We'll wire properly in the second pass.
            for attr in attrs:
                try:
                    # Attempt a dry-run setpointer with a harmless pointer (we’ll override later)
                    # We need a SectionRef to call _ref_v; grab this cell's soma or any section
                    sec_item = None
                    for secname, secdict in post.secs.items():
                        sec_item = secdict['hObj'](0.5)
                        break
                    if sec_item is None:
                        continue
                    h.setpointer(sec_item._ref_v, attr, syn)
                    class2field[cname] = attr
                    break
                except Exception:
                    continue

    # Second pass: wire to *presynaptic* V using discovered field
    for post in sim.net.cells:
        for conn in post.conns:
            pre_gid = conn.get('preGid', None)
            nc = conn.get('hObj', None)
            if pre_gid is None or nc is None:
                continue
            syn = _syn_from_nc(nc)
            if syn is None:
                continue

            cname = syn.hname() if hasattr(syn,'hname') else type(syn).__name__
            attr = class2field.get(cname, None)
            if attr is None:
                continue  # no pointer field for this class

            # Find presyn cell and section
            try:
                pre = sim.net.cells[sim.net.gid2lid[pre_gid]]
            except Exception:
                missing_pre += 1
                continue
            pre_sec_name, pre_x = presyn_loc
            if pre_sec_name not in pre.secs:
                pre_sec_name, pre_x = ('soma', 0.5) if 'soma' in pre.secs else (next(iter(pre.secs.keys())), 0.5)
            pre_sec = pre.secs[pre_sec_name]['hObj']

            try:
                h.setpointer(pre_sec(pre_x)._ref_v, attr, syn)
                wired += 1
            except Exception:
                skipped += 1

    # Print summary
    print("[autowire] class -> pointer field:")
    for cname, cnt in sorted(class_counts.items(), key=lambda kv: -kv[1]):
        field = class2field.get(cname, None)
        print(f"   {cname:20s} : {field}   (n={cnt})")
    print(f"[autowire] wired={wired}, missing_pre={missing_pre}, skipped={skipped}")

# ---- call after connect & seeding ----
# autodiscover_and_wire_ptr()




sim.net.addStims()
sim.setupRecording()
sim.runSim()
sim.gatherData()
sim.saveData()
sim.analysis.plotData()
sim.close()
