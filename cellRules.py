# cellRules.py
# NetPyNE cell rules that mirror your cell_classes.py mechanisms + geometry,
# with NO sleep-state / init hooks. Only standard fixed params set; anything
# you want to modulate later (sleep, NO, etc.) can be layered on top.

import math

def build_cell_rules():
    rules = {}

    # ---------- Helper geometry ----------
    # Cortical soma area was ~100 µm^2 in your code:
    d_cx   = math.sqrt(100.0 / math.pi)          # ≈ 5.642 µm
    L_soma = d_cx
    # Dend lengths matched your code multipliers:
    L_pyr_dend = 165.0 * L_soma
    L_inh_dend =  50.0 * L_soma

    # Thalamic soma areas (from your code: S_RE=1.43e-4 cm^2, S_TC=2.9e-4 cm^2)
    # Convert cm^2 → µm^2 (×1e8), then to diameter for a sphere-equivalent cylinder
    d_re = math.sqrt((1.43e-4 * 1e8) / math.pi)  # ≈ 67.5 µm
    d_tc = math.sqrt((2.90e-4 * 1e8) / math.pi)  # ≈ 96.2 µm

    # Axial resistances (from your geom_factor math in cell_classes)
    # If you later want to tune, you can — these are fixed here.
    Ra_pyr = 1000.0 / (332.0 / (math.pi * d_cx))   # same derivation you had
    Ra_inh = 1000.0 / (102.0 / (math.pi * d_cx))
    Ra_th  = 100.0                                 # RE/TC used 100 in your code

    # =========================
    # PYR (2 sections: soma/dend)
    # =========================
    rules['Pyr'] = {
        'conds': {'cellType': 'Pyr'},
        'secs': {
            'soma': {
                'geom': {'L': L_soma, 'diam': d_cx, 'nseg': 1, 'Ra': Ra_pyr, 'cm': 0.75},
                'mechs': {
                    # Fast K, Na, and persistent Na (names as in your MODs)
                    'kdr': {'gkbar': 0.2},
                    'naf': {'gnabar': 3.0},
                    'nap': {'gnabar': 0.0003},
                },
                'ions': {'na': {'ena': 50}, 'k': {'ek': -90}},
            },
            'dend': {
                'geom': {'L': L_pyr_dend, 'diam': d_cx, 'nseg': 1, 'Ra': Ra_pyr, 'cm': 0.75},
                'mechs': {
                    # Dendritic Na (low), persistent Na, HVA Ca, M current, Ca dyn, K(Ca)
                    'naf': {'gnabar': 0.0008},
                    'nap': {'gnabar': 0.000042},
                    'hva': {'gcabar': 0.000012},
                    'km':  {'gkbar': 0.00002},
                    'cad': {'taur': 165.0, 'depth': 1.0, 'cainf': 2.4e-4},
                    'kca': {'gkcabar': 0.00005},
                    # K leak channel present in your class; leave at default (no sleep/init)
                    'kL':  {},           # no gkL set; uses MOD’s default (often 0)
                    # Passive as in your dend
                    'pas': {'g': 1.1e-05, 'e': -67.0},
                },
                'ions': {'na': {'ena': 50}, 'k': {'ek': -90}, 'ca': {'eca': 140}},
                'topol': {'parentSec': 'soma', 'parentX': 1.0, 'childX': 0.0},
            }
        }
    }

    # =========================
    # INH (2 sections: soma/dend)
    # =========================
    rules['Inh'] = {
        'conds': {'cellType': 'Inh'},
        'secs': {
            'soma': {
                'geom': {'L': L_soma, 'diam': d_cx, 'nseg': 1, 'Ra': Ra_inh, 'cm': 0.75},
                'mechs': {
                    'kdr': {'gkbar': 0.2},
                    'naf': {'gnabar': 2.5},
                },
                'ions': {'na': {'ena': 50}, 'k': {'ek': -90}},
            },
            'dend': {
                'geom': {'L': L_inh_dend, 'diam': d_cx, 'nseg': 1, 'Ra': Ra_inh, 'cm': 0.75},
                'mechs': {
                    'naf': {'gnabar': 0.0008},
                    'hva': {'gcabar': 0.000012},
                    'km':  {'gkbar': 0.000015},
                    'cad': {'taur': 165.0, 'depth': 1.0, 'cainf': 2.4e-4},
                    'kca': {'gkcabar': 0.00005},
                    'kL':  {},           # present, no gkL set
                    'pas': {'g': 9.0e-06, 'e': -70.0},
                },
                'ions': {'na': {'ena': 50}, 'k': {'ek': -90}, 'ca': {'eca': 140}},
                'topol': {'parentSec': 'soma', 'parentX': 1.0, 'childX': 0.0},
            }
        }
    }

    # =========================
    # RE (1 section: soma)
    # =========================
    rules['RE'] = {
        'conds': {'cellType': 'RE'},
        'secs': {
            'soma': {
                'geom': {'L': d_re, 'diam': d_re, 'nseg': 1, 'Ra': Ra_th, 'cm': 1.0},
                'mechs': {
                    'pas':   {'e': -77.0, 'g': 5.0e-05},
                    'kL':    {},            # present, no gkL set
                    'naf_re': {'gnabar': 0.100},
                    'kdr_re': {'gkbar': 0.010},
                    'it_re':  {'gcabar': 0.0022},
                    'cad':   {'depth': 1.0, 'taur': 5.0, 'cainf': 2.4e-4},
                },
                'ions': {'na': {'ena': 50}, 'k': {'ek': -95}, 'ca': {'cao': 2.0}},
            }
        }
    }

    # =========================
    # TC (1 section: soma)
    # =========================
    rules['TC'] = {
        'conds': {'cellType': 'TC'},
        'secs': {
            'soma': {
                'geom': {'L': d_tc, 'diam': d_tc, 'nseg': 1, 'Ra': Ra_th, 'cm': 1.0},
                'mechs': {
                    'pas':    {'e': -70.0, 'g': 1.0e-05},
                    'kL':     {},           # present, no gkL set
                    'naf_tc': {'gnabar': 0.090},
                    'kdr_tc': {'gkbar': 0.012},
                    'it_tc':  {'gcabar': 0.0025, 'qm': 3.55, 'qh': 3.0},
                    'iar':    {},           # H-current present; leave ghbar at MOD default
                    'cad':    {'depth': 2.0, 'taur': 5.0, 'cainf': 2.4e-4},
                },
                'ions': {'na': {'ena': 50}, 'k': {'ek': -95}, 'ca': {'cai': 2.4e-4, 'cao': 2.0, 'eca': 120}},
            }
        }
    }

    return rules
