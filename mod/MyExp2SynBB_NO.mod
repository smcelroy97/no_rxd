TITLE Exp2Syn with NO modulation (drop-in for MyExp2SynBB)

NEURON {
    POINT_PROCESS MyExp2SynBB_NO
    RANGE tau1, tau2, e, i, g, gmax_base, alpha, K
    RANGE no_local
    NONSPECIFIC_CURRENT i
}

UNITS {
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (uS) = (microsiemens)
    (molar) = (1/liter)
    (nM) = (nanomolar)
}

PARAMETER {
    e          = -80 (mV)    : GABAA reversal, match your model
    tau1       = 0.07 (ms)   : rise time
    tau2       = 18.2 (ms)   : decay time
    gmax_base  = 0.001 (uS)  : peak conductance per weight=1 (baseline)
    alpha      = 2e-3 (/nM)  : NO sensitivity (linear/Hill mix below)
    K          = 100 (nM)    : half-saturation (set large to approximate linear)
}

ASSIGNED {
    v (mV)
    i (nA)
    g (uS)
    no_local (nM)             : comes from voxel via POINTER
    scale
}

STATE { A (uS)  B (uS) }     : Exp2 states

INITIAL {
}

BREAKPOINT {
    SOLVE states METHOD cnexp

    : --- NO modulation law ---
    : Saturating gain: scale = 1 + alpha * (no_local / (K + no_local))
    : If you want linear small-signal, set K very large (e.g., 1e9 nM)
    scale = 1 + alpha * (no_local / (K + no_local))

    g = (B - A) * scale
    i = g * (v - e)
}

DERIVATIVE states {
    A' = -A/tau1
    B' = -B/tau2
}

NET_RECEIVE (w) {
    : Keep same weight semantics as your MyExp2SynBB:
    : bump both states by w*gmax_base on each event
    A = A + w * gmax_base
    B = B + w * gmax_base
}