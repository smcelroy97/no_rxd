
TITLE Nitric Oxide Diffusion

COMMENT Equations from
        Pablo Fernandez Lopez, Patricio Garcia Baez, and Carmen Paz Suarez Araujo. Nitric Oxide Diﬀusion and Multi-compartmental
        Systems: Modeling and Implications DOI: 10.1007/978-3-319-26555-1 59

        Instituto Universitario de Ciencias y Tecnologias Cibernticas,
        Universidad de Las Palmas de Gran Canaria

        Adapted to be used in NetPyNE by Scott McElroy
ENDCOMMENT

NEURON {
    SUFFIX no_voxel
    RANGE conc, dx_pos, dx_neg, dy_pos, dy_neg, dz_pos, dz_neg, lam, F
    POINTER conc_xp, conc_xn, conc_yp, conc_yn, conc_zp, conc_zn
}

UNITS {
    (nM) = (nanomolar)
    (s) = (second)
}

PARAMETER {
    dx_pos = 0 (1/s)
    dx_neg = 0 (1/s)
    dy_pos = 0 (1/s)
    dy_neg = 0 (1/s)
    dz_pos = 0 (1/s)
    dz_neg = 0 (1/s)
    lam    = 0 (1/s)
    F      = 0 (nM/s)
}

ASSIGNED {
    conc_xp (nM)
    conc_xn (nM)
    conc_yp (nM)
    conc_yn (nM)
    conc_zp (nM)
    conc_zn (nM)
}

STATE {
    conc (nM)
}

INITIAL {
    conc = 0
}

DERIVATIVE diffusion {
    conc' = - (dx_pos + dx_neg + dy_pos + dy_neg + dz_pos + dz_neg + lam) * conc + F
          + (dx_pos * conc_xp) + (dx_neg * conc_xn)
          + (dy_pos * conc_yp) + (dy_neg * conc_yn)
          + (dz_pos * conc_zp) + (dz_neg * conc_zn)
}
