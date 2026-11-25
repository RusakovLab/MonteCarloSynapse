"""Receptor kinetic models (ODE right-hand sides).

Verbatim Python ports of the verified MATLAB kinetics: AMPA.m, AMPA1.m,
AMPA2.m and NMDA_basic.m. Each function returns dy/dt for use with
scipy.integrate.solve_ivp. Glutamate concentration is supplied in uM.

All four schemes conserve total receptor probability (sum of states = const);
this is checked in tests/test_kinetics.py.
"""
from __future__ import annotations

import numpy as np


# --- Model metadata: state layout and how to read open / desensitized ----
MODELS = {
    "AMPA":  {"n": 6,  "y0": [1, 0, 0, 0, 0, 0],
              "open": lambda Y: Y[:, 5],
              "desens": lambda Y: Y[:, 3] + Y[:, 4]},
    "AMPA1": {"n": 6,  "y0": [1, 0, 0, 0, 0, 0],
              "open": lambda Y: Y[:, 4] + Y[:, 5],
              "desens": lambda Y: Y[:, 3]},
    "AMPA2": {"n": 12, "y0": [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              "open": lambda Y: Y[:, 5] + Y[:, 6] + Y[:, 7],
              "desens": lambda Y: Y[:, 8] + Y[:, 9] + Y[:, 10] + Y[:, 11]},
    "NMDA":  {"n": 5,  "y0": [1, 0, 0, 0, 0],
              "open": lambda Y: Y[:, 3],
              "desens": lambda Y: Y[:, 4]},
}


def ampa(t, y, glu):
    """6-state AMPA (Patneau & Mayer / Destexhe). y=[C0,C1,C2,D1,D2,O]."""
    g = float(glu(t))
    Rb, Ru1, Ru2, Rd, Rr, Ro, Rc = 13.0, 0.0059, 86.0, 0.9, 0.064, 2.7, 0.2
    N = 0.5
    rb = Rb * g
    d = np.empty(6)
    d[0] = N * (-rb * y[0] + Ru1 * y[1])
    d[1] = N * (rb * y[0] - (Ru1 + rb + Rd) * y[1] + Rr * y[3] + Ru2 * y[2])
    d[2] = N * (rb * y[1] - (Ru2 + Rd + Ro) * y[2] + Rr * y[4] + Rc * y[5])
    d[3] = N * (Rd * y[1] - Rr * y[3])
    d[4] = N * (Rd * y[2] - Rr * y[4])
    d[5] = N * (Ro * y[2] - Rc * y[5])
    return d


def ampa1(t, y, glu):
    """6-state AMPA (Raman & Trussell). y=[C0,C1,C2,D,O1,O2]."""
    g = float(glu(t))
    Rb, Ru1, Ru2, Rd, Rr = 13.0, 0.3, 200.0, 30.0, 0.02
    Ro1, Rc1, Ro2, Rc2 = 100.0, 2.0, 2.0, 0.25
    N = 0.5
    rb = Rb * g
    d = np.empty(6)
    d[0] = N * (-rb * y[0] + Ru1 * y[1])
    d[1] = N * (rb * y[0] - (Ru1 + rb) * y[1] + Ru2 * y[2])
    d[2] = N * (rb * y[1] - (Ru2 + Rd + Ro1 + Ro2) * y[2]
                + Rr * y[3] + Rc1 * y[4] + Rc2 * y[5])
    d[3] = N * (Rd * y[2] - Rr * y[3])
    d[4] = N * (Ro1 * y[2] - Rc1 * y[4])
    d[5] = N * (Ro2 * y[2] - Rc2 * y[5])
    return d


def ampa2(t, y, glu):
    """12-state GluR1 AMPA. y=[R0..R4,O2,O3,O4,D1..D4]."""
    g = float(glu(t))
    C = 1000.0
    alpha = 3100.0 / C
    beta = 8000.0 / C
    k_f = 2e7 / C**2
    k_r = 9000.0 / C
    gamma = 7.6 / C
    delta = 1800.0 / C
    N = 0.5
    rb = k_f * g
    d = np.empty(12)
    d[0] = N * (k_r * y[1] - 4 * rb * y[0])
    d[1] = N * (4 * rb * y[0] - (3 * rb + k_r) * y[1] + 2 * k_r * y[2]
                - delta * y[1] + gamma * y[8])
    d[2] = N * (3 * rb * y[1] - (2 * rb + 2 * k_r) * y[2] + 3 * k_r * y[3]
                - delta * y[2] + gamma * y[9] + alpha * y[5] - 2 * beta * y[2])
    d[3] = N * (2 * rb * y[2] - (rb + 3 * k_r) * y[3] + 4 * k_r * y[4]
                - delta * y[3] + gamma * y[10] + alpha * y[6] - 3 * beta * y[3])
    d[4] = N * (rb * y[3] - 4 * k_r * y[4] - delta * y[4]
                + gamma * y[11] + alpha * y[7] - 4 * beta * y[4])
    d[5] = N * (2 * beta * y[2] - alpha * y[5])
    d[6] = N * (3 * beta * y[3] - alpha * y[6])
    d[7] = N * (4 * beta * y[4] - alpha * y[7])
    d[8] = N * (delta * y[1] - gamma * y[8] - 3 * rb * y[8] + k_r * y[9])
    d[9] = N * (delta * y[2] - gamma * y[9] - k_r * y[9] - 2 * rb * y[9]
                + 3 * rb * y[8] + 2 * k_r * y[10])
    d[10] = N * (delta * y[3] - gamma * y[10] + 2 * rb * y[9] - 2 * k_r * y[10]
                 - rb * y[10] + 3 * k_r * y[11])
    d[11] = N * (delta * y[4] - gamma * y[11] + rb * y[10] - 3 * k_r * y[11])
    return d


def nmda(t, y, glu):
    """5-state NMDA (Destexhe x5). y=[unbound,single,double,open,desens]."""
    g = float(glu(t))
    Rb, Ru, Ro, Rc, Rd, Rr = 25e-3, 47.5e-3, 80e-3, 65e-3, 12.5e-3, 29.5e-3
    N = 0.5
    rb1 = Rb * 1e3 * g   # x1e3 -> working coefficient 25 /uM/ms (matches AMPA scale)
    d = np.empty(5)
    d[0] = N * (Ru * y[1] - rb1 * y[0])
    d[1] = N * (rb1 * y[0] - (Ru + rb1) * y[1] + Ru * y[2])
    d[2] = N * (rb1 * y[1] - (Ru + Rd + Ro) * y[2] + Rc * y[3] + Rr * y[4])
    d[3] = N * (Ro * y[2] - Rc * y[3])
    d[4] = N * (Rd * y[2] - Rr * y[4])
    return d


RHS = {"AMPA": ampa, "AMPA1": ampa1, "AMPA2": ampa2, "NMDA": nmda}
