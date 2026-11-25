"""Numerical checks for the ported models (mirror the README verification)."""
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from glutamate import kinetics, default_config
from glutamate.concentration import molecule_per_um3_to_uM


def _pulse(peak=1.0, t0=0.5, t1=1.5):
    return lambda t: peak if (t0 <= t <= t1) else 0.0


def test_probability_conservation():
    for name, rhs in kinetics.RHS.items():
        meta = kinetics.MODELS[name]
        y0 = np.array(meta["y0"], float)
        sol = solve_ivp(rhs, (0, 200), y0, method="LSODA",
                        args=(_pulse(),), rtol=1e-9, atol=1e-12, max_step=0.5)
        sums = sol.y.sum(axis=0)
        assert np.max(np.abs(sums - 1.0)) < 1e-6, name


def test_sum_of_derivatives_zero():
    rng = np.random.default_rng(0)
    for name, rhs in kinetics.RHS.items():
        n = kinetics.MODELS[name]["n"]
        worst = max(abs(rhs(0.7, rng.random(n), _pulse()).sum()) for _ in range(20))
        assert worst < 1e-10, name


def test_relaxation_to_unbound():
    for name, rhs in kinetics.RHS.items():
        n = kinetics.MODELS[name]["n"]
        y0 = np.full(n, 1.0 / n)
        T = 200000 if name == "AMPA" else 5000  # AMPA relaxes slowly
        sol = solve_ivp(rhs, (0, T), y0, method="LSODA",
                        args=(lambda t: 0.0,), rtol=1e-9, atol=1e-12)
        assert sol.y[0, -1] > 0.999, (name, sol.y[0, -1])


def test_concentration_constant():
    assert abs(molecule_per_um3_to_uM() - 1.6605e-3) < 1e-6


def test_config_links():
    cfg = default_config()
    assert cfg.ampa.TimeGlu_ms == cfg.basic.TotalTime
    assert cfg.ampa.BinGlu == cfg.computation.DistanceBins


if __name__ == "__main__":
    test_probability_conservation()
    test_sum_of_derivatives_zero()
    test_relaxation_to_unbound()
    test_concentration_constant()
    test_config_links()
    print("All kinetics tests passed.")
