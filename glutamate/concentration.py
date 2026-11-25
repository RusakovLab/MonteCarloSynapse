"""Convert shell molecule counts to glutamate concentration in uM.

Common to both receptor drivers: one molecule per um^3 corresponds to
``1e21 / N_A`` uM. The concentration handed to the kinetics is in uM.
"""
from __future__ import annotations

import numpy as np


def molecule_per_um3_to_uM(avogadro: float = 6.02214076e23) -> float:
    """uM produced by one molecule per um^3 (~1.6605e-3)."""
    return 1e21 / avogadro


def counts_to_uM(molecule_count, shell_volume_um3, avogadro: float = 6.02214076e23):
    """Concentration (uM) = factor * counts / shell_volume_um3.

    ``molecule_count`` is (n_shells, n_time); ``shell_volume_um3`` is (n_shells,).
    """
    molecule_count = np.asarray(molecule_count, dtype=float)
    shell_volume_um3 = np.asarray(shell_volume_um3, dtype=float).reshape(-1, 1)
    return molecule_per_um3_to_uM(avogadro) * molecule_count / shell_volume_um3
