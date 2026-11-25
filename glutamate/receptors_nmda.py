"""Spatial NMDA analysis driver (port of NMDA_SpaceSR_Optimized.m).

Reconstructs shell concentration in uM from radial free-glutamate histograms
and runs the 5-state NMDA model per shell. Real shell volumes / coordinates can
be supplied directly (e.g. from a diffusion run); otherwise the run's
SimulationMetadata (JSON written here, or MATLAB .mat) is used when present,
and finally the legacy geometry.
"""
from __future__ import annotations

import glob
import os
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from . import io_files, kinetics
from .concentration import counts_to_uM


@dataclass
class NmdaResult:
    time_ms: np.ndarray
    distance_um: np.ndarray
    open_state: np.ndarray          # (n_shells, n_out)
    desensitization: np.ndarray     # (n_shells, n_out)
    concentration_uM: np.ndarray    # (n_shells, n_time)
    glu_time_ms: np.ndarray
    model_type: str = "NMDA"


def run_nmda_spatial(cfg, data_folder=None, initial_distribution=None,
                     shell_volume_um3=None, glu_time_ms=None, distance_um=None):
    if data_folder is None and initial_distribution is None:
        raise ValueError(
            "run_nmda_spatial needs either initial_distribution=<matrix> or "
            "data_folder=<folder with DistanceFree*.txt> (see run.py nmda --data)")

    meta = None
    if initial_distribution is None:
        files = sorted(glob.glob(os.path.join(data_folder, "DistanceFree*.txt")),
                       key=os.path.getmtime, reverse=True)
        if not files:
            raise FileNotFoundError("No DistanceFree*.txt in %s" % data_folder)
        chosen = files[0]
        initial_distribution = np.loadtxt(chosen)
        meta = io_files.load_metadata(data_folder, chosen)

    distance_bins = initial_distribution.shape[0]
    glu_bins = initial_distribution.shape[1] - 1

    if shell_volume_um3 is None:
        if meta is not None:
            shell_volume_um3 = meta["ShellVolume_um3"]
            if glu_time_ms is None:
                glu_time_ms = meta["TimePoints_ms"]
            if distance_um is None:
                distance_um = meta["OuterRadius_um"]
        else:
            shell_volume_um3 = _legacy_shell_volumes(distance_bins)
    if glu_time_ms is None:
        glu_time_ms = np.arange(1, glu_bins + 1) * (cfg.nmda.LegacyMasterTotalTime_ms / glu_bins)
    if distance_um is None:
        distance_um = np.linspace(cfg.cleft.RadiusofCleft / 1000.0,
                                  cfg.nmda.LegacyMaximumDistance_um, distance_bins)

    nmda_total = cfg.nmda.TotalTime_ms
    counts = initial_distribution[:, 1:]
    conc_uM = counts_to_uM(counts, shell_volume_um3, cfg.analysis.AvogadroNumber)

    out_time = np.linspace(0.0, nmda_total, glu_bins)
    y0 = np.array(kinetics.MODELS["NMDA"]["y0"], dtype=float)
    open_read = kinetics.MODELS["NMDA"]["open"]
    desens_read = kinetics.MODELS["NMDA"]["desens"]
    selected = np.zeros((distance_bins, glu_bins))
    desens = np.zeros((distance_bins, glu_bins))

    for shell in range(distance_bins):
        trace = conc_uM[shell, :]
        in_t = np.concatenate(([0.0], glu_time_ms))
        in_c = np.concatenate(([trace[0]], trace))
        step = np.median(np.diff(glu_time_ms)) if glu_bins > 1 else glu_time_ms[0]
        first_zero = min(glu_time_ms[-1] + step, nmda_total)
        if first_zero < nmda_total:
            in_t = np.append(in_t, [first_zero, nmda_total])
            in_c = np.append(in_c, [0.0, 0.0])
        else:
            in_t = np.append(in_t, nmda_total)
            in_c = np.append(in_c, 0.0)

        def glu(t, _t=in_t, _c=in_c):
            return np.interp(t, _t, _c)

        sol = solve_ivp(kinetics.nmda, (0.0, nmda_total), y0, method="RK45",
                        t_eval=out_time, args=(glu,),
                        rtol=cfg.analysis.ODE_RelTol, atol=cfg.analysis.ODE_AbsTol)
        Y = sol.y.T
        selected[shell, :] = open_read(Y)
        desens[shell, :] = desens_read(Y)

    return NmdaResult(out_time, distance_um, selected, desens, conc_uM, glu_time_ms)
