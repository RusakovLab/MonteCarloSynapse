"""Spatial AMPA analysis driver (port of Unified_AMPA_SpaceSR.m).

Given radial free-glutamate histograms (``DistanceFree*.txt`` or an in-memory
matrix), reconstruct the shell concentration in uM and run one of the three
AMPA kinetic schemes independently for every radial shell. Real shell volumes
(e.g. from a diffusion run) may be supplied to avoid the legacy volume formula.
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
class AmpaResult:
    model_type: str
    time_ms: np.ndarray
    distance_um: np.ndarray
    desensitization: np.ndarray
    open_state: np.ndarray
    concentration_uM: np.ndarray
    glu_time_ms: np.ndarray


def load_distance_free(folder):
    files = sorted(glob.glob(os.path.join(folder, "DistanceFree*.txt")))
    if not files:
        raise FileNotFoundError("No DistanceFree*.txt in %s" % folder)
    mats = [np.loadtxt(f) for f in files]
    return np.mean(mats, axis=0)


def reconstruct_concentration(initial_distribution, cfg):
    """Legacy normalized-radius shell volumes (as in the original Unified)."""
    bins = cfg.ampa.BinGlu
    distance_real = np.arange(1, bins + 2) / bins
    outer = (2 * distance_real[1:]) ** 3
    inner = (2 * distance_real[:-1]) ** 3
    volume = (4.0 / 3.0) * np.pi * (outer - inner)
    con = 1000.0 * 1.66e-6
    counts = initial_distribution[:, 1:]
    return con * counts / volume[:, None]


def run_ampa_spatial(cfg, model_type="AMPA", data_folder=None,
                     initial_distribution=None, shell_volume_um3=None, distance_um=None):
    if model_type not in ("AMPA", "AMPA1", "AMPA2"):
        raise ValueError("model_type must be 'AMPA', 'AMPA1', or 'AMPA2'")
    if initial_distribution is None:
        if data_folder is None:
            raise ValueError(
                "run_ampa_spatial needs either initial_distribution=<matrix> or "
                "data_folder=<folder with DistanceFree*.txt> (see run.py ampa --data)")
        initial_distribution = load_distance_free(data_folder)
        # Use the geometry actually written by the diffusion run, so the
        # command line reproduces what the GUI computes; the legacy
        # normalized-radius formula is only a fallback for older folders.
        meta = io_files.load_metadata(data_folder)
        if meta is not None:
            if shell_volume_um3 is None:
                shell_volume_um3 = meta["ShellVolume_um3"]
            if distance_um is None:
                distance_um = meta["OuterRadius_um"]

    bins = cfg.ampa.BinGlu
    time_glu = cfg.ampa.TimeGlu_ms
    horizon = cfg.ampa.TotalTime_ms
    out_step = cfg.ampa.OutputTimeStep_ms
    max_dist = cfg.ampa.MaxDistance_um

    if shell_volume_um3 is not None:
        # Real shell volumes (from diffusion): physically correct concentration.
        counts = initial_distribution[:, 1:]
        concentration = counts_to_uM(counts, shell_volume_um3, cfg.analysis.AvogadroNumber)
    else:
        concentration = reconstruct_concentration(initial_distribution, cfg)
    n_time = concentration.shape[1]

    bin_total = int(round(bins * horizon / time_glu))
    ampa_time_bin = np.linspace(0.0, horizon, bin_total)

    sample_times = np.arange(0.0, horizon + 0.5 * out_step, out_step)
    if sample_times[-1] < horizon:
        sample_times = np.append(sample_times, horizon)

    meta = kinetics.MODELS[model_type]
    rhs = kinetics.RHS[model_type]
    y0 = np.array(meta["y0"], dtype=float)
    method = "BDF" if model_type == "AMPA2" else "RK45"
    atol = cfg.analysis.ODE_AbsTol_Stiff if model_type == "AMPA2" else cfg.analysis.ODE_AbsTol

    desens = np.zeros((bins, sample_times.size))
    openm = np.zeros((bins, sample_times.size))

    for shell in range(bins):
        trace = concentration[shell, :]
        if trace.size < bin_total:
            padded = np.zeros(bin_total)
            padded[:trace.size] = trace
        else:
            padded = trace[:bin_total]

        def glu(t, _tb=ampa_time_bin, _f=padded):
            return np.interp(t, _tb, _f, left=0.0, right=0.0)

        sol = solve_ivp(rhs, (0.0, horizon), y0, method=method,
                        t_eval=sample_times, args=(glu,),
                        rtol=cfg.analysis.ODE_RelTol, atol=atol)
        Y = sol.y.T
        desens[shell, :] = meta["desens"](Y)
        openm[shell, :] = meta["open"](Y)

    if distance_um is None:
        distance_um = np.linspace(cfg.cleft.RadiusofCleft / 1000.0, max_dist, bins)
    glu_time = np.linspace(0.0, time_glu, n_time)
    return AmpaResult(model_type, sample_times, distance_um, desens, openm,
                      concentration, glu_time)
