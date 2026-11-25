"""Read/write the project's data files (kept separate from computation).

Mirrors the DistanceFree / DistanceBound / SimulationMetadata / PD snapshot
outputs of DiffusionGlutamateBalls.m so the Python and MATLAB pipelines can
exchange data.
"""
from __future__ import annotations

import glob
import json
import os
from datetime import datetime

import numpy as np


def _timestamp() -> str:
    return datetime.now().strftime("%Y%m%d-%H%M%S-%f")[:-3]


def save_receptor_input(result, folder: str) -> dict:
    """Write DistanceFree/DistanceBound matrices + metadata JSON. Returns paths."""
    os.makedirs(folder, exist_ok=True)
    run_id = _timestamp()
    free_path = os.path.join(folder, f"DistanceFree {run_id}.txt")
    bound_path = os.path.join(folder, f"DistanceBound {run_id}.txt")
    np.savetxt(free_path, result.average_free, fmt="%.6g")
    np.savetxt(bound_path, result.average_bound, fmt="%.6g")

    meta = dict(result.metadata)
    meta["RunID"] = run_id
    meta["DistanceFreeFile"] = os.path.basename(free_path)
    meta["DistanceBoundFile"] = os.path.basename(bound_path)
    meta_path = os.path.join(folder, f"SimulationMetadata {run_id}.json")
    with open(meta_path, "w") as fh:
        json.dump({k: (v.tolist() if isinstance(v, np.ndarray) else v)
                   for k, v in meta.items()}, fh, indent=2)
    return {"free": free_path, "bound": bound_path, "metadata": meta_path}


def save_snapshots(result, folder: str) -> list:
    os.makedirs(folder, exist_ok=True)
    paths = []
    for label, arr in result.snapshots.items():
        path = os.path.join(folder, f"PD {label} ms {_timestamp()}.txt")
        # columns: X Y Z State  (rows = particles)
        np.savetxt(path, arr.T, fmt="%10.5f %10.5f %10.5f %10.0f")
        paths.append(path)
    return paths


def load_distance_matrix(path: str) -> np.ndarray:
    return np.loadtxt(path)


_GEOMETRY_KEYS = ("ShellVolume_um3", "TimePoints_ms", "OuterRadius_um")


def load_metadata(folder: str, distance_file: str = None):
    """Return a run's geometry metadata as a dict, or None if there is none.

    Looks for the JSON written by this package (preferring the file whose run
    id matches ``distance_file``), then for MATLAB ``SimulationMetadata*.mat``
    written by the original implementation. Only the keys needed by the
    receptor drivers are returned, as float arrays.
    """
    run_id = None
    if distance_file:
        base = os.path.basename(distance_file)
        run_id = (base.replace("DistanceFree", "").replace(".txt", "").strip()
                  or None)

    candidates = []
    if run_id:
        candidates.append(os.path.join(folder, "SimulationMetadata %s.json" % run_id))
    candidates += sorted(glob.glob(os.path.join(folder, "SimulationMetadata*.json")))
    for path in candidates:
        if not os.path.isfile(path):
            continue
        try:
            with open(path) as fh:
                raw = json.load(fh)
        except (OSError, ValueError):
            continue
        if all(k in raw for k in _GEOMETRY_KEYS):
            return {k: np.asarray(raw[k], float).ravel() for k in _GEOMETRY_KEYS}

    mat_candidates = []
    if run_id:
        mat_candidates.append(os.path.join(folder, "SimulationMetadata %s.mat" % run_id))
    mat_candidates += sorted(glob.glob(os.path.join(folder, "SimulationMetadata*.mat")))
    for path in mat_candidates:
        if not os.path.isfile(path):
            continue
        try:
            from scipy.io import loadmat
            m = loadmat(path, squeeze_me=True, struct_as_record=False)["Metadata"]
            return {k: np.asarray(getattr(m, k), float).ravel() for k in _GEOMETRY_KEYS}
        except Exception:
            continue
    return None
