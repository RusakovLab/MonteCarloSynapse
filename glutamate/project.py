"""Project / results-folder management.

  - save / load the parameter set as JSON (Parameters.json);
  - run one or many trials, writing each trial's outputs into one folder;
  - save receptor (AMPA/AMPA1/AMPA2/NMDA) open/desensitization maps;
  - reload a previous run's files (diffusion AND receptor) so results can be
    viewed without recomputing.
"""
from __future__ import annotations

import dataclasses
import glob
import json
import os

import numpy as np

from .config import (Config, BasicParameters, Balls, ComputationSet, Cleft,
                     SizeCube, Adhesion, AMPAConfig, NMDAConfig, Analysis)
from .diffusion import run_diffusion


_GROUPS = {
    "basic": BasicParameters, "balls": Balls, "computation": ComputationSet,
    "cleft": Cleft, "cube": SizeCube, "adhesion": Adhesion,
    "ampa": AMPAConfig, "nmda": NMDAConfig, "analysis": Analysis,
}


def config_to_dict(cfg):
    return dataclasses.asdict(cfg)


def _build_group(cls, d):
    init_names = {f.name for f in dataclasses.fields(cls) if f.init}
    return cls(**{k: v for k, v in (d or {}).items() if k in init_names})


def config_from_dict(d):
    groups = {name: _build_group(cls, d.get(name, {})) for name, cls in _GROUPS.items()}
    return Config(**groups)


def save_config(cfg, folder):
    os.makedirs(folder, exist_ok=True)
    path = os.path.join(folder, "Parameters.json")
    with open(path, "w") as fh:
        json.dump(config_to_dict(cfg), fh, indent=2)
    return path


def load_config(path):
    if os.path.isdir(path):
        path = os.path.join(path, "Parameters.json")
    with open(path) as fh:
        return config_from_dict(json.load(fh))


def _jsonable(v):
    if isinstance(v, np.ndarray):
        return v.tolist()
    if isinstance(v, (np.floating, np.integer)):
        return v.item()
    return v


def save_balls(result, folder, trial):
    os.makedirs(folder, exist_ok=True)
    data = np.column_stack([result.balls_xyz, result.balls_radius,
                            result.balls_adhesive.astype(int)])
    path = os.path.join(folder, "Balls Trial %02d.txt" % trial)
    np.savetxt(path, data, fmt=["%.3f", "%.3f", "%.3f", "%.3f", "%d"],
               header="X_nm Y_nm Z_nm Radius_nm Adhesive")
    return path


def save_trial(result, folder, trial):
    os.makedirs(folder, exist_ok=True)
    tag = "Trial %02d" % trial
    np.savetxt(os.path.join(folder, "DistanceFree %s.txt" % tag), result.average_free, fmt="%.6g")
    np.savetxt(os.path.join(folder, "DistanceBound %s.txt" % tag), result.average_bound, fmt="%.6g")
    save_balls(result, folder, trial)
    spill = np.column_stack([result.time_ms, result.number_spillover,
                             result.fraction_spillover, result.number_free_spillover,
                             result.fraction_free_spillover])
    np.savetxt(os.path.join(folder, "Spillover %s.txt" % tag), spill, fmt="%.6g",
               header="time_ms number fraction number_free fraction_free")
    meta = {k: _jsonable(v) for k, v in result.metadata.items()}
    meta["Trial"] = trial
    meta["OccupiedFraction"] = result.occupied_fraction
    meta["MeanDiffusion_um2_per_ms"] = result.mean_diffusion
    with open(os.path.join(folder, "SimulationMetadata %s.json" % tag), "w") as fh:
        json.dump(meta, fh, indent=2)
    for label, arr in result.snapshots.items():
        path = os.path.join(folder, "PD %s ms %s.txt" % (label, tag))
        np.savetxt(path, arr.T, fmt="%10.5f %10.5f %10.5f %10.0f")


def save_receptor_result(res, folder, model):
    """Save a receptor result's open/desensitization maps + axes grid."""
    os.makedirs(folder, exist_ok=True)
    np.savetxt(os.path.join(folder, "%s Open State.txt" % model), res.open_state, fmt="%.6g")
    np.savetxt(os.path.join(folder, "%s Desensitization State.txt" % model),
               res.desensitization, fmt="%.6g")
    grid = dict(time_ms=np.asarray(res.time_ms).tolist(),
                distance_um=np.asarray(res.distance_um).tolist(),
                glu_time_ms=np.asarray(res.glu_time_ms).tolist())
    with open(os.path.join(folder, "%s Grid.json" % model), "w") as fh:
        json.dump(grid, fh)


def run_trials(cfg, folder=None, seed=None, save=True, progress_cb=None, frame_cb=None):
    n_trials = max(1, int(cfg.basic.Trials))
    rng = np.random.default_rng(seed)
    if save and folder:
        os.makedirs(folder, exist_ok=True)
        save_config(cfg, folder)
    results = []
    free_stack = []
    bound_stack = []
    for t in range(1, n_trials + 1):
        def _cb(p, _t=t):
            if progress_cb is not None:
                progress_cb(((_t - 1) * 100.0 + p) / n_trials)
        res = run_diffusion(cfg, rng=rng, progress_cb=_cb,
                            frame_cb=(frame_cb if t == 1 else None),
                            graphics_every=cfg.basic.GraphicsUpdateEveryNSteps)
        results.append(res)
        free_stack.append(res.average_free)
        bound_stack.append(res.average_bound)
        if save and folder:
            save_trial(res, folder, t)
    avg_free = np.mean(free_stack, axis=0)
    avg_bound = np.mean(bound_stack, axis=0)
    occ = float(np.mean([r.occupied_fraction for r in results]))
    dmean = float(np.nanmean([r.mean_diffusion for r in results]))
    out = dict(results=results, average_free=avg_free, average_bound=avg_bound,
               metadata=results[0].metadata, occupied_fraction=occ,
               mean_diffusion=dmean, n_trials=n_trials)
    if save and folder:
        np.savetxt(os.path.join(folder, "DistanceFree Mean.txt"), avg_free, fmt="%.6g")
        np.savetxt(os.path.join(folder, "DistanceBound Mean.txt"), avg_bound, fmt="%.6g")
        manifest = dict(n_trials=n_trials, occupied_fraction=occ,
                        mean_diffusion_um2_per_ms=dmean,
                        release_source=cfg.basic.ReleaseSource,
                        particles=cfg.basic.NumberOfParticles,
                        obstacles=cfg.basic.Numbersphere)
        with open(os.path.join(folder, "RunManifest.json"), "w") as fh:
            json.dump(manifest, fh, indent=2)
    return out


def _mean_or_trials(folder, prefix):
    mean = os.path.join(folder, prefix + " Mean.txt")
    if os.path.isfile(mean):
        return np.loadtxt(mean)
    files = sorted(glob.glob(os.path.join(folder, prefix + " Trial*.txt")))
    if not files:
        return None
    return np.mean([np.loadtxt(f) for f in files], axis=0)


def load_results(folder):
    """Load a previous run's files (diffusion + receptor) for re-viewing."""
    out = {}
    pcfg = os.path.join(folder, "Parameters.json")
    out["cfg"] = load_config(pcfg) if os.path.isfile(pcfg) else None
    out["average_free"] = _mean_or_trials(folder, "DistanceFree")
    out["average_bound"] = _mean_or_trials(folder, "DistanceBound")

    metas = sorted(glob.glob(os.path.join(folder, "SimulationMetadata*.json")))
    meta = {}
    if metas:
        with open(metas[0]) as fh:
            raw = json.load(fh)
        for k, v in raw.items():
            meta[k] = np.asarray(v) if isinstance(v, list) else v
    out["metadata"] = meta

    sp_files = sorted(glob.glob(os.path.join(folder, "Spillover Trial*.txt")))
    if sp_files:
        a = np.loadtxt(sp_files[0])
        out["spillover"] = dict(time_ms=a[:, 0], fraction=a[:, 2], fraction_free=a[:, 4])
    else:
        out["spillover"] = None

    snaps = {}
    snap_files = (sorted(glob.glob(os.path.join(folder, "PD * ms Trial 01.txt")))
                  or sorted(glob.glob(os.path.join(folder, "PD * ms Trial*.txt"))))
    for f in snap_files:
        base = os.path.basename(f)
        try:
            label = base.split("PD ", 1)[1].split(" ms", 1)[0]
        except Exception:
            label = base
        a = np.loadtxt(f)
        if a.ndim == 1:
            a = a[None, :]
        snaps[label] = a.T
    out["snapshots"] = snaps

    receptors = {}
    for gridf in sorted(glob.glob(os.path.join(folder, "* Grid.json"))):
        model = os.path.basename(gridf)[:-len(" Grid.json")]
        openf = os.path.join(folder, "%s Open State.txt" % model)
        desf = os.path.join(folder, "%s Desensitization State.txt" % model)
        if not os.path.isfile(openf):
            continue
        with open(gridf) as fh:
            g = json.load(fh)
        receptors[model] = dict(
            open_state=np.loadtxt(openf),
            desensitization=(np.loadtxt(desf) if os.path.isfile(desf) else None),
            time_ms=np.asarray(g["time_ms"]),
            distance_um=np.asarray(g["distance_um"]),
            glu_time_ms=np.asarray(g.get("glu_time_ms", g["time_ms"])))
    out["receptors"] = receptors

    man = os.path.join(folder, "RunManifest.json")
    if os.path.isfile(man):
        with open(man) as fh:
            m = json.load(fh)
        out["n_trials"] = int(m.get("n_trials", 1))
        out["occupied_fraction"] = float(m.get("occupied_fraction", float("nan")))
        out["mean_diffusion"] = float(m.get("mean_diffusion_um2_per_ms", float("nan")))
    else:
        out["n_trials"] = 1
        occ = meta.get("OccupiedFraction", float("nan")) if meta else float("nan")
        out["occupied_fraction"] = float(occ)
        out["mean_diffusion"] = float("nan")
    return out
