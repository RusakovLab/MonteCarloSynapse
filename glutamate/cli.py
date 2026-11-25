"""Command-line interface for the glutamate package.

Examples:
  montecarlosynapse diffusion --out results/run1
  montecarlosynapse ampa  --data results/run1 --model AMPA2
  montecarlosynapse nmda  --data results/run1
  montecarlosynapse all   --out results/run1

Paths are always given on the command line; nothing is hard-coded.

All parameters come from glutamate.config.default_config(); edit that or pass
--seed for reproducibility.
"""
from __future__ import annotations

import argparse

import numpy as np

from .config import default_config, quick_config
from . import diffusion, receptors_ampa, receptors_nmda, io_files


def _cfg(args):
    """Full-scale parameters, or the small --quick set for a fast trial run."""
    return quick_config() if getattr(args, "quick", False) else default_config()


def cmd_diffusion(args):
    cfg = _cfg(args)
    rng = np.random.default_rng(args.seed)
    print("Running diffusion Monte Carlo...")
    res = diffusion.run_diffusion(cfg, rng=rng, progress=True)
    print(f"Occupancy: occupied={100*res.occupied_fraction:.3f}%  "
          f"free={100*res.free_fraction:.3f}%  SE={100*res.occupancy_se:.3f}%")
    print(f"Effective diffusion: mean={res.mean_diffusion:.4g}  std={res.std_diffusion:.4g}")
    if args.out:
        paths = io_files.save_receptor_input(res, args.out)
        io_files.save_snapshots(res, args.out)
        print("Saved:", paths)
    return res


def cmd_ampa(args):
    cfg = _cfg(args)
    res = receptors_ampa.run_ampa_spatial(cfg, model_type=args.model, data_folder=args.data)
    print(f"{args.model}: open max={res.open_state.max():.4f}  "
          f"desens max={res.desensitization.max():.4f}")
    return res


def cmd_nmda(args):
    cfg = _cfg(args)
    res = receptors_nmda.run_nmda_spatial(cfg, data_folder=args.data)
    print(f"NMDA: open max={res.open_state.max():.4f}")
    return res


def cmd_all(args):
    args.data = args.out
    cmd_diffusion(args)
    for m in ("AMPA", "AMPA1", "AMPA2"):
        args.model = m
        cmd_ampa(args)
    cmd_nmda(args)


def main():
    p = argparse.ArgumentParser(description="Glutamate diffusion + receptor kinetics")
    sub = p.add_subparsers(dest="cmd", required=True)

    quick_help = ("small, fast configuration for trying the workflow "
                  "(80 particles, 60 obstacles, 0.5 ms) instead of the "
                  "full-scale defaults, which take a long time")

    pd = sub.add_parser("diffusion"); pd.add_argument("--out", default=None)
    pd.add_argument("--quick", action="store_true", help=quick_help)
    pd.add_argument("--seed", type=int, default=None); pd.set_defaults(func=cmd_diffusion)

    pa = sub.add_parser("ampa"); pa.add_argument("--data", required=True)
    pa.add_argument("--model", default="AMPA", choices=["AMPA", "AMPA1", "AMPA2"])
    pa.add_argument("--quick", action="store_true", help=quick_help)
    pa.add_argument("--seed", type=int, default=None)
    pa.set_defaults(func=cmd_ampa)

    pn = sub.add_parser("nmda"); pn.add_argument("--data", required=True)
    pn.add_argument("--quick", action="store_true", help=quick_help)
    pn.add_argument("--seed", type=int, default=None); pn.set_defaults(func=cmd_nmda)

    pl = sub.add_parser("all"); pl.add_argument("--out", default="results/run")
    pl.add_argument("--quick", action="store_true", help=quick_help)
    pl.add_argument("--seed", type=int, default=None); pl.set_defaults(func=cmd_all)

    args = p.parse_args()
    args.func(args)
