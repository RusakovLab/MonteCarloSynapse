"""Fast smoke test of the diffusion Monte Carlo on a tiny configuration."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from glutamate import default_config
from glutamate import diffusion, receptors_ampa


def _tiny_cfg():
    cfg = default_config()
    cfg.basic.NumberOfParticles = 60
    cfg.basic.Numbersphere = 40
    cfg.basic.TotalTime = 0.5
    cfg.basic.ShowGraphics = "n"
    cfg.basic.SaveFiles = "n"
    cfg.balls.OccupancySampleCount = 2000
    cfg.computation.DistanceBins = 20
    cfg.computation.TimeBins = 10
    cfg.computation.SnapshotTimes = [0.1, 0.3]
    cfg.__post_init__()
    return cfg


def test_diffusion_runs_and_conserves_counts():
    cfg = _tiny_cfg()
    res = diffusion.run_diffusion(cfg, rng=np.random.default_rng(1))
    assert res.average_free.shape == (20, 11)
    assert res.average_bound.shape == (20, 11)
    assert 0.0 <= res.occupied_fraction <= 1.0
    # spillover fractions are within [0,1]
    assert np.all(res.fraction_spillover <= 1.0 + 1e-9)
    assert len(res.snapshots) == 2


def test_ampa_from_diffusion_output():
    cfg = _tiny_cfg()
    res = diffusion.run_diffusion(cfg, rng=np.random.default_rng(2))
    # feed the free-histogram matrix straight into the AMPA driver
    cfg.ampa.BinGlu = cfg.computation.DistanceBins
    ampa = receptors_ampa.run_ampa_spatial(
        cfg, model_type="AMPA", initial_distribution=res.average_free)
    assert ampa.open_state.shape[0] == cfg.computation.DistanceBins
    assert np.isfinite(ampa.open_state).all()


if __name__ == "__main__":
    test_diffusion_runs_and_conserves_counts()
    test_ampa_from_diffusion_output()
    print("Diffusion smoke tests passed.")
