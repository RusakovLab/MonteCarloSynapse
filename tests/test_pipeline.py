"""End-to-end check of the file-based pipeline used by the command line.

A diffusion run is saved to a folder, then the AMPA and NMDA drivers are run
from that folder alone -- the path a user takes with
``python run.py diffusion --out ...`` followed by ``run.py ampa|nmda --data``.
"""
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from glutamate import quick_config
from glutamate import diffusion, io_files, project, receptors_ampa, receptors_nmda


def _run(tmp):
    cfg = quick_config()
    cfg.basic.TotalTime = 0.2
    cfg.basic.NumberOfParticles = 40
    cfg.basic.Numbersphere = 30
    cfg.balls.OccupancySampleCount = 1000
    cfg.computation.__post_init__()
    cfg.__post_init__()
    res = diffusion.run_diffusion(cfg, rng=np.random.default_rng(0))
    io_files.save_receptor_input(res, tmp)
    return cfg, res


def test_metadata_round_trip():
    with tempfile.TemporaryDirectory() as tmp:
        cfg, res = _run(tmp)
        meta = io_files.load_metadata(tmp)
        assert meta is not None, "SimulationMetadata JSON not found"
        assert np.allclose(meta["ShellVolume_um3"], res.metadata["ShellVolume_um3"])
        assert np.allclose(meta["OuterRadius_um"], res.metadata["OuterRadius_um"])


def test_receptor_drivers_read_the_run_geometry():
    """Driving from the folder must match driving from the in-memory arrays."""
    with tempfile.TemporaryDirectory() as tmp:
        cfg, res = _run(tmp)

        from_folder = receptors_ampa.run_ampa_spatial(cfg, model_type="AMPA",
                                                      data_folder=tmp)
        in_memory = receptors_ampa.run_ampa_spatial(
            cfg, model_type="AMPA", initial_distribution=res.average_free,
            shell_volume_um3=res.metadata["ShellVolume_um3"],
            distance_um=res.metadata["OuterRadius_um"])
        assert np.allclose(from_folder.concentration_uM, in_memory.concentration_uM)
        assert np.allclose(from_folder.distance_um, in_memory.distance_um)

        n_folder = receptors_nmda.run_nmda_spatial(cfg, data_folder=tmp)
        n_memory = receptors_nmda.run_nmda_spatial(
            cfg, initial_distribution=res.average_free,
            shell_volume_um3=res.metadata["ShellVolume_um3"],
            glu_time_ms=res.metadata["TimePoints_ms"],
            distance_um=res.metadata["OuterRadius_um"])
        assert np.allclose(n_folder.open_state, n_memory.open_state)
        assert np.allclose(n_folder.desensitization, n_memory.desensitization)


def test_nmda_results_save_and_reload():
    with tempfile.TemporaryDirectory() as tmp:
        cfg, res = _run(tmp)
        n = receptors_nmda.run_nmda_spatial(cfg, data_folder=tmp)
        project.save_receptor_result(n, tmp, "NMDA")
        assert os.path.isfile(os.path.join(tmp, "NMDA Open State.txt"))
        assert os.path.isfile(os.path.join(tmp, "NMDA Desensitization State.txt"))
        back = project.load_results(tmp)["receptors"]["NMDA"]
        assert np.allclose(back["open_state"], n.open_state, rtol=1e-5, atol=1e-9)
        assert np.allclose(back["desensitization"], n.desensitization,
                           rtol=1e-5, atol=1e-9)


if __name__ == "__main__":
    test_metadata_round_trip()
    test_receptor_drivers_read_the_run_geometry()
    test_nmda_results_save_and_reload()
    print("Pipeline tests passed.")
