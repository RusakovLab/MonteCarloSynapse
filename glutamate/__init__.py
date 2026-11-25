"""MonteCarloSynapse -- glutamate diffusion and receptor kinetics.

Layered Python package:
  config          - single parameter source for every stage
  kinetics        - AMPA/AMPA1/AMPA2/NMDA ODE right-hand sides
  concentration   - molecule counts -> uM
  diffusion       - Monte Carlo glutamate simulation
  receptors_ampa  - spatial AMPA driver (per radial shell)
  receptors_nmda  - spatial NMDA driver (per radial shell)
  project         - parameter save/load, multi-trial runs, results folders
  plotting        - figures, separate from computation
  io_files        - read/write radial-histogram / metadata / snapshot files
  video           - optional MP4/GIF export of the 3D field
"""
from .config import Config, default_config, quick_config

__all__ = ["Config", "default_config", "quick_config"]
__version__ = "0.2.0"
