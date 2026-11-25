"""Command-line entry point: ``python run.py <command>``.

Examples:
  python run.py diffusion --out results/run1
  python run.py ampa  --data results/run1 --model AMPA2
  python run.py nmda  --data results/run1
  python run.py all   --out results/run1

Equivalent to the installed ``montecarlosynapse`` command and to
``python -m glutamate``. Paths are always given on the command line;
nothing is hard-coded. Pass --seed for reproducibility.
"""
import os
import sys

# Work from any working directory, including an unzipped download.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from glutamate.cli import main

if __name__ == "__main__":
    main()
