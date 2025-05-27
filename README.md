
# NatureComm  
## Computational Model of Glutamate Diffusion and Receptor Activation  
[![MATLAB 2022b](https://img.shields.io/badge/MATLAB-2022b-blue.svg)](https://uk.mathworks.com/products/matlab.html)  
[![License](https://img.shields.io/badge/License-UCL-red.svg)](https://uk.mathworks.com/academia/tah-portal/university-college-london-649021.html)

This repository contains MATLAB code to simulate **glutamate diffusion in the extracellular space** and the subsequent **activation of NMDA and AMPA receptors**, incorporating biophysically realistic kinetics. The models support a study submitted to *Nature Communications*.

---

## 1. System Requirements

- **Operating System:** Windows 10/11 (64-bit recommended)
- **MATLAB:** R2022b (required for running `.m` files)
- **Processor (tested):**  
  - Intel Core i9-12900K @ 3.2 GHz (16 Cores, 24 Threads)  
  - UCL High-Performance Computing Cluster (Slurm)
- **RAM:** 16 GB or more recommended
- **Disk Space:** ~1 GB for full output set

---

## 2. Installation Guide

### MATLAB Users

1. Clone or download this repository.
2. Open MATLAB R2022b.
3. Add the main directory to your MATLAB path.
4. Edit `InputParametersSR.m` and `statisticSR.txt` to set parameters.

### Non-MATLAB Users (Demo Only)

1. Go to the `Demo/` directory.
2. Unzip and install:
   - `MATLAB_Runtime_R2022b_Update_10_win64.zip`
   - Run `setup.exe` inside the ZIP.
3. Launch:
   - `DiffusionGlutamateBalls.exe`

---

## 3. Demo

The standalone executable `DiffusionGlutamateBalls.exe` demonstrates glutamate diffusion for 100 particles.

### Runtime (Benchmark):
| System                        | Task                | Particles | Time      |
|------------------------------|---------------------|-----------|-----------|
| Local PC (i9-12900K)         | Demo (EXE)          | 100       | ~20 mins  |
| Local PC (i9-12900K)         | Full MATLAB Sim     | 5000      | ~8 hours  |
| UCL HPC Cluster              | Full MATLAB Sim     | 5000      | ~2 hours  |

---

## 4. Instructions for Use

### Step 1: Configure Simulation

- `InputParametersSR.m`: Defines domain size, molecule count, release radius, etc.
- `statisticSR.txt`: Defines number of trials and adhesion probabilities.

### Step 2: Run Diffusion

```matlab
DiffusionGlutamateBalls.m
````

#### Output Files:

| File Name                  | Description                                         |
| -------------------------- | --------------------------------------------------- |
| `Balls distribution *.txt` | Coordinates and radii of spherical obstacles        |
| `DistanceFree *.txt`       | Free glutamate distribution over space and time     |
| `DistanceBound *.txt`      | Adhesion-bound glutamate distribution               |
| `PD 0.1/0.3/1/3 ms *.txt`  | Molecular positions at different times (bound/free) |

### Step 3: Receptor Dynamics

#### NMDA Activation

```matlab
NMDA_SpaceSR.m
```

#### AMPA Activation (Unified)

```matlab
RunAmpa.m
```

This calls `Unified_AMPA_SpaceSR.m` and supports 3 kinetic models:

| Model          | File      | Description                          |
| -------------- | --------- | ------------------------------------ |
| Patneau-Mayer  | `AMPA.m`  | 6-state (Neuron 1991, Destexhe 1996) |
| Raman-Trussell | `AMPA1.m` | Dual open state (Neuron 1992)        |
| GluR1-based    | `AMPA2.m` | 12-state ligand model                |

To switch between free vs bound glutamate:

```matlab
filePattern = fullfile(myFolder, 'DistanceFree*.txt');  % or DistanceBound
```

### Step 4: Visualization

```matlab
PlotDataControl.m
```

Also, each `*_SpaceSR.m` script automatically plots receptor state maps and time series.

---

## 5. License

This software is provided under the **University College London academic use license**.
For use beyond academic or research purposes, contact the author.

---

## 6. Reproducing Manuscript Results

To reproduce quantitative results from the manuscript:

1. Run `DiffusionGlutamateBalls.m` with:

   * `N = 5000`, `TimeSteps = 1000`
2. Simulate receptor activation with:

   * `RunAmpa.m` or `NMDA_SpaceSR.m`
3. Compare results using:

   * `PlotDataControl.m` and visual output scripts embedded in `*_SpaceSR.m`

All required simulation parameters and example outputs are included in the repository.

---

## Citation

If you use this code, please cite the associated paper submitted to *Nature Communications*.

---

## Repository

📦 GitHub: [https://github.com/RusakovLab/NatureComm](https://github.com/RusakovLab/NatureComm)
📄 Corresponding author: Lesha Savchenko (UCL Institute of Neurology)

```

---

Let me know if you want me to generate a `LICENSE` file or a zipped version of the entire structure for submission.
```
