
# NatureComm
## Computational Model of Glutamate Diffusion and Receptor Activation

[![MATLAB 2022b](https://img.shields.io/badge/MATLAB-2022b-blue.svg)](https://uk.mathworks.com/products/matlab.html)
[![Windows](https://img.shields.io/badge/OS-Windows-blue.svg)](https://www.microsoft.com/windows/)
[![License](https://img.shields.io/badge/License-UCL-red.svg)](https://www.ucl.ac.uk/legal-services/licensing-and-intellectual-property)

This repository contains MATLAB code to simulate **glutamate diffusion** in the extracellular space and the subsequent **activation of NMDA and AMPA receptors**, incorporating biophysically realistic kinetics. The models support a study submitted to *Nature Communications*.

---

## 1. System Requirements

* **Operating System:** Windows 10/11 (64-bit recommended)
* **MATLAB:** R2022b (required for `.m` files)
* **Processor (tested):**
    * Intel Core i9-12900K @ 3.2 GHz (16 Cores, 24 Threads)
    * UCL High-Performance Computing Cluster (Slurm)
* **RAM:** 16 GB or more recommended
* **Disk Space:** ~1 GB for full output set

---

## 2. Getting Started

### For MATLAB Users

1.  **Clone** or download this repository.
2.  Open **MATLAB R2022b**.
3.  Add the main directory to your MATLAB path.
4.  Edit `InputParametersSR.m` and `statisticSR.txt` to configure your simulation parameters.

### For Non-MATLAB Users (Demo Only)

1.  Navigate to the `Demo/` directory.
2.  Unzip and install the **MATLAB Runtime R2022b**:
    * `MATLAB_Runtime_R2022b_Update_10_win64.zip`
    * Run `setup.exe` inside the ZIP archive.
3.  Launch `DiffusionGlutamateBalls.exe` to run the demo.

---

## 3. Demo Overview

The standalone executable `DiffusionGlutamateBalls.exe` provides a demonstration of glutamate diffusion for 100 particles.

### Estimated Runtimes

| System                  | Task              | Particles | Time      |
| :---------------------- | :---------------- | :-------- | :-------- |
| Local PC (i9-12900K)    | Demo (EXE)        | 100       | ~20 mins  |
| Local PC (i9-12900K)    | Full MATLAB Sim   | 5000      | ~8 hours  |
| UCL HPC Cluster (Slurm) | Full MATLAB Sim   | 5000      | ~2 hours  |

---

## 4. Usage Instructions

### Step 1: Configure Simulation Parameters

* **`InputParametersSR.m`**: Defines core simulation settings like domain size, molecule count, and release radius.
* **`statisticSR.txt`**: Sets statistical parameters, including the number of trials and adhesion probabilities.

### Step 2: Run Glutamate Diffusion

Execute the primary diffusion simulation script:

```matlab
DiffusionGlutamateBalls.m
```

This script generates the following output files in your working directory:

| File Name                | Description                                         |
| :----------------------- | :-------------------------------------------------- |
| `Balls distribution *.txt` | Coordinates and radii of spherical obstacles        |
| `DistanceFree *.txt`     | Space-time distribution of free glutamate           |
| `DistanceBound *.txt`    | Space-time distribution of adhesion-bound glutamate |
| `PD 0.1/0.3/1/3 ms *.txt` | Molecular positions at specific times (0=free, 1=bound) |

### Step 3: Simulate Receptor Dynamics

#### NMDA Receptor Activation

Run the following script to simulate NMDA receptor kinetics:

```matlab
NMDA_SpaceSR.m
```

#### AMPA Receptor Activation (Unified)

Use this script to run all AMPA receptor models:

```matlab
RunAmpa.m
Unified_AMPA_SpaceSR.m
```

This script calls `Unified_AMPA_SpaceSR.m`, supporting three kinetic models:

| Model          | File      | Description                            |
| :------------- | :-------- | :------------------------------------- |
| Patneau-Mayer  | `AMPA.m`  | 6-state model (Neuron 1991, Destexhe 1996) |
| Raman-Trussell | `AMPA1.m` | Dual open-state model (Neuron 1992)    |
| GluR1-based    | `AMPA2.m` | 12-state ligand binding model          |

To switch between simulating with **free** versus **bound** glutamate molecules, modify the `filePattern` in `Unified_AMPA_SpaceSR.m`:

```matlab
filePattern = fullfile(myFolder, 'DistanceFree*.txt');  % Use 'DistanceBound*.txt' for adhered molecules
```

### Step 4: Visualize Results

For additional visualization of spatial glutamate fields and NMDA activation:

```matlab
PlotDataControl.m
```

Note that each `*_SpaceSR.m` script automatically generates plots of receptor state maps and time series.

---

## 5. License

This MATLAB is used under the **University College London academic use license**. 

---

## 6. Reproducing Manuscript Results

To reproduce the quantitative results from the submitted manuscript:

1.  Run DiffusionGlutamateBalls.m using particle counts of N = 1000 and N = 5000, specified in the **statisticSR.txt** file along with the astrocyte density probability (ProbabilityofAstrocytes = 0.1).
2.  Simulate receptor activation using `RunAmpa.m` or `NMDA_SpaceSR.m`. Test data can be found in the `DataforAMPAandNMDATesting` folder on the GitHub repository.
3.  Compare your results using `PlotDataControl.m` and the visual output scripts embedded in `*_SpaceSR.m`.

All necessary simulation parameters and example outputs are included in this repository.

---

## Citation

If you use this code in your research, please cite the associated paper submitted to *Nature Communications*.

---

## Repository & Contact

* **GitHub Repository:** [https://github.com/RusakovLab/NatureComm](https://github.com/RusakovLab/NatureComm)
* **Corresponding Author:** Leonid Savchenko (UCL Institute of Neurology)
