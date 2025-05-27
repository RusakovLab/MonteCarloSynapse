Absolutely! Adding a badge for MATLAB 2022b is a nice touch. For a license badge, it's common to state the type of license (e.g., MIT, Apache 2.0). Since your project is associated with UCL and a Nature Communications paper, you might not have an open-source license defined yet. I'll add a placeholder that indicates it's under "Proprietary" or "UCL License" for now, and you can change it if you decide on a specific open-source license later.

Here's the updated `README.md` with the badges at the top:

```markdown
# NatureComm
## Computational Model of Glutamate Diffusion and Receptor Activation

[![MATLAB 2022b](https://img.shields.io/badge/MATLAB-2022b-blue.svg)](https://uk.mathworks.com/products/matlab.html)
[![License](https://img.shields.io/badge/License-UCL-red.svg)](https://www.ucl.ac.uk/legal-services/licensing-and-intellectual-property)

This repository provides MATLAB code for simulating the **diffusion of glutamate molecules** in the extracellular space and modeling the subsequent **activation of NMDA and AMPA receptors**. The simulations incorporate biophysically detailed kinetics, reflecting conditions of glutamate release from astrocytes and its interaction with cellular structures. These models support a study submitted to *Nature Communications*.

---

## Workflow

### 1. Configure Parameters

Before running simulations, adjust the following files to match your experimental setup:

* `InputParametersSR.m`: Defines simulation parameters such as domain size, release volume, and particle count.
* `statisticSR.txt`: Sets statistical parameters like the number of trials and adhesion probabilities.

---

### 2. Run Glutamate Diffusion Simulation

Execute `DiffusionGlutamateBalls.m` to simulate glutamate diffusion and adhesion within a 3D extracellular environment. This script generates several output files in your working directory, each providing different insights into the simulation results:

* **`Balls distribution <timestamp>.txt`**: Contains the 3D coordinates and radii of the spherical obstacles, which represent astrocyte proxies.
* **`DistanceFree <timestamp>.txt`**: Records the space-time distribution of free glutamate molecules in concentric bins.
* **`DistanceBound <timestamp>.txt`**: Shows the space-time distribution of glutamate molecules adhered to obstacles.
* **`PD 0.1 ms <timestamp>.txt`**, **`PD 0.3 ms <timestamp>.txt`**, **`PD 1 ms <timestamp>.txt`**, **`PD 3 ms <timestamp>.txt`**: These files provide the coordinates of all glutamate molecules at specific time points (0.1 ms, 0.3 ms, 1 ms, and 3 ms). The last column indicates whether the molecule is free (0) or bound (1).

These output files offer a comprehensive view of glutamate movement and interaction in both space and time, serving as crucial input for the subsequent receptor models.

---

### 3. Compute Receptor Activation Dynamics

#### NMDA Receptors

Simulate NMDA receptor kinetics in response to the generated glutamate profiles by running:

```matlab
NMDA_SpaceSR.m
```

#### AMPA Receptors (Unified Interface)

All AMPA receptor models are accessible through a single script:

```matlab
RunAmpa.m
```

This script calls `Unified_AMPA_SpaceSR.m`, which supports multiple AMPA receptor models:

| Model Name       | Files    | Description                                   |
| :--------------- | :------- | :-------------------------------------------- |
| Patneau–Mayer    | `AMPA.m` | 6-state model (Neuron 1991, Destexhe 1996)    |
| Raman–Trussell   | `AMPA1.m` | Dual-open-state model (Neuron 1992)           |
| GluR1-based      | `AMPA2.m` | 12-state model with multiple ligand bindings |

To select the glutamate source (free or bound molecules) for AMPA simulations, modify the `filePattern` variable within `Unified_AMPA_SpaceSR.m`:

```matlab
filePattern = fullfile(myFolder, 'DistanceFree*.txt');  % for free molecules
% or
filePattern = full(myFolder, 'DistanceBound*.txt'); % for adhered molecules
```

---

### 4. Visualization

Each `*_SpaceSR.m` script automatically generates:

* **Contour plots** of desensitization and open-state probability.
* **Time-resolved graphs** for different bins and distances.

For additional visualization of spatial glutamate fields and NMDA activation, use:

```matlab
PlotDataControl.m
```

---

## Included Files

### Diffusion
* `DiffusionGlutamateBalls.m`
* `InputParametersSR.m`
* `statisticSR.txt`

### NMDA Model
* `NMDA.m`
* `NMDA_SpaceSR.m`

### AMPA Models
* `AMPA.m` – Patneau–Mayer
* `AMPA1.m` – Raman–Trussell
* `AMPA2.m` – GluR1 multistate
* `Unified_AMPA_SpaceSR.m` – Core processor
* `RunAmpa.m` – Unified entry point

### Visualization
* `PlotDataControl.m`

### Example Output Files
* `Balls distribution *.txt`
* `DistanceFree *.txt`
* `DistanceBound *.txt`
* `PD 0.1 ms *.txt`
* `PD 0.3 ms *.txt`
* `PD 1 ms *.txt`
* `PD 3 ms *.txt`

---

## Run Without MATLAB

A standalone version of the glutamate diffusion simulator is provided in the `Demo/` folder.

**Files:**
* `DiffusionGlutamateBalls.exe`: Precompiled executable.
* `MATLAB_Runtime_R2022b_Update_10_win64.zip`: Required MATLAB Runtime.

**Instructions:**
1.  Unzip and install the MATLAB Runtime (run `setup.exe`).
2.  Execute `DiffusionGlutamateBalls.exe` to launch the simulation.

---

## Citation

If you use this code in your research, please cite the associated paper submitted to *Nature Communications*.
```
