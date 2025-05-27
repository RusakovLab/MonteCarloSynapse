# NatureComm
Computational Model of Glutamate Diffusion and Receptor Activation

## Overview

This repository contains MATLAB code and associated files used for simulating glutamate diffusion in the extracellular space and its interaction with both NMDA and AMPA receptors. The results support a study submitted to *Nature Communications*.

The model captures the release of glutamate from astrocytes, diffusion through a 3D domain filled with cellular structures, and subsequent receptor activation using detailed biophysical kinetics.

## Components

### Glutamate Diffusion
- **`FixBallsAstrogliaSurfaceRelease.m`**: Simulates 3D Brownian motion of glutamate among spherical obstacles representing astrocytes, including adhesion effects.
- **`InputParametersSR.m`**: Sets simulation parameters (particle number, domain size, surface properties).
- **`statisticSR.txt`**: Contains statistical settings like adhesion probability and number of trials.

### NMDA Receptor Dynamics
- **`NMDA.m`**, **`NMDA_SpaceSR.m`**: Implement kinetic models of NMDA receptor activation, including desensitization and glutamate-triggered transitions.

### AMPA Receptor Dynamics *(NEW)*
The following six files simulate AMPA receptor kinetics under different biophysical scenarios and kinetic schemes:

- **`AMPA.m`**, **`AMPA_SpaceSR.m`**: Implements a 6-state model from Patneau & Mayer (Neuron, 1991) and Destexhe et al. (1996). Used to assess desensitization and open probability across distances.
- **`AMPA1.m`**, **`AMPA1_SpaceSR.m`**: Implements an alternative AMPA model based on Raman & Trussell (Neuron, 1992), with dual open states and fast desensitization.
- **`AMPA2.m`**, **`AMPA2_SpaceSR.m`**: A more detailed multistate GluR1-based AMPA model with 12 dynamic states capturing multiple ligand-bound, open, and desensitized configurations.

Each `*_SpaceSR.m` script processes glutamate distributions and numerically solves the receptor kinetics using `ode45`, outputting receptor states over time and distance.

### Visualization
- **`PlotDataControl.m`**: Plots 3D glutamate distributions and NMDA receptor states.
- AMPA scripts generate additional plots (contour, line) showing receptor open probabilities and desensitization over time and space.

## Usage

1. **Parameter Setup**
   - Edit `InputParametersSR.m` and `statisticSR.txt` for diffusion parameters.
   - Define glutamate distributions as `DistanceFree*.txt` in the working directory.

2. **Run Diffusion Simulation**
   - Execute `FixBallsAstrogliaSurfaceRelease.m` to simulate glutamate release and space sampling.

3. **Receptor Activation Modeling**
   - For NMDA: Run `NMDA_SpaceSR.m`.
   - For AMPA:
     - Use `AMPA_SpaceSR.m` (Patneau-Mayer model)
     - Use `AMPA1_SpaceSR.m` (Raman-Trussell model)
     - Use `AMPA2_SpaceSR.m` (GluR1 multistate model)

4. **Visualization**
   - Use built-in figures in each script or run `PlotDataControl.m` to visualize receptor states and glutamate fields.

## Included Files

- `FixBallsAstrogliaSurfaceRelease.m`
- `InputParametersSR.m`
- `statisticSR.txt`
- `NMDA.m`, `NMDA_SpaceSR.m`
- `AMPA.m`, `AMPA_SpaceSR.m`
- `AMPA1.m`, `AMPA1_SpaceSR.m`
- `AMPA2.m`, `AMPA2_SpaceSR.m`
- `PlotDataControl.m`

## Citation

If you use this code in your research, please cite the associated paper submitted to *Nature Communications*.
```


