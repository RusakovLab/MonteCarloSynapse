# NatureComm: Glutamate Diffusion and Receptor Activation Model

[![MATLAB 2022b](https://img.shields.io/badge/MATLAB-2022b-blue.svg)](https://uk.mathworks.com/products/matlab.html)
[![Windows](https://img.shields.io/badge/OS-Windows-blue.svg)](https://www.microsoft.com/windows/)
[![License](https://img.shields.io/badge/License-UCL-red.svg)][(https://www.ucl.ac.uk/legal-services/licensing-and-intellectual-property)](https://uk.mathworks.com/academia/tah-portal/university-college-london-649021.html)

> **Computational model for simulating glutamate diffusion in extracellular space and subsequent NMDA/AMPA receptor activation with biophysically realistic kinetics.**

This repository contains MATLAB code supporting a study submitted to *Nature Communications*.

## 🚀 Quick Start

### For MATLAB Users
```bash
git clone https://github.com/RusakovLab/NatureComm.git
cd NatureComm
```
1. Open MATLAB R2022b
2. Add the main directory to your MATLAB path
3. Configure `InputParametersSR.m` and `statisticSR.txt`
4. Run `DiffusionGlutamateBalls.m`

### For Non-MATLAB Users (Demo)
1. Navigate to `Demo/` directory
2. Install MATLAB Runtime R2022b from `MATLAB_Runtime_R2022b_Update_10_win64.zip`
3. Run `DiffusionGlutamateBalls.exe`

## 📋 System Requirements

| Component | Requirement |
|-----------|-------------|
| **OS** | Windows 10/11 (64-bit recommended) |
| **MATLAB** | R2022b (required for `.m` files) |
| **RAM** | 16 GB+ recommended |
| **Storage** | ~1 GB for full output |
| **Processor** | Tested on Intel i9-12900K & UCL HPC Cluster |

## ⏱️ Performance Benchmarks

| System | Configuration | Particles | Runtime |
|--------|---------------|-----------|---------|
| Local PC (i9-12900K) | Demo executable | 100 | ~20 min |
| Local PC (i9-12900K) | Full MATLAB simulation | 5,000 | ~8 hours |
| UCL HPC Cluster | Full MATLAB simulation | 5,000 | ~2 hours |

## 🔬 Model Components

### Glutamate Diffusion
- **Primary script**: `DiffusionGlutamateBalls.m`
- Simulates molecular diffusion in extracellular space
- Accounts for spherical obstacles and adhesion binding

### Receptor Models

#### NMDA Receptors
- **Script**: `NMDA_SpaceSR.m`
- Biophysically realistic kinetics

#### AMPA Receptors
- **Main script**: `RunAmpa.m` → calls `Unified_AMPA_SpaceSR.m`
- **Three kinetic models**:
  - **Patneau-Mayer** (`AMPA.m`): 6-state model
  - **Raman-Trussell** (`AMPA1.m`): Dual open-state model  
  - **GluR1-based** (`AMPA2.m`): 12-state ligand binding model

## 📊 Output Files

| File Pattern | Content |
|--------------|---------|
| `Balls distribution *.txt` | Spherical obstacle coordinates & radii |
| `DistanceFree *.txt` | Free glutamate space-time distribution |
| `DistanceBound *.txt` | Bound glutamate space-time distribution |
| `PD 0.1/0.3/1/3 ms *.txt` | Molecular positions at specific timepoints |

## 🛠️ Usage Workflow

### Step 1: Configure Parameters
Edit simulation settings in:
- `InputParametersSR.m` - Core parameters (domain size, molecule count, release radius)
- `statisticSR.txt` - Statistical parameters (trials, adhesion probabilities)

### Step 2: Run Diffusion Simulation
```matlab
DiffusionGlutamateBalls.m
```

### Step 3: Simulate Receptor Dynamics
```matlab
% For NMDA receptors
NMDA_SpaceSR.m

% For AMPA receptors (all models)
RunAmpa.m
```

### Step 4: Visualize Results
```matlab
PlotDataControl.m
```

## ⚙️ Configuration Options

### Switching Between Free/Bound Molecules
In `Unified_AMPA_SpaceSR.m`, modify the `filePattern`:

```matlab
% For free glutamate molecules
filePattern = fullfile(myFolder, 'DistanceFree*.txt');

% For adhesion-bound molecules  
filePattern = fullfile(myFolder, 'DistanceBound*.txt');
```

## 📈 Reproducing Manuscript Results

1. **Run diffusion simulations** with N = 1,000 and N = 5,000 particles:
   ```matlab
   DiffusionGlutamateBalls.m
   ```

2. **Simulate receptor activation**:
   ```matlab
   RunAmpa.m        % For AMPA receptors
   NMDA_SpaceSR.m   % For NMDA receptors
   ```

3. **Analyze and compare** using test data from `DataforAMPAandNMDATesting/` folder

4. **Generate visualizations**:
   ```matlab
   PlotDataControl.m
   ```

## 📁 Repository Structure

```
NatureComm/
├── Demo/
│   ├── DiffusionGlutamateBalls.exe
│   └── MATLAB_Runtime_R2022b_Update_10_win64.zip
├── DataforAMPAandNMDATesting/
├── InputParametersSR.m
├── statisticSR.txt
├── DiffusionGlutamateBalls.m
├── NMDA_SpaceSR.m
├── RunAmpa.m
├── Unified_AMPA_SpaceSR.m
├── AMPA.m / AMPA1.m / AMPA2.m
└── PlotDataControl.m
```

## 📝 License

This MATLAB 2022b was used under the **University College London academic use license**. 

## 📚 Citation

If you use this code in your research, please cite the associated paper submitted to *Nature Communications*.

## 🤝 Contributing

This repository supports a specific research publication. For questions or collaboration inquiries, please contact the corresponding author.

## 📞 Contact

**Corresponding Author**: Leonid Savchenko  
**Institution**: UCL Institute of Neurology  
**Repository**: [https://github.com/RusakovLab/NatureComm](https://github.com/RusakovLab/NatureComm)

---

<div align="center">

**⭐ Star this repository if you find it useful!**

</div>
