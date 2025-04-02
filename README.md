# Glassy Adhesion Dynamics Based Cell Migration

#### Contained within this repository are codes and data for generating migration trajectories (for both conventional model and glassy model) on varying levels of viscosity and elasticity of the substrate and their analysis. 
### This repository relates to the research publication "Sharma V. et al., Glassy adhesion dynamics govern transitions between sub-diffusive and super-diffusive cell migration on viscoelastic substrates, BioRxiv, 2025" (https: ). Please cite our paper if this code is used or modified for further use.
## System Requirements
### Operating system:
#### This package is supported for macOS and windows. The package has been tested on the following systems:
#### •	macOS Monterey 12.4, Processor: 2.3 GHz 8-Core Intel i9, RAM: 32 GB
#### •	Windows 11 Home, Processor: AMD Ryzen 5 5500U Hexa-Core processor, RAM: 16 GB
### MATLAB dependencies
#### This package requires MATLAB and has been tested on the following software versions:
#### •	MATLAB_R2024
## Installation guide
### 1.	Install MATLAB (https://www.mathworks.com/products/matlab.html). The installation usually takes approximately 1 hour.

## Simulation Instructions

### 1. Code Overview  
All MATLAB scripts required to generate migration trajectories for both the **glassy model** and the **conventional (constant τ<sub>off</sub>) model** are provided in the **`Codes and Data`** folder.

### 2. Running the Simulations  
Download and run either `glassy_model_xD.m` or `constant_tau_off_model.m` using MATLAB. Instructions for modifying parameters are provided within the scripts—most notably the parameter `eta`, which controls the substrate viscosity.

### 3. Sample Outputs  
Sample outputs are included for reference. You can compare migration trajectories of cells on **slow-** and **fast-relaxing substrates** using both models.

> ⚠️ Note: Since the simulations are stochastic, individual runs may vary slightly, but the key trends will remain consistent.  

**Example behaviors:**  
- In the **glassy model**, cells on slow-relaxing substrates tend to become trapped, while on fast-relaxing substrates, they escape trapping and take longer steps.  
- In the **constant τ<sub>off</sub> model**, cells display uniform behavior across substrate types, with no trapping observed.

### 4. Data Analysis  
Various quantities of interest (e.g., displacement, speed, persistence) can be extracted from the simulation outputs for further analysis, as described in the manuscript and supplementary information.

