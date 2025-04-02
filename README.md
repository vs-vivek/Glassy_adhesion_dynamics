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

## Simulation Instructions:
### 1. All the code required to produce migration trajectories for both glassy model and conventional model are provided in the "Codes and Data" folder. 
### 2. Simply download the matlab scripts titled "glassy_model_xD.m" or "constant_tau_off_model.m" and run it using MATLAB. Instructions to change parameters are provided in the scripts but main parameter to change is "eta" i.e viscosity. 
### 3. Sample code runs have also been provided for reference and the migration trajectories of cells for slow and fast relaxing substrates using both glassy and constant_tau_off model can be compared to see the drastic difference.
