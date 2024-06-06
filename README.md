# SRM-is-projected-to-increase-land-carbon-storage-and-to-protect-the-Amazon-rainforest

This repository contains the necessary Python code to analyse the data and produce the figures as in the manuscript "Solar Radiation Modification is projected to increase land carbon storage and to protect the Amazon rainforest". 

CMIP6 data for use with this code is publicly available at: https://doi.org/10.22033/ESGF/CMIP6.10034 (Danabasoglu, 2019), https://doi.org/10.22033/ESGF/CMIP6.3907 (Séférian, 2019), https://doi.org/10.22033/ESGF/CMIP6.5059 (Boucher et al., 2020), https://doi.org/10.22033/ESGF/CMIP6.6448 (Niemeier et al., 2019), https://doi.org/10.22033/ESGF/CMIP6.5822 (Jones, 2019). 

## Repository Contents
- Interpolation scripts: Python code for interpolating CMIP6 data
- Figure scripts: Python code for creating figures, as found in the manuscript "Solar Radiation Modification is projected to increase land carbon storage and to protect the Amazon rainforest"

# System Requirements 
## Hardware Requirements
The code in this repository requires only a standard computor with enough RAM to support the running of the code. The code was developed on a computor with the following specs:
  RAM: 16 GB
  CPU: 4 cores, 2.4 GHz/core 

## Software Requirements
### OS Requirements 
The code included in this repository is tested on Windows operating systems and has been tested on Windows 10. It should however be compatible with Windows, Mac and Linux operating systems. 

Users should have Python installed before trying to run any code. The code was developed using Python version 3.10. 

# Installation Guide
## Development version
### Package dependencies
Users should install the following packages prior to running any of the code in this repository: 
  - numpy V_1.23.5
  - iris V_3.6.1
  - matplotlib V_3.7.0
  - netCDF4 V_1.6.4
  - cartopy V_0.21.1
  - pathlib V_1.0.1

# Demo
## Interpolation Scripts
These scripts are designed to interpolate the CMIP6 data from different models onto a common grid with 1 degree resolution. 
  

