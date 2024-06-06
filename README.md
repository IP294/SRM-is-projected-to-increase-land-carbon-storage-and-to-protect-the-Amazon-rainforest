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
The following scripts are designed to interpolate single or multiple files of CMIP6 data onto a common world grid with 1 degree resolution, producing a single output file of interpolated data which may then be used for either further interpolation/concatenation or analysis.
  - World_interpolation_Geo.py - designed to run on G6sulfur datasets
  - World_interpolation_hist.py - designed to run on historical datasets
  - World_interpolation_landuse.py - designed to run on land use change data
  - World_interpolation_ssp245.py - designed to run on ssp245 datasets
  - World_interpolation_ssp585.py - designed to run on ssp585 datasets
  - World_interpolation_co2.py - designed to run on co2 datasets

The following scripts interpolate the data from files output by the World_interpolation scripts onto a regional grid of the Amazon with 1 degree resolution. A single file of regionally interpolated data is output which may then be used in either concatenation scripts or analysis. 
  - Regional_interpolation_Geo.py - designed to run on G6sulfur datasets
  - Regional_interpolation_hist.py - designed to run on historical datasets
  - Regional_interpolation_landuse.py - designed to run on land use change data
  - Regional_interpolation_co2.py - designed to run on co2 datasets
  - Regional_interpolation_ssp585.py - designed to run on ssp585 datasets
  - Regional_interpolation_ssp245.py - designed to run on ssp245 datasets

The following scripts concatenate historical datasets with ssp245, ssp585 and G6sulfur datasets to create a continuous run of data from 1850-2100 for each scenario. 
  - concatenate_G6_world.py - concatenates worldwide historical datasets with G6sulfur datasets, for some models ssp585 data is used to fill the gap between the end of the       historical run and the beginning of G6sulfur runs
  - concatenate_ssp245_hist_world.py - concatenates worldwide historical datasets with ssp245 datasets
  - concatenate_ssp585_hist_world.py - concatenates worldwide histotical datasets with ssp585 datasets 
  - concatenate_ssp585_hist_region.py - concatenates regional (Amazon) historical datasets with ssp585 or ssp245 datasets 
  
  

  

