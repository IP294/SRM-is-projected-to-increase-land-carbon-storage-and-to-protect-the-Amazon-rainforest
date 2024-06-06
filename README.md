# SRM-is-projected-to-increase-land-carbon-storage-and-to-protect-the-Amazon-rainforest

This repository contains the necessary Python code to analyse the data and produce the figures as in the manuscript "Solar Radiation Modification is projected to increase land carbon storage and to protect the Amazon rainforest". 

CMIP6 data for use with this code is publicly available at: https://doi.org/10.22033/ESGF/CMIP6.10034 (Danabasoglu, 2019), https://doi.org/10.22033/ESGF/CMIP6.3907 (Séférian, 2019), https://doi.org/10.22033/ESGF/CMIP6.5059 (Boucher et al., 2020), https://doi.org/10.22033/ESGF/CMIP6.6448 (Niemeier et al., 2019), https://doi.org/10.22033/ESGF/CMIP6.5822 (Jones, 2019), https://doi.org/10.22033/ESGF/input4MIPs.10468 (Hurtt et al, 2019). 

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
All scripts are expected to run within few minutes or less on a normal desktop computor. 

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

## Figure scripts
The following scripts are dependent on the above interpolation and concatenation scripts being run in order for the neccessary files to be available for analysis. 

- Time_series_plot_geo_hist_World_figure1_final.py - A script which uses  interpolated data to produce Figure 1 of the manuscipt (Timeseries showing the evolution of the decadal means of surface temperature, precipitation, net primary productivity, and land carbon anomalies relative to the pre-industrial period with time from 1900 to 2100.)
- Maps_Geo_Figures2&3_final.py - A script which uses interpolated data to produce Figures 2 & 3 of the manuscript (Maps showing the difference between the 2090-2100 means of G6sulfur and SSP585/SSP245 for surface temperature (oC), precipitation (mm day-1), net primary productivity (kgC m-2 yr-1), and land carbon (kgC m-2).)
- Time_series_plot_geo_World_figure4_final.py - A script which uses interpolated data to produce Figure 4 of the manuscript (evolution of the global mean changes in NPP and land carbon anomalies relative to the pre-industrial period against temperature anomaly and CO2.)
-Time_series_plot_geo_World_SM1.py - A script which uses interpolated data to produce Supplementary figure 1 (A timeseries showing the evolution of land carbon relative to the pre-industrial period for all models analyses in G6sulfur, SSP585, and SSP245.)
- Maps_Geo_comp_Amazon.py - A script which uses interpolated data to produce Supplementary figures 2&4 (maps showing the difference between the 2090-2100 means of G6sulfur and SSP585/SSP245 in the Amazon for surface temperature (oC), precipitation (mm day-1), net primary productivity (kgC m-2 yr-1, and land carbon (kgC m-2).)
- Plotting_changes_World_landuse_SM3.py - A script which uses interpolated land use data to produce Supplementary figure 3 (map showing the difference between prescribed forested land (primf+secdf) in SSP245 versus SSP585.)
- Time_series_plot_geo_Amazon_SM5.py - A script which uses interpolated data to produce Supplementary figure 5 (timeseries showing the evolution of the decadal mean of surface temperature, precipitation, net primary productivity, and land carbon anomalies relative to the pre-industrial period in the Amazon with time from 1900 to 2100.) 

