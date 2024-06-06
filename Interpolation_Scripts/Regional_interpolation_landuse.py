# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 11:41:58 2019

@author: pdlr201

Script to interpolate world interpolated monthly CMIP6 data for specific region
"""
### import libraries ###
import iris
import numpy as np

##### Specify experiment of interest ######
experiment = 'ssp245'

###### Specify variable of interest #########
var = 'primf'
var_long_name = 'potentially forested secondary land' #longform name of variable 

################## Specify path to processed data directory #####################
path = 'C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/land_use/Processed_data/'


# Define latitude and longitude (uses 1 degree by 1 degree grid) for region
my_lats=['latitude',np.arange(-20,14,step=1)]
my_lons=['longitude',np.arange(-83,-34,step=1)]

############ Specify name of directory interpolated World data is stored in #################
region1 = 'World'

############ Specify name of directory interpolated Amazon data will be stored in #################
region2 = 'Amazon'

# File name of interpolated World data
fname = path +region1+'/'+var+'_'+experiment+'_states_landuse.nc'

# File name for new interpolated S. America data
outname = path+region2+'/'+var+'_'+experiment+'_states_landuse.nc'

# Load in world data into iris cube
world_cube = iris.load_cube(fname,var)

# Interpolate cube onto S. America grid
region_cube = world_cube.interpolate([my_lats, my_lons],iris.analysis.Linear(extrapolation_mode='nan'))

# Ensure new cube has correct name
region_cube.rename(var)

# Save data to file name
iris.save(region_cube, outname, unlimited_dimensions=[])
