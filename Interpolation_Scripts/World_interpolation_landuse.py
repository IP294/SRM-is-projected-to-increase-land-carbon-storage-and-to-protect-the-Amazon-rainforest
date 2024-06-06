#- coding: utf-8 -*-
"""
Created on Mon Dec 16 09:01:27 2019

@author: pdlr201

Script to interpolate annual 1pctCO2 CMIP6 data onto 1 degree world grid
"""

import numpy as np
from netCDF4 import Dataset
import iris

def shiftlon(lon,lon_0):
    """returns original sequence of longitudes (in degrees) recentered
    in the interval [lon_0-180,lon_0+180]"""
    lon_shift = np.asarray(lon)
    lon_shift = np.where(lon_shift > lon_0+180, lon_shift-360 ,lon_shift)
    lon_shift = np.where(lon_shift < lon_0-180, lon_shift+360 ,lon_shift)
    itemindex = len(lon)-np.where(lon_shift[0:-1]-lon_shift[1:] >= 180)[0]
    return (np.roll(lon_shift,itemindex-1), itemindex-1)


######################## Specify path to data ##########################
# (I have 2 directories a raw data directory and a processed data directory)
# path = 'C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/land_use/Original_data/'   
# ######################## Specify directory of raw data ##################
# raw_data = '/Original_data/'

######################## Specify directory of processed data ##################
# processed_data = '/Processed_data/'
variable = 'primf'
var_long_name = 'potentially forested secondary land'

# variable2 = ''se

experiment = 'ssp585'

############ Specify directory within processed data directory for the region interpolated #################
region = 'World'

# Define latitude and longitude (uses 1 degree by 1 degree grid) 
my_lats=['latitude',np.arange(-90,91,step=1)]
my_lons=['longitude',np.arange(-180,181,step=1)]


# Specify path to data
path = 'C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/land_use/Original_data/'

# File name of data
fname = path + 'multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-'+experiment+'-2-1-f_gn_2015-2100.nc'
fname2 = path + 'multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-'+experiment+'-2-1-f_gn_2015-2100.nc'

# Check if file exists
import os
if not os.path.isfile(fname):
    raise FileNotFoundError(f"File {fname} not found.")

total_years = 86 

# Load NetCDF file
try:
    with Dataset(fname, 'r') as nc:
        # Load data and coordinates
        data = nc.variables[variable][:]
        lat = nc.variables['lat'][:]
        lon = nc.variables['lon'][:]

    # Create Iris cube
    cube = iris.cube.Cube(data, long_name=var_long_name, units='1')

    # Add latitude and longitude coordinates
    time_coord = iris.coords.DimCoord(range(total_years), standard_name='time', long_name='time', units=1)
    
    lon_coord = iris.coords.DimCoord(lon, standard_name='longitude', units='degrees')
    lat_coord = iris.coords.DimCoord(lat, standard_name='latitude', units='degrees')
    cube.add_dim_coord(lat_coord, 1)
    cube.add_dim_coord(lon_coord, 2)
    cube.add_dim_coord(time_coord, 0)
    
    # Print information about the cube
    print("Cube loaded successfully:")
    print(cube)

except Exception as e:
    print("Error loading cube:", e)

# Centre coordinates about longitude 0
cube.coord('longitude').points, rollindex = shiftlon(cube.coord('longitude').points, 0)

# Roll data to be consistent with longitude 0
cube.data = np.roll(cube.data, rollindex, axis=2)

# Interpolate data onto new coordinate system
region_cube = cube.interpolate([my_lats, my_lons],iris.analysis.Linear(extrapolation_mode='nan'))

print(region_cube)

region_cube.rename(variable)

# New file name for data
path2 = 'C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/land_use/Processed_data/'
outname = path2+region+'/'+variable+'_'+experiment+'_states_landuse.nc'

# Save data to file name
iris.save(region_cube, outname, unlimited_dimensions=[])
