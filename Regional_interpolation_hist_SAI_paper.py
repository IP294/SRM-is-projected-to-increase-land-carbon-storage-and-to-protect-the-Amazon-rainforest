# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 11:41:58 2019

@author: pdlr201

Script to interpolate world interpolated monthly CMIP6 data for specific region
"""
### import libraries ###
import iris
import numpy as np

#### Specify area of interest ####
experiment = 'historical'

######### Specify variable of interest #####
var = 'cLand' # or treeFrac

################## Specify path to processed data directory #####################
path = 'C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/'+var+'/'+experiment+'/Processed_data/'

################# Specify model here ################################## 
models = ['CESM2-WACCM', 'CNRM-ESM2-1', 'IPSL-CM6A-LR', 'MPI-ESM1-2-LR', 'MPI-ESM1-2-HR', 'UKESM1-0-LL']

for m in range (np.size(models)): #loop through models 
    model = models[m] #select model

    ##### define details associated with each model #####
    if model == 'EC-Earth3-Veg':
        date_range = '1850-2000' #date range in file name
        variant_id = 'r1i1p1f1' #variant ID
        
    elif model == 'CNRM-ESM2-1':
        date_range = '1850-2014'
        variant_id = 'r1i1p1f2'
        
    elif model == 'MPI-ESM1-2-LR':
        date_range = '1850-2014'
        variant_id = 'r1i1p1f1'
        
    elif model == 'MPI-ESM1-2-HR':
        date_range = '1850-2014'
        variant_id = 'r1i1p1f1'
        
    elif model == 'UKESM1-0-LL':
        date_range = '1850-2014'
        variant_id = 'r1i1p1f2'
        
    elif model == 'CESM2-WACCM':
        date_range = '1850-2014'
        variant_id = 'r1i1p1f1'
        
    elif model == 'IPSL-CM6A-LR':
        date_range = '1850-2014'
        variant_id = 'r1i1p1f1'
        
    elif model == 'GFDL-ESM4':
        date_range = '1850-1999'
        variant_id = 'r1i1p1f1'
        
    elif model == 'NorCPM1':
        date_range = '1850-2013'
        variant_id = 'r1i1p1f1'
        
    elif model == 'TaiESM1':
        date_range = '1850-2000'
        variant_id = 'r1i1p1f1'
        
    elif model == 'SAM0-UNICON':
        date_range = '1850-1999'
        variant_id = 'r1i1p1f1'
        
    elif model == 'CanESM5':
        date_range = '0000-0000'
        variant_id = 'r1i1p1f1'
        
    elif model == 'BCC-ESM1':
        date_range = '1850-2300'
        variant_id = 'r1i1p1f1'
        
    elif model == 'CESM2':
        date_range = '0000-0000'
        variant_id = 'r1i1p1f1'
        
    elif model == 'AWI-ESM-1-1-LR':
        date_range = '1855-1954'
        variant_id = 'r1i1p1f1'
        
    elif model == 'ACCESS-ESM1-5':
        date_range = '100-1100'
        variant_id = 'r1i1p1f1'
        
    elif model == 'BCC-CSM2-MR':
        date_range = '1850-2449'
        variant_id = 'r1i1p1f1'
        
    elif model == 'CMCC-CM2-SR5':
        date_range = '1850-2349'
        variant_id = 'r1i1p1f1'
        
    elif model == 'CMCC-ESM2':
        date_range = '1850-2349'
        variant_id = 'r1i1p1f1'
        
    elif model == 'MPI-ESM-1-2-HAM':
        date_range = '1850-2849'
        variant_id = 'r1i1p1f1'
        
    elif model == 'NorESM2-LM':
        date_range = '1600-2100'
        variant_id = 'r1i1p1f1'
        
    elif model == 'MIROC-ES2L':
        date_range = '1850-2350'
        variant_id = 'r1i1p1f1'
        
    # Define latitude and longitude (uses 1 degree by 1 degree grid) of region
    my_lats=['latitude',np.arange(-20,14,step=1)]
    my_lons=['longitude',np.arange(-83,-34,step=1)]
    
    ############ Specify name of directory interpolated World data is stored in #################
    region1 = 'World'
    
    ############ Specify name of directory interpolated Amazon data will be stored in #################
    region2 = 'Amazon'
    
    # File name of interpolated World data
    fname = path+region1+'/'+var+'_'+model+'_'+experiment+'_'+variant_id+'_'+date_range+'_'+region1+'.nc'
    
    # File name for new interpolated S. America data
    outname = path+region2+'/'+var+'_'+model+'_'+experiment+'_'+variant_id+'_'+date_range+'_'+region2+'.nc'
    
    # Load in world data into iris cubeme,var)
    world_cube = iris.load_cube(fname,var)
    
    # Interpolate cube onto S. America grid
    region_cube = world_cube.interpolate([my_lats, my_lons],iris.analysis.Linear(extrapolation_mode='nan'))
    
    # Ensure new cube has correct name
    region_cube.rename(var)
    
    # Save data to file name
    iris.save(region_cube, outname, unlimited_dimensions=[])
