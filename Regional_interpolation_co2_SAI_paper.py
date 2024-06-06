o# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 11:41:58 2019

@author: pdlr201

Script to interpolate world interpolated co2 CMIP6 data for specific region
"""
### import libraries ###
import iris
import numpy as np

#### specify experiment of interest ####
experiment = 'ssp585'

#### Specify variable of interest #####
var = 'co2' 

################## Specify path to processed data directory #####################
path = 'C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/'+var+'/'+experiment+'/Processed_data/'

################# Specify model here ################################## 
models = ['ACCESS-ESM1-5', 'CanESM5', 'CESM2', 'MIROC-ES2L', 'NorESM2-LM', 'UKESM1-0-LL']

for m in range(np.size(models)): #loop through models
    model=models[m] #select model 
    
    ##### define details associated with each model #####
    if model == 'NorESM2-LM':
        date_range = '141-399' #date range of file
        variant_id = 'r1i1p1f1' #varaint ID
        
    elif model == 'UKESM1-0-LL':
        date_range = '1999-2639'
        variant_id = 'r1i1p1f2'
        
    elif model == 'CanESM5':
        date_range = '1991-2290'
        variant_id = 'r1i1p2f1'
    
    elif model == 'ACCESS-ESM1-5':
        date_range = '241-1000'
        variant_id = 'r1i1p1f1'
    
    elif model == 'MIROC-ES2L':
        date_range = '1989-2388'
        variant_id = 'r1i1p1f2'
    
    elif model == 'CESM2':
        date_range = '1-199'
        variant_id = 'r1i1p1f1'
    
    elif model == 'CNRM-ESM2-1':
        date_range = '1990-2190'
        variant_id = 'r1i1p1f2'
    
    elif model == 'GFDL-ESM4':
        date_range = '140-340'
        # date_range2 = '01-150'
        variant_id = 'r1i1p1f1'
    
    elif model == 'BCC-CSM2-MR':
        date_range = '2015-2100'
        variant_id = 'r1i1p1f1'
    
    # Define latitude and longitude (uses 1 degree by 1 degree grid) for Amazon 
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
    
    # Load in world data into iris cube
    world_cube = iris.load_cube(fname,var)
    
    # Interpolate cube onto S. America grid
    region_cube = world_cube.interpolate([my_lats, my_lons],iris.analysis.Linear(extrapolation_mode='nan'))
    
    # Ensure new cube has correct name
    region_cube.rename(var)
    
    # Save data to file name
    iris.save(region_cube, outname, unlimited_dimensions=[])
