# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 11:41:58 2019

@author: pdlr201

Script to interpolate world interpolated monthly CMIP6 data for specific region
"""

import iris
import numpy as np

experiment = 'ssp585'

################# Specify variable of interest ######################
var = 'cLand' # or treeFrac

################## Specify path to processed data directory #####################
path = 'C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/'+var+'/'+experiment+'/Processed_data/'

################# Specify model here ################################## 
models = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CESM2-WACCM', 'CMCC-CM2-SR5', 'CNRM-ESM2-1', 'EC-Earth3-Veg', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'TaiESM1'] 

for i in range (np.size(models)):
    model = models[i]

    # Specific entries for the 3 models initially used
    if model == 'CESM2-WACCM':
        date_range = '2015-2300'
        variant_id = 'r1i1p1f1'
        
    if model == 'CNRM-ESM2-1':
        date_range = '2015-2100'
        variant_id = 'r1i1p1f2'
    
    if model == 'IPSL-CM6A-LR':
        date_range = '2015-2100'
        variant_id = 'r1i1p1f1'
    
    if model == 'MPI-ESM1-2-LR':
        date_range = '2015-2100'
        variant_id = 'r1i1p1f1'
    
    if model == 'MPI-ESM1-2-HR':
        date_range = '2015-2100'
        variant_id = 'r1i1p1f1'
    
    if model =='UKESM1-0-LL':
        date_range = '2015-2100'  
        variant_id = 'r1i1p1f2'
    
    if model =='ACCESS-ESM1-5':
        date_range = '2015-2300'  
        variant_id = 'r1i1p1f1'  

    if model == 'BCC-CSM2-MR':
        date_range = '2015-2100'
        variant_id = 'r1i1p1f1'        
    
    if model == 'CanESM5':
        date_range = '2015-2100'
        variant_id = 'r1i1p1f1'
    
    if model == 'CanESM5-1':
        date_range = '2015-2100'
        variant_id = 'r1i1p1f1'
    
    if model == 'CAS-ESM2-0':
        date_range = '2015-2100'
        variant_id = 'r1i1p1f1'
    
    if model == 'CMCC-CM2-SR5':
        date_range = '2015-2100'
        variant_id = 'r1i1p1f1'
    
    if model == 'CMCC-ESM2':
        date_range = '2015-2100'
        variant_id = 'r1i1p1f1'
    
    if model == 'CNRM-CM6-1':
        date_range = '2015-2100'
        variant_id = 'r1i1p1f2'
    
    if model == 'EC-Earth3-Veg-LR':
        date_range = '2015-2100'
        variant_id = 'r1i1p1f1'
    
    if model == 'EC-Earth3-Veg':
        date_range = '2015-2100'
        variant_id = 'r1i1p1f1'
        
    if model == 'GISS-E2-1-H':
        date_range = '2015-2100'
        variant_id = 'r1i1p1f2'
    
    if model == 'MIROC-ES2L':
        date_range = '2015-2100'
        variant_id = 'r1i1p1f2'
    
    if model == 'TaiESM1':
        date_range = '2015-2100'
        variant_id = 'r1i1p1f1'
    
    
    # Define latitude and longitude (uses 1 degree by 1 degree grid and these are
    # the coordinates I currently use for the Amazon) 
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
