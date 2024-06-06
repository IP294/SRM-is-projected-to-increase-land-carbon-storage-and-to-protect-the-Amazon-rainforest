# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 11:41:58 2019

@author: pdlr201

Script to interpolate world interpolated monthly CMIP6 data for specific region
"""

import iris
import numpy as np

experiment = 'ssp245'

################# Specify variable of interest ######################
var = 'cLand' # or treeFrac

################## Specify path to processed data directory #####################
path = 'C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/'+var+'/'+experiment+'/Processed_data/'

################# Specify model here ################################## 
models = ['CNRM-ESM2-1', 'IPSL-CM6A-LR', 'MPI-ESM1-2-LR', 'UKESM1-0-LL', 'CESM2-WACCM']
for i in range (np.size(models)):
    model = models[i]

    # Specific entries for the 3 models initially used
    if model == 'CESM2-WACCM':
        date_range = '2015-2100'
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
