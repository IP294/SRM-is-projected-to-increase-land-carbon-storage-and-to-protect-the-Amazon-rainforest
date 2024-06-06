# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 11:19:55 2023

Concatenate historical data with ssp585 data 

@author: ip294
"""

### import libraries ###
import iris
import numpy as np

### define models for analysis ###
models = ['CESM2-WACCM', 'CNRM-ESM2-1', 'IPSL-CM6A-LR', 'MPI-ESM1-2-LR', 'UKESM1-0-LL']
### variable for concatentation ###
var = 'pr' 

### read in the ssp585 and historical cubes 
path1 = 'C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/'+var+'/historical/Processed_data/World/'
path2 = 'C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/'+var+'/ssp585/Processed_data/World/'

for m in range (np.size(models)): #loop through models
    model = models[m] #select models 
   
    ### define date ranges and variant ids 
    hist_range = '1850-2014'        
    
    if model == 'CESM2-WACCM':
        date_range = '2015-2300'
        variant_id = 'r1i1p1f1'
        
    if model == 'CNRM-ESM2-1':
        date_range = '2015-2100'
        variant_id = 'r1i1p1f2'
    
    if model == 'IPSL-CM6A-LR':
        date_range = '2015-2300'
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
    
    #define full file names 
    fname1 = path1+var+'_'+model+'_historical_'+variant_id+'_'+hist_range+'_World.nc'
    fname2 = path2+var+'_'+model+'_ssp585_'+variant_id+'_'+date_range+'_World.nc'
    
    #load in data as iris cubes
    hist_cube = iris.load_cube(fname1,var)
    ssp_cube = iris.load_cube(fname2, var)
    
    #edit ssp cube times to be after hist cube
    ssp_cube.remove_coord('time')
    newcoord = iris.coords.DimCoord(np.arange(165, np.size(ssp_cube, 0)+165,step=1), standard_name='time', units=1, long_name='time', var_name='time')
    ssp_cube.add_dim_coord(newcoord,0)

    if model == 'CESM2-WACCM' or model == 'IPSL-CM6A-LR': #these models have longer datasets, so cut off at length 86 
        ssp_cube = ssp_cube[:86,:,:]
        
    cubes = iris.cube.CubeList([hist_cube, ssp_cube]) #put cubes ito cubelist
    new_cube = cubes.concatenate_cube() #concatenate cubes
    
    #define filepath for saving concatenated data
    outname = 'C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/'+var+'/historical/concatenated_data/World/'+var+'_'+model+'_hist_ssp585_'+variant_id+'_World.nc'
    #save new cube 
    iris.save(new_cube, outname, unlimited_dimensions=[])
    
    
