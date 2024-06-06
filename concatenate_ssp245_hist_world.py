# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 11:19:55 2023

Concatenate historical data with ssp245 data 

@author: ip294
"""
### import libraries ###
import iris
import numpy as np

### define models for analysis ###
models = ['CESM2-WACCM', 'CNRM-ESM2-1', 'IPSL-CM6A-LR', 'MPI-ESM1-2-LR', 'UKESM1-0-LL']
### data variable to be concatenated 
var = 'pr'

### read in the ssp245 and historical cubes 
path1 = 'C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/'+var+'/historical/Processed_data/World/'
path2 = 'C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/'+var+'/ssp245/Processed_data/World/'

for m in range (np.size(models)): #loop through models
    model = models[m] #select model 
   
    # define date ranges and variant_ids
    hist_range = '1850-2014'        
    
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
    
    #define full filepaths to the historical and ssp245 data
    fname1 = path1+var+'_'+model+'_historical_'+variant_id+'_'+hist_range+'_World.nc'
    fname2 = path2+var+'_'+model+'_ssp245_'+variant_id+'_'+date_range+'_World.nc'
    
    #load in data as iris cubes 
    hist_cube = iris.load_cube(fname1,var)
    ssp_cube = iris.load_cube(fname2, var)
    
    #edit ssp cube times to be after hist cube
    ssp_cube.remove_coord('time')
    newcoord = iris.coords.DimCoord(np.arange(165, np.size(ssp_cube, 0)+165,step=1), standard_name='time', units=1, long_name='time', var_name='time')
    ssp_cube.add_dim_coord(newcoord,0)

    cubes = iris.cube.CubeList([hist_cube, ssp_cube]) #put cubes into a cubelist
    new_cube = cubes.concatenate_cube() #concatenate cubes
    
    if model == 'IPSL-CM6A-LR':
        new_cube = new_cube[:251] #IPSL has a longer dataset for ssp245 so cut off at length 251
    
    #define filepath for saving concatenated data
    outname = 'C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/'+var+'/historical/concatenated_data/World/'+var+'_'+model+'_hist_ssp245_'+variant_id+'_World.nc'
    #save new cube 
    iris.save(new_cube, outname, unlimited_dimensions=[])
    
    