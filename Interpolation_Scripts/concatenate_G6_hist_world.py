# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 11:19:55 2023

Concatenate historical data with G6 data 

@author: ip294
"""
### import libraries ###
import iris
import numpy as np

### define models for analysis ###
models = ['CESM2-WACCM', 'CNRM-ESM2-1', 'IPSL-CM6A-LR', 'MPI-ESM1-2-LR', 'UKESM1-0-LL']
### data variable to be concatenated ###
var = 'co2'

g6_exp = 'G6solar' #define G6 experiment

### read in the ssp245, historical and G6 cubes 
path1 = 'C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/'+var+'/historical/Processed_data/World/'
path2 = 'C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/'+var+'/ssp245/Processed_data/World/'
path3 = 'C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/'+var+'/'+g6_exp+'/Processed_data/World/'

for m in range (np.size(models)): #loop through models 
    model = models[m] #select models 
   
    ### define date_ranges and variant_id ###
    hist_range = '1850-2014'        
    
    if model == 'CESM2-WACCM':
        ssp_date_range = '2015-2100'
        if g6_exp == 'G6solar':
            g6_date_range = '2019-2100'
            variant_id2 = 'r1i1p1f1'
        elif g6_exp == 'G6sulfur':
            g6_date_range = '2020-2100'
            variant_id2 = 'r1i1p1f2'
        variant_id = 'r1i1p1f1'
                
    if model == 'CNRM-ESM2-1':
        ssp_date_range = '2015-2100'
        variant_id = 'r1i1p1f2'
        variant_id2 = 'r1i1p1f2'
        g6_date_range = '2015-2100'
    
    if model == 'IPSL-CM6A-LR':
        ssp_date_range = '2015-2100'
        g6_date_range = '2020-2100'
        variant_id = 'r1i1p1f1'
        variant_id2 = 'r1i1p1f1'
    
    if model == 'MPI-ESM1-2-LR':
        ssp_date_range = '2015-2100'
        g6_date_range = '2015-2099'
        variant_id = 'r1i1p1f1'
        variant_id2 = 'r1i1p1f1'
    
    if model =='UKESM1-0-LL':
        ssp_date_range = '2015-2100'  
        g6_date_range = '2020-2100'
        variant_id = 'r1i1p1f2'  
        variant_id2 = 'r1i1p1f2'
    
    #define full filepaths to the historical and ssp245 data
    fname1 = path1+var+'_'+model+'_historical_'+variant_id+'_'+hist_range+'_World.nc'
    fname2 = path2+var+'_'+model+'_ssp245_'+variant_id+'_'+ssp_date_range+'_World.nc'
    fname3 = path3+var+'_'+model+'_'+g6_exp+'_'+variant_id2+'_'+g6_date_range+'_World.nc'

    #load in data as iris cubes 
    hist_cube = iris.load_cube(fname1,var)
    ssp_cube = iris.load_cube(fname2, var)
    g6_cube = iris.load_cube(fname3, var)
    
    #edit ssp cube times to be after hist cube according to differences in models 
    if model == 'CNRM-ESM2-1' or model == 'MPI-ESM1-2-LR':
        #edit g6 cube times to be after hist cube
        g6_cube.remove_coord('time')
        newcoord = iris.coords.DimCoord(np.arange(165, np.size(g6_cube, 0)+165,step=1), standard_name='time', units=1, long_name='time', var_name='time')
        g6_cube.add_dim_coord(newcoord,0)
       
        cubes = iris.cube.CubeList([hist_cube, g6_cube]) #put cubes into a cubelist
        new_cube = cubes.concatenate_cube() #concatenate cubes
       
    elif model == 'IPSL-CM6A-LR' or model == 'UKESM1-0-LL':
        ssp_cube.remove_coord('time')
        newcoord = iris.coords.DimCoord(np.arange(165, np.size(ssp_cube, 0)+165,step=1), standard_name='time', units=1, long_name='time', var_name='time')
        ssp_cube.add_dim_coord(newcoord,0)
        ssp_cube = ssp_cube[:5,:,:]
        
        g6_cube.remove_coord('time')
        newcoord = iris.coords.DimCoord(np.arange(170, np.size(g6_cube, 0)+170,step=1), standard_name='time', units=1, long_name='time', var_name='time')
        g6_cube.add_dim_coord(newcoord,0)
      
        cubes = iris.cube.CubeList([hist_cube, ssp_cube, g6_cube]) #put cubes into a cubelist
        new_cube = cubes.concatenate_cube() #concatenate cubes
    
    elif model == 'CESM2-WACCM' and g6_exp == 'G6solar':
        ssp_cube.remove_coord('time')
        newcoord = iris.coords.DimCoord(np.arange(165, np.size(ssp_cube, 0)+165,step=1), standard_name='time', units=1, long_name='time', var_name='time')
        ssp_cube.add_dim_coord(newcoord,0)
        ssp_cube = ssp_cube[:4,:,:]
        
        g6_cube.remove_coord('time')
        newcoord = iris.coords.DimCoord(np.arange(169, np.size(g6_cube, 0)+169,step=1), standard_name='time', units=1, long_name='time', var_name='time')
        g6_cube.add_dim_coord(newcoord,0)
      
        cubes = iris.cube.CubeList([hist_cube, ssp_cube, g6_cube]) #put cubes into a cubelist
        new_cube = cubes.concatenate_cube() #concatenate cubes
    
    elif model == 'CESM2-WACCM' and g6_exp == 'G6sulfur':
        ssp_cube.remove_coord('time')
        newcoord = iris.coords.DimCoord(np.arange(165, np.size(ssp_cube, 0)+165,step=1), standard_name='time', units=1, long_name='time', var_name='time')
        ssp_cube.add_dim_coord(newcoord,0)
        ssp_cube = ssp_cube[:5,:,:]
        
        #edit ssp cube times to be after hist cube
        g6_cube.remove_coord('time')
        newcoord = iris.coords.DimCoord(np.arange(170, np.size(g6_cube, 0)+170,step=1), standard_name='time', units=1, long_name='time', var_name='time')
        g6_cube.add_dim_coord(newcoord,0)
      
        cubes = iris.cube.CubeList([hist_cube, ssp_cube, g6_cube]) #put cubes into a cubelist
        new_cube = cubes.concatenate_cube() #concatenate cubes
    
    #define filepath for saving concatenated data
    outname = 'C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/'+var+'/historical/concatenated_data/World/'+var+'_'+model+'_hist_'+g6_exp+'_'+variant_id2+'_World.nc'
    #save new cube
    iris.save(new_cube, outname, unlimited_dimensions=[])
    
    
