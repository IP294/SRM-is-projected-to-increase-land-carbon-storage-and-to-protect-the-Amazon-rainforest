# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 12:55:01 2021

@author: impy2

Timeseries showing the evolution of the decadal mean of (a) surface temperature,
(b) precipitation, (c) net primary productivity, and (d) land carbon anomalies
relative to the pre-industrial period in the Amazon.

"""
### Import libraries ###
from pathlib import Path
import iris.coords
import iris.cube
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

################# Specify the experiment and models ###############
experiments = ['G6sulfur', 'ssp585', 'ssp245']
colours = ['royalblue', 'firebrick','goldenrod']#, 'forestgreen', 'indigo'] # specify colours for each model in plot 
models = ['CESM2-WACCM', 'CNRM-ESM2-1', 'IPSL-CM6A-LR', 'MPI-ESM1-2-LR', 'UKESM1-0-LL']
letters = ['a','b', 'c', 'd'] #specify letter labels for panels of figure     

################# Specify variable of interest ######################
variables = ['tas', 'pr', 'NPP', 'cLand'] #variables to be plotted
variables2 = ['tas', 'pr', 'npp','cland'] #lower case version of vars for file path

names = ['Near Surface Air Temperature', 'Precipitation Rates', 'Net Primary Productivity','Land Carbon'] #longform variable names 

### initialize ensemble arrays ###
ens = np.zeros((np.size(experiments),np.size(models), 241))
temp_ens = np.zeros((np.size(experiments),np.size(models), 241))
a = np.zeros((5,34,49))

### initialize figure ###
fig, axes = plt.subplots(2,2, figsize=(15, 9))
axes = axes.flatten()

for v in range (np.size(variables)): #loop through variables 
    var = variables[v] #select variable 
    var2 = variables2[v]
    variable_name = names[v] #select variable name 
    
    ############## Specify units of variable ##################
    if var == 'cVeg':
        units = 'kgC m$^{-2}$'
    elif var == 'tas':
        units = '$^\circ$C'
    elif var == 'pr':
        units = 'mm day$^{-1}$'
    elif var == 'NPP':
        units = 'PgC yr$^{-1}$'
    elif var == 'cSoil' or var == 'cLand':
        units = 'PgC'
        
        
    axes[v].set_xlabel('Year', fontsize=13) # set x axis label 
    axes[v].set_ylabel(variable_name+' ('+units+')') # set y axis label
    axes[v].set_title('Change in '+ variable_name) # set panel title
    axes[v].text(.02, 1.07, letters[v], ha='left', va='top', transform=axes[v].transAxes, fontsize='large') #label the panel 
    
    for i in range (np.size(models)): # loop through all models 
        model = models[i] # select model
        
        for k in range (np.size(experiments)): # loop through each experiment
            experiment = experiments[k] # select experiment 
            colour = colours[k] # select colour for the experiment line

            #### select date_range and variant_id corresponding to model ####   
            if model == 'CESM2-WACCM':
                if experiment == 'G6solar':
                    variant_id = 'r1i1p1f1'
                elif experiment == 'G6sulfur':
                    variant_id = 'r1i1p1f2'
                elif experiment == 'ssp585':
                    variant_id = 'r1i1p1f1'
                elif experiment == 'ssp245':
                    variant_id = 'r1i1p1f1'
                hist_dates = '1850-2014'
                hist_id = 'r1i1p1f1'
                
            elif model == 'CNRM-ESM2-1':
                variant_id = 'r1i1p1f2'
                hist_dates = '1850-2350'
                hist_id = 'r1i1p1f2'
        
            elif model == 'IPSL-CM6A-LR':
                variant_id = 'r1i1p1f1'
                hist_dates = '1850-2014'
                hist_id = 'r1i1p1f1'
        
            elif model == 'MPI-ESM1-2-LR':
                variant_id = 'r1i1p1f1'
                hist_dates = '1850-2014'
                hist_id = 'r1i1p1f1'
        
            elif model =='UKESM1-0-LL':
                variant_id = 'r1i1p1f2' 
                hist_dates = '1850-2014'
                hist_id = 'r1i1p1f2'

            ############ Specify name of directory of data which you want to plot #################
            region = 'Amazon'

            ################## Specify path to processed data directory ###########
            path = 'C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/'+var+'/historical/concatenated_data/'+region+'/'
            
            # File name of interpolated data
            fname = path+var+'_'+model+'_hist_'+experiment+'_'+variant_id+'_'+region+'.nc'
            fname = Path(fname)

            # Load in data
            try:
                f = nc.Dataset(fname,'r')
                
            except FileNotFoundError as e: #exception if file fails to load
                print(model+' failed: {}'.format(e))
        
            else:
                x = f.variables[var2][:]
                
                #close dataset
                f.close()
                
                ##### Define latitude and longitude (uses 1 degree by 1 d) ######
                my_lats=np.arange(-20,14,step=1)
                my_lons=np.arange(-83,-34,step=1)
                
                #define dimensions to be used in iris cube manipulations
                ny =  np.size(my_lats)                 
                nx = np.size(my_lons)
                
                #create new iris coordinates for latitude and longitude 
                coord1 = iris.coords.DimCoord(my_lats,bounds=np.array([my_lats-0.5,my_lats+0.5]).T, standard_name='latitude', units='degrees', var_name='lat', attributes={'title': 'Latitude', 'type': 'double'}, circular=True)
                coord2 = iris.coords.DimCoord(my_lons,bounds=np.array([my_lons-0.5,my_lons+0.5]).T, standard_name='longitude', units='degrees', var_name='lon', attributes={'title': 'Longitude', 'type': 'double', 'valid_max': 180.0, 'valid_min': -180.0}, circular=True)
                cube = iris.cube.Cube(np.zeros((ny,nx),np.float32),dim_coords_and_dims=[(coord1,0),(coord2,1)])
            
                #extract longitude and latitude
                latitude = coord1.points
                longitude = coord2.points
            
                # Obtain dimension sizes of the dataset
                nt = int(len(x[:,0,0]))
                ny = int(len(x[0,:,0]))
                nx = int(len(x[0,0,:]))
                
                if var == 'cLand' or var == 'NPP': #load in land mask for cland and NPP data 
                ############ Specify path to land mask ######################
                    fname2 = 'C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/Masks/'+region+'.nc'
                    fname2 = Path(fname2)
                    
                    # Load in land mask and remove any ocean values
                    f2 = nc.Dataset(fname2,'r')
                    sftlf = f2.variables['sftlf'][:]
                    
                    sftlf_mask = np.broadcast_to(sftlf, (nt, ny,nx))
                    x[np.where(sftlf_mask <= 0)] = np.nan  
                    x[np.where(~np.isfinite(sftlf_mask))] = np.nan #remove ocean values 
                    
                    #create areas array for weighting data 
                    areas = iris.analysis.cartography.area_weights(cube, normalize=False)
                    
                    land_area = np.sum(areas[np.where((sftlf>0)&(np.isfinite(sftlf)))])
                    
                    #close land mask data 
                    f2.close()
                    
                #ensure data is all of the same length
                data = np.empty((251, ny, nx)) #create empty dataset of the required dimensions
                data.fill(np.nan) #fill with nans
                data[:np.size(x,axis=0), :,:] = x[:251,:,:] #fill array with data
                data[data > 1e36] = np.nan
                
                #convert units
                if var == 'pr':
                    data = data * 86400
                elif var == 'tas':
                    data = data -273
                elif var == 'NPP':
                    data = data * 3.154e7 #*1e-12 #multiplied to make per year and converted to Pg
    
                t = np.arange(1855, np.size(data,0)-10 +1855) # create time array for x axis

                areas = iris.analysis.cartography.area_weights(cube, normalize=False)
                
                if var == 'cLand' or var == 'NPP':
                    #### sum variable across region ####
                    vsum = np.zeros(np.size(data, 0)) #initialize zero array
                    
                    for j in range (np.size(data, 0)):
                        v_point = data[j,:,:] #select data from single timepoint
                        vsum[j] = np.nansum(v_point*areas)*1e-12 #convert to petagrams
                        
                        if vsum[j] == 0:
                            vsum[j] = np.nan #if global sum is 0, relplace with nan
                    
                    #### calculate the anomaly of the y variable relative to the first 50 years of the hist run(if cveg or cland)
                    v_anom = vsum - np.nanmean(vsum[0:50])

                elif var == 'tas' or var == 'pr': 
                    ## if variable temp or precip find average across region ###
                    vavr = np.nansum((data*areas),axis=1)/np.nansum(areas) # finding the weighted average 
                    vavr = np.nansum(vavr, axis=1)
                    vavr[vavr==0] = np.nan
                        
                    v_anom = vavr - np.nanmean(vavr[0:50]) #find anomaly relative to preindustrial period
                     
            smoothed_anom = np.zeros(len(v_anom)-10) #initialize array for smoothed data 
            
            #find average across 10 year window 
            for j in range (len(v_anom)-10):
                m1 = j #start of window               
                m2 = j+10 #end of window 
                smoothed_anom[j] = np.mean(v_anom[m1:m2]) #find mean across 10 year window 
                
            #### save data into ens for later averaging ####       
            ens[k,i,:] = smoothed_anom   
            
        axes[k].tick_params(labelcolor="black", bottom=True, left=True)
        fig.subplots_adjust(hspace=0.5) #adjust spacing of figures 

    ### calculate ensemble averages ###
    ens_avr = np.mean(ens, axis = 1)
    temp_avr = np.mean(temp_ens, axis=1)
    
    err = np.zeros(np.size(ens_avr, axis = 1)) #create array for errors 

    for k in range (np.size(experiments)):
        for j in range (np.size(ens_avr, 1)):
            err[j] = np.std(ens[k,:,j]) #calculate error
        up_err = (ens_avr[k,:])+err #calculate upper and lower bounds of error
        low_err = (ens_avr[k,:])-err
        axes[v].fill_between(t[50:], low_err[50:], up_err[50:], alpha=0.25, color=colours[k]) #plot shading of error
        axes[v].plot(t[50:], ens_avr[k,50:], colours[k], label = experiments[k]) #plot ensemble mean lines 
        
    axes[0].legend(fontsize=13) #plot legend 

#save figure 
fig.savefig('C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/Figures/Geoeng_paper/all_var_timeseries_Amazon_V1.png', bbox_inches='tight') 

