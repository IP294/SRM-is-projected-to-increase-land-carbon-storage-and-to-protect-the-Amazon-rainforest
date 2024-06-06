# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 12:55:01 2021

@author: impy2

Script to create time-series plots for specifc tipping points 
"""
### import libraries ###
from pathlib import Path
import iris.coords
import iris.cube
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

################# Specify the experiment and models ###############
experiments = ['G6sulfur', 'ssp585', 'ssp245']
colours = ['royalblue', 'firebrick','goldenrod']
models = ['CESM2-WACCM', 'CNRM-ESM2-1', 'IPSL-CM6A-LR', 'MPI-ESM1-2-LR', 'UKESM1-0-LL']
    
################# Specify variable of interest ######################
variables = ['tas', 'pr', 'NPP', 'cLand']
variables2 = ['tas', 'pr', 'npp','cland'] #lowercase variable names 

names = ['Near Surface Air Temperature', 'Precipitation Rates', 'Net Primary Productivity','Land Carbon'] #longform variable names
letters = ['a','b', 'c', 'd'] #panel labels 

### initialize ensemble zero arrays ###
ens = np.zeros((np.size(experiments),np.size(models), 241))
temp_ens = np.zeros((np.size(experiments),np.size(models), 241))

### initialize figure ###
fig, axes = plt.subplots(2,2, figsize=(15, 9))
axes = axes.flatten()

for v in range (np.size(variables)): #loop through variables
    var = variables[v] #select variables 
    var2 = variables2[v]
    variable_name = names[v] #seelect longform variable names 
    
    ########## Specify units of variable ############
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
        
    axes[v].set_xlabel('Year', fontsize=13) #set xlabel
    axes[v].set_ylabel('Change in '+variable_name+' ('+units+')') #select ylabel
    axes[v].set_title('Change in '+ variable_name) #set panel title
    axes[v].text(.02, 1.07, letters[v], ha='left', va='top', transform=axes[v].transAxes, fontsize='large') #label panel    
    
    for i in range (np.size(models)): #loop through models
        model = models[i] #select model
        
        for k in range (np.size(experiments)): #loop through experiments
            experiment = experiments[k] #select experiment 
            colour = colours[k] #select colour of experiment line 

            #### select date_range and variant_id corresponding to model ####   
            if model == 'CESM2-WACCM':
                if experiment == 'G6solar':
                    date_range = '2019-2100' ## geo
                    variant_id = 'r1i1p1f1'
                elif experiment == 'G6sulfur':
                    date_range = '2020-2100'
                    variant_id = 'r1i1p1f2'
                elif experiment == 'ssp585':
                    date_range = '2015-2300'   ## ssp585
                    variant_id = 'r1i1p1f1'
                elif experiment == 'ssp245':
                    date_range = '2015-2100'  ##ssp245
                    variant_id = 'r1i1p1f1'
                hist_dates = '1850-2014'
                hist_id = 'r1i1p1f1'
                
            elif model == 'CNRM-ESM2-1':
                if experiment == 'G6solar':
                    date_range = '2015-2100' ## geo
                elif experiment == 'G6sulfur':
                    date_range = '2015-2100'
                elif experiment == 'ssp585':
                    date_range = '2015-2100'   ## ssp585
                elif experiment == 'ssp245':
                    date_range = '2015-2100'  ##ssp245
                variant_id = 'r1i1p1f2'
                hist_dates = '1850-2350'
                hist_id = 'r1i1p1f2'
        
            elif model == 'IPSL-CM6A-LR':
                if experiment == 'G6solar':
                    date_range = '2020-2100' ## geo
                elif experiment == 'G6sulfur':
                    date_range = '2020-2100'
                elif experiment == 'ssp585':
                    if var =='cLand':
                        date_range = '2015-2100' ##ssp585
                    elif var== 'NPP' or var=='cVeg':
                        date_range = '2015-2300'
                elif experiment == 'ssp245':
                    date_range = '2015-2100'  ##ssp245
                variant_id = 'r1i1p1f1'
                hist_dates = '1850-2014'
                hist_id = 'r1i1p1f1'
        
            elif model == 'MPI-ESM1-2-LR':
                if experiment == 'G6solar':
                    date_range = '2015-2099' ## geo
                elif experiment == 'G6sulfur':
                    date_range = '2015-2099'
                elif experiment == 'ssp585':
                    date_range = '2015-2100'   ## ssp585
                elif experiment == 'ssp245':
                    date_range = '2015-2100'  ##ssp245
                variant_id = 'r1i1p1f1'
                hist_dates = '1850-2014'
                hist_id = 'r1i1p1f1'
        
            elif model =='UKESM1-0-LL':
                if experiment == 'G6solar':
                    date_range = '2020-2100' ## geo
                elif experiment == 'G6sulfur':
                    date_range = '2020-2100'
                elif experiment == 'ssp585':
                    date_range = '2015-2100'   ## ssp585
                elif experiment == 'ssp245':
                    date_range = '2015-2100'  ##ssp245'
                variant_id = 'r1i1p1f2' 
                hist_dates = '1850-2014'
                hist_id = 'r1i1p1f2'
                          
            ############ Specify name of directory of data which you want to plot #################
            region = 'World'
                      
            ################## Specify path to processed data directory ###########
            path = 'C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/'+var+'/historical/concatenated_data/'+region+'/'
           
            # File name of interpolated data
            fname = path+var+'_'+model+'_hist_'+experiment+'_'+variant_id+'_'+region+'.nc'
            fname = Path(fname)
            
            # Load in data
            try:
                f = nc.Dataset(fname,'r')
                
            except FileNotFoundError as e:
                print(model+' failed: {}'.format(e)) #exception if file fails to load
        
            else:
                x = f.variables[var2][:]
                #close dataset
                f.close()

                ##### Define latitude and longitude (uses 1 degree by 1 d) ######
                my_lats=np.arange(-90,91,step=1)
                my_lons=np.arange(-180,181,step=1)
            
                ny =  np.size(my_lats)
                nx = np.size(my_lons)
            
                #create new coordinates for iris cube 
                coord1 = iris.coords.DimCoord(my_lats,bounds=np.array([my_lats-0.5,my_lats+0.5]).T, standard_name='latitude', units='degrees', var_name='lat', attributes={'title': 'Latitude', 'type': 'double', 'valid_max': 90.0, 'valid_min': -90.0}, circular=True)
                coord2 = iris.coords.DimCoord(my_lons,bounds=np.array([my_lons-0.5,my_lons+0.5]).T, standard_name='longitude', units='degrees', var_name='lon', attributes={'title': 'Longitude', 'type': 'double', 'valid_max': 180.0, 'valid_min': -180.0}, circular=True)
                #create new iris cube
                cube = iris.cube.Cube(np.zeros((ny,nx),np.float32),dim_coords_and_dims=[(coord1,0),(coord2,1)])
                #extract latitude and longitude 
                latitude = coord1.points
                longitude = coord2.points
            
                # Obtain dimension sizes of data
                nt = int(len(x[:,0,0]))
                ny = int(len(x[0,:,0]))
                nx = int(len(x[0,0,:]))
                
                if var == 'cLand' or var == 'NPP': #load in land mask if variable is npp or cland
                ### Specify path to land mask ###
                    fname2 = 'C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/Masks/'+region+'.nc'
                    fname2 = Path(fname2)
                    
                    # Load in land mask and remove any ocean values
                    f2 = nc.Dataset(fname2,'r')
                    sftlf = f2.variables['sftlf'][:]
                    
                    sftlf_mask = np.broadcast_to(sftlf, (nt, ny,nx))
                    x[np.where(sftlf_mask <= 0)] = np.nan  
                    x[np.where(~np.isfinite(sftlf_mask))] = np.nan  #remove ocean values 
                    
                    areas = iris.analysis.cartography.area_weights(cube, normalize=False)
                    
                    land_area = np.sum(areas[np.where((sftlf>0)&(np.isfinite(sftlf)))]) #find total land area 
    
                    f2.close()
                    
                #ensure data is correct length
                data = np.empty((251, 181, 361)) #initialize empty data array of the correct size
                data.fill(np.nan) #fill with nans
                data[:np.size(x,axis=0), :,:] = x[:251,:,:] #fill array with data
                data[data > 1e36] = np.nan #remove fill values 
               
                ### convert units ###
                if var == 'pr':
                    data = data * 86400
                elif var == 'tas':
                    data = data -273
                elif var == 'NPP':
                    data = data * 3.154e7 #multiplied to per year 
                elif var == 'cVeg':
                    data = data * 1e-12 # convert to petagrams
                elif var == 'cLand':
                    data = data
                    
                t = np.arange(1855, np.size(data,0)-10 +1855) #create time array for x axis

                areas = iris.analysis.cartography.area_weights(cube, normalize=False)

                if var == 'cLand' or var == 'NPP':
                #### sum variable across the globe #####
                    vsum = np.zeros(np.size(data, 0)) #initialize zero array
                    
                    for j in range (np.size(data, 0)):
                        v_point = data[j,:,:] #select data from single time index                                                
                        vsum[j] = np.nansum(v_point*areas)*1e-12 #sum accross region and convert to petagrams 
                        
                        if vsum[j] == 0:
                            vsum[j] = np.nan #if sum is zero, replace with nan 
                    
                    #### calculate the anomaly of the y variable relative to the first 50 years of the hist run(if cveg or cland)
                    v_anom = vsum - np.nanmean(vsum[0:50])

                elif var == 'tas' or var == 'pr': #if variable is temp or precip, average across region
                    vavr = np.nansum((data*areas),axis=1)/np.nansum(areas) # finding the weighted average 
                    vavr = np.nansum(vavr, axis=1)
                    vavr[vavr==0] = np.nan #if sum is zero, replace with nan
                        
                    v_anom = vavr - np.nanmean(vavr[0:50]) #find anomaly relative to the preindustrial period
                     
            ### taking a running decadal mean - centred on the middle of the window  ###
            smoothed_anom = np.zeros(len(v_anom)-10) #initialize zero array
            
            for j in range (len(v_anom)-10):
                m1 = j #start of window                
                m2 = j+10 #end of window
                smoothed_anom[j] = np.mean(v_anom[m1:m2]) #find average of sliding window
                
            #### save data into ens for later averaging ####      
            ens[k,i,:] = smoothed_anom
    
        axes[k].tick_params(labelcolor="black", bottom=True, left=True)
        fig.subplots_adjust(hspace=0.5) #adjust figure spacing

    ### calculate ensemble averages ###
    ens_avr = np.mean(ens, axis = 1)
    temp_avr = np.mean(temp_ens, axis=1)
    
    err = np.zeros(np.size(ens_avr, axis = 1)) #initialize zero array for error

    for k in range (np.size(experiments)):
        for j in range (np.size(ens_avr, 1)):
            err[j] = np.std(ens[k,:,j]) #calculate error
        up_err = (ens_avr[k,:])+err #upper bound of error
        low_err = (ens_avr[k,:])-err #lower bound of error
        axes[v].fill_between(t[50:], low_err[50:], up_err[50:], alpha=0.25, color=colours[k]) #plot error as shaded error
        axes[v].plot(t[50:], ens_avr[k,50:], colours[k], label = experiments[k]) #plot ensemble means 
        axes[v].set_xlim(2000,2100) #set x limits 
    axes[0].legend(fontsize=13) #plot legend

#save figure 
fig.savefig('C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/Figures/Geoeng_paper/all_var_timeseries_V2.png', bbox_inches='tight') 

