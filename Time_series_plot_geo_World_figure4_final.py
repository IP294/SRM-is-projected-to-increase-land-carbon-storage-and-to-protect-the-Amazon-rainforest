# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 12:55:01 2021

@author: impy2

Script to create time-series plots for specifc tipping points 
"""

######### import libraries #############
from pathlib import Path
import iris.coords
import iris.cube
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from pandas import read_csv

def find_nearest(array, value):
    """
    Function to find the index of the closest value in an array to a chosen number
    Input: Array - to be searched for closest value, Value - number searched for
    Returns: index of the nearest value in array
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx
    
################# Specify the experiment and models ################
experiments = ['G6sulfur', 'ssp585', 'ssp245']
colours = ['royalblue', 'firebrick','goldenrod', 'forestgreen', 'indigo']
letters = ['a', 'b', 'c', 'd'] # letters for labelling panels

################# Specify model here ################################## 
models = ['CESM2-WACCM', 'CNRM-ESM2-1', 'IPSL-CM6A-LR', 'MPI-ESM1-2-LR', 'UKESM1-0-LL'] #models for analysis
names = ['net primary productivity','land carbon'] #longform names of variables
names2 = ['NPP','land carbon'] #shortform names of variables 

variables = ['NPP','cLand'] #variables for analysis
variables2 = ['npp', 'cland'] # same variables 

# x axis variables
var4 = 'tas'
var5 = 'co2'

###### initlaize zero arrays ########
ens = np.zeros((np.size(experiments),np.size(models), 251))
sorted_ens = np.zeros((np.size(experiments),np.size(models), 251))
indexes = np.zeros((np.size(experiments),np.size(models), 251))
temp_ens = np.zeros((np.size(experiments),np.size(models), 241))

######## initialize figure ###########
fig, axes = plt.subplots(2,2, figsize=(15, 9))
axes = axes.flatten()

for v in range (np.size(variables)): # loop through each variable
    var = variables[v] #select variable
    var2 = variables2[v]
    
    for i in range (np.size(models)): #loop through each model
        model = models[i] # select model
        
        for k in range (np.size(experiments)): #loop through experiments 
            experiment = experiments[k] # select experiment
            
            ############## Specify units of variable #################
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
                
            ############ Specify name of directory of data which you want to plot #################
            region = 'World'
            
            ############ Specify name of directory of that stores the figures #####
            region2 = 'World'
            
            # Window length used to calculate abrupt shift over
            wl = 15
            
            ################## Specify path to processed data directory ###########
            #define file path to data 
            path = 'C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/'+var+'/historical/concatenated_data/'+region+'/'
            
            # File name of interpolated data
            fname = path+var+'_'+model+'_hist_'+experiment+'_'+variant_id+'_'+region+'.nc' # name of file path to data
            fname = Path(fname)
            
            ############ Load in data ###########
            try:
                f = nc.Dataset(fname,'r') 
                
            except FileNotFoundError as e: #exception if file fails to load
                print(model+' failed: {}'.format(e))
        
            else:
                # Extract data and if measured in per second convert to per year
                x = f.variables[var2][:]
                
                #Close dataset
                f.close()
                
                ### if the variable is cland or NPP then load in a land mask ###
                if var == 'cLand' or var == 'NPP':
                
                    # Obtain dimension sizes
                    nt = int(len(x[:,0,0]))
                    ny = int(len(x[0,:,0]))
                    nx = int(len(x[0,0,:]))
                    
                ############ Specify path to land mask ######################
                    fname2 = 'C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/Masks/'+region+'.nc'
                    fname2 = Path(fname2)
                    
                    # Load in land mask 
                    f2 = nc.Dataset(fname2,'r')
                    sftlf = f2.variables['sftlf'][:] 
                    sftlf_mask = np.broadcast_to(sftlf, (nt, ny,nx))
                    
                    # remove ocean values from data 
                    x[np.where(sftlf_mask <= 0)] = np.nan  
                    x[np.where(~np.isfinite(sftlf_mask))] = np.nan
    
                    f2.close() # close land mask dataset
       
                data = x # rename x as data
                
                # ensure dataset is 252 in length
                app = 251 - np.size(data,0)
                if app > 0:
                    data = np.insert(data, -1, np.nan, axis=0)
                    data[-1,:,:] = np.nan
                    data[data>1e36] = np.nan
               
                ##### Define latitude and longitude (uses 1x1 degree grid) ######
                my_lats=np.arange(-90,91,step=1)
                my_lons=np.arange(-180,181,step=1)
                
                #define dimensions to be used in iris cube manipulations
                ny =  np.size(my_lats) 
                nx = np.size(my_lons)
                
                #create new iris coordinates for latitude and longitude
                coord1 = iris.coords.DimCoord(my_lats,bounds=np.array([my_lats-0.5,my_lats+0.5]).T, standard_name='latitude', units='degrees', var_name='lat', attributes={'title': 'Latitude', 'type': 'double', 'valid_max': 90.0, 'valid_min': -90.0}, circular=True)
                coord2 = iris.coords.DimCoord(my_lons,bounds=np.array([my_lons-0.5,my_lons+0.5]).T, standard_name='longitude', units='degrees', var_name='lon', attributes={'title': 'Longitude', 'type': 'double', 'valid_max': 180.0, 'valid_min': -180.0}, circular=True)
                cube = iris.cube.Cube(np.zeros((ny,nx),np.float32),dim_coords_and_dims=[(coord1,0),(coord2,1)]) #create new iris cube from coordinates
            
                areas = iris.analysis.cartography.area_weights(cube, normalize=False) #create areas array for weighting data
               
                ###### convert units depending on variable ######
                if var == 'pr':
                    data = data * 86400
                elif var == 'tas':
                    data = data -273
                elif var == 'NPP':
                    data = data * 3.154e7 *1e-12 #multiplied to make per year and converted to Pg
                elif var == 'cVeg':
                    data = data * 1e-12 # convert to petagrams
                elif var == 'cLand':
                    data = data* 1e-12
                    
                ############## calcualte the temperature anomolies over time ################
                if experiment == 'ssp245' or experiment == 'G6sulfur':
                    x_temps = np.linspace(0, 3.5, 251) # create an array of temperatures to act as the x axis
                elif experiment == 'ssp585':
                    x_temps = np.linspace(0, 6, 251) # create an array of temperatures to act as the x axis
                
                ##### specify file path to temperature data #####
                path4 = 'C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/'+var4+'/historical/concatenated_data/World'
                fname = path4+'/'+var4+'_'+model+'_hist_'+experiment+'_'+variant_id+'_World.nc'
                
                ### Load in world data into iris cube ###
                world_cube = iris.load_cube(fname,var4)
                world_data = world_cube.data
                
                ### ensure data is 250 in length ###
                app = 251 - np.size(world_data,0)
                if app > 0:
                    world_data = np.insert(world_data, -1, np.nan, axis=0)
                    world_data[-1,:,:] = np.nan
                
                world_data[world_data>1e36] = np.nan # remove filler values in data 
                
                ### area weight the temperature ###
                weighted_tas = np.nansum((world_data*areas),axis=1)/np.nansum(areas)
                weighted_tas = np.nansum(weighted_tas, axis=1)


                ### cut off longer datasets at 250 years ###
                if model =='CESM2-WACCM' and experiment == 'ssp585':
                    weighted_tas = weighted_tas[:252]
                elif model == 'IPSL-CM6A-LR' and experiment == 'ssp585':
                    weighted_tas = weighted_tas[:251]    

                ####find average temperature from 1850-1899 ####
                weighted_hist = weighted_tas[0:50] #select temp reference period 
                hist_ref = np.nanmean(weighted_hist) # average for temp reference 

                anom = np.zeros(np.size(weighted_tas)-10) # define zero array for temperature anomaly 
                
                ### calculate the temperaure anomaly for each 10 year window going forward
                for j in range (0, np.size(weighted_tas)-10): #looping through each model year
                    temp = np.average(weighted_tas[j:j+10])
                    anom[j] = temp - hist_ref
              
                #### find the index of temperatures closest to the x_temp array
                for j in range (len(x_temps)):
                    temp_index = find_nearest(anom, x_temps[j])
                    indexes[k,i,j]  = temp_index
  
                weighted_data = data*areas ##### area weight the variable data #####
                
                #################### total of variable across the globe ####################
                vsum = np.zeros(np.size(weighted_data, 0)) #initialise zero array
                
                for j in range (np.size(weighted_data, 0)):
                    v_point = weighted_data[j,:,:]
                    v_point[v_point > 1e36] = np.nan #remove filler values in the data 
                    
                    vsum[j] = np.nansum(v_point) # calculate global sum
                    if vsum[j] == 0: #remove zero values from the sum
                        vsum[j] = np.nan
                
                #### calculate the anomaly of the y variable relative to the first 50 years of the hist run(if cveg or cland)
                vsum_anom = vsum - np.nanmean(vsum[0:50])
                
                ens[k,i,:] = vsum_anom # save anomaly into ensemble array
                indexes = indexes.astype('int32') 
                sorted_ens[k,i,:] = vsum_anom[indexes[k,i,:]] # order data according to temperature indexes

    ### calculate ensemble averages ###
    ens_avr = np.mean(ens, axis = 1)
    sorted_ens_avr = np.mean(sorted_ens, axis=1)
    err = np.zeros(np.size(ens_avr, axis = 1))
    s_err = np.zeros(np.size(ens_avr, axis = 1))

    
    fig.subplots_adjust(hspace=0.3) # adjust spacing of panels 
    
    ########## plot figure ###########
    for k in range (np.size(experiments)):
        experiment = experiments[k] #select experiment for plotting
        
        if experiment == 'ssp245' or experiment == 'G6sulfur':
            x_temps = np.linspace(0, 3.5, 251) # create an array of temperatures to act as the x axis
        elif experiment == 'ssp585':
            x_temps = np.linspace(0, 6, 251) # create an array of temperatures to act as the x axis
        
        ################ Importing co2 arrays for x axis ###################
        if experiment == 'ssp585' or experiment == 'G6sulfur':
            co2_csv = read_csv("C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/co2/historical/Processed_data/World/ssp585_CO2.csv")
        elif experiment == 'ssp245':
            co2_csv = read_csv("C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/co2/historical/Processed_data/World/ssp245_CO2.csv")
        
        co2 = np.array(co2_csv['co2']) # co2 x axis array
        
        ###### calculate errors ######
        for j in range (np.size(ens_avr, 1)):
            err[j] = np.std(ens[k,:,j]) # error is one standard deviation 
        
        #calculate upper and lower bounds of error for co2 array
        up_err = (ens_avr[k,:])+err
        low_err = (ens_avr[k,:])-err
        
        #calculate upper and lower bpunds of error for sorted temp array
        s_up_err = (sorted_ens_avr[k,:])+err
        s_low_err = (sorted_ens_avr[k,:])-err
        
        # plot data onto figures 
        axes[2*v].fill_between(x_temps, s_low_err, s_up_err, alpha=0.25, color=colours[k])
        axes[2*v].plot(x_temps, sorted_ens_avr[k,:], colours[k], label = experiments[k])   
        
        axes[2*v+1].fill_between(co2[60:], low_err[60:], up_err[60:], alpha=0.25, color=colours[k])
        axes[2*v+1].plot(co2[50:], ens_avr[k,50:], colours[k], label= experiments[k])
            
    axes[2*v].text(.02, 1.07, letters[v], ha='left', va='top', transform=axes[v].transAxes, fontsize='large') #label panels 
    axes[2*v+1].text(.02, 1.07, letters[v+2], ha='left', va='top', transform=axes[v+2].transAxes, fontsize='large')

    axes[0].legend(fontsize=13) # create legend
    
    ## label axes ##
    axes[2*v].set_xlabel('Temperature anomaly ($^\circ$C)')
    axes[2*v].set_ylabel('Change in total '+names2[v]+' ('+units+')')
    axes[2*v].set_title('Change in '+ names2[v]+' versus Global Warming')
    
    axes[2*v+1].set_xlabel('Atmospheric CO$_{2}$ concentration (ppm)')
    axes[2*v+1].set_ylabel('Change in total '+names2[v]+' ('+units+')')
    axes[2*v+1].set_title('Change in '+ names2[v]+' versus CO$_{2}$')

### save figures ###
fig.savefig('C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/Figures/Geoeng_paper/temp_co2_anom_timeseries_V1-alt.png', bbox_inches='tight') 