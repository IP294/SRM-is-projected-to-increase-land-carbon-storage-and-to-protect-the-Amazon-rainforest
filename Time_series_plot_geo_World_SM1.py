# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 12:55:01 2021

@author: impy2

Script to create time-series plots for specifc tipping points 
"""
### Import libraries ###
from pathlib import Path
import iris.coords
import iris.cube
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

################# Specify the experiment and models ################
experiments = ['G6sulfur','ssp585', 'ssp245']
colours = ['royalblue', 'goldenrod', 'forestgreen', 'firebrick', 'indigo'] # specify colours for each model in plot 
models = ['CESM2-WACCM', 'CNRM-ESM2-1', 'IPSL-CM6A-LR', 'MPI-ESM1-2-LR', 'UKESM1-0-LL']
letters = ['a', 'b', 'c'] #specify letter labels for panels of figure 


############## specify variable #####################
var = 'cLand' #variable to be plotted
var2 = 'cland' #lower case version of var for file path

############## Specify units of variable ################
if var == 'cVeg':
    units = '$Pg$'
elif var == 'tas':
    units = '$^\circ$C'
elif var == 'pr':
    units = '$mm day^{-1}$'
elif var == 'NPP':
    units = '$Kg[C] yr^{-1}$'
elif var == 'cSoil' or var == 'cLand':
    units = '$Pg$'

### initialize ensemble array ###
ens = np.zeros((np.size(experiments),np.size(models), 81))

### initialize figure ###
fig, axes = plt.subplots(1,3, figsize=(21, 4))
axes = axes.flatten()

for i in range (np.size(models)): # loop through all models 
    colour = colours[i] # select colour for the model line
    model = models[i] # select model
    
    for k in range (np.size(experiments)): # loop through each experiment
        experiment = experiments[k] # select experiment 
    
        axes[k].set_xlabel('Year') # set x axis label 
        axes[k].text(.02, 1.07, letters[k], ha='left', va='top', transform=axes[k].transAxes, fontsize='large') #label the panel 
        axes[k].set_title(experiment) # set panel title
        axes[k].set_ylabel('Change in Land Carbon ('+units+')') # set y axis label
               
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
        
        ############ Specify name of directory of that stores the figures #################
        region2 = 'World'
        
        # Window length used to calculate abrupt shift over
        wl = 15
        nyr = 86 
        
        ################## Specify path to processed data directory #####################
        path = 'C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/'+var+'/'+experiment+'/Processed_data/'+region+'/'
        # File name of interpolated data
        fname = path+var+'_'+model+'_'+experiment+'_'+variant_id+'_'+date_range+'_'+region+'.nc'
        fname = Path(fname)
        
        # Load in data
        try:
            f = nc.Dataset(fname,'r')
            
        except FileNotFoundError as e: #exception if file fails to load
            print(model+' failed: {}'.format(e))
    
        else:
            x = f.variables[var2][:]
            
            #cut off data for models that run for longer than 86 years
            if model == 'CESM2-WACCM' and experiment =='ssp585':
                x = x[:86,:,:]
            elif model == 'IPSL-CM6A-LR' and experiment == 'ssp245':
                x = x[:86,:,:]
                
            # Close dataset
            f.close()
        
            data = x # rename x variable as 'data'
            
            ##### Define latitude and longitude (uses 1 degree by 1 d) ######
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
            elif var == 'npp':
                data = data * 3.154e7 #multiplied to per year
            elif var == 'cVeg':
                data = data * 1e-12 # convert to petagrams
            elif var == 'cLand':
                data = data* 1e-12

            ##### area weight the variable data #####
            weighted_data =np.zeros((np.size(data,0), np.size(data,1), np.size(data,2))) # initialise array for area weighted data
            
            for j in range (0, np.size(data,0)):
                weighted_data[j] = data[j,:,:]*areas # weight data by areas 
            
            #################### Find the average of variable across the globe ####################
            vsum = np.zeros(np.size(weighted_data, 0)) # initliase arrays for loop
            
            for j in range (np.size(weighted_data, 0)):
                v_point = weighted_data[j,:,:]
                v_point[v_point > 1e36] = np.nan # remove fill values from array
                
                vsum[j] = np.nansum(v_point) # sum array for global total
            
            #### calculate the anomaly of the y variable relative to the first 10 years (if cveg or cland)
            vsum_anom = vsum - np.mean(vsum[0:10])
                    
            t = np.arange(2015, np.size(data,0)+2015) #create time array for x axis 
          
            ###### plot line for this experiment ######
            axes[k].plot(t, vsum_anom, colour, label=model)
            
        #### save data into ens for later averaging ####
        ens[k,i,:] = vsum_anom[:81]

    axes[0].legend() # create legend for models 
    
    axes[k].tick_params(labelcolor="black", bottom=True, left=True)
    fig.subplots_adjust(hspace=0.5) #adjust spacing of plots 
    
### save figure ###
fig.savefig('C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/Figures/Geoeng_paper/all_model_timeseries.png', bbox_inches='tight')    
