# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 10:14:22 2020

@author: ip294

Create Figures 2&3
"""

### import libraries ###
from pathlib import Path
import iris.coords
import iris.cube
import netCDF4 as nc
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.mpl.ticker as cticker
import matplotlib.colors as mcolors

experiment1 = 'ssp245' #may be set as ssp245 or ssp585
experiment2 = 'G6sulfur' #SAI experiment compared to 
experiments = ['ssp585', 'ssp245']

################# Specify model here ##################################
models = ['CESM2-WACCM', 'CNRM-ESM2-1', 'IPSL-CM6A-LR', 'MPI-ESM1-2-LR', 'UKESM1-0-LL'] #models to be analysed 

### initialize zero arrays ###
futrs1 = np.zeros((np.size(models), 34,49))
futrs2 = np.zeros((np.size(models), 34,49))
rels = np.zeros((np.size(models), 34,49))
abss = np.zeros((np.size(models), 34,49))
abss2 = np.zeros((np.size(models), 34,49))
lands = np.zeros((2,86,34,49))

################# Specify variable of interest ######################
varss = ['tas', 'pr','NPP','cLand'] #variables to plot
vars2 = ['tas', 'pr', 'npp', 'cland']# var2 is the same as var but all lower case

land_var = 'primf' #land use variable for land use data
land_var2 = 'secdf' # second land use variable for land use data

letters = ['a', 'b', 'c', 'd'] #labels for panels 
    
var_names = ['near surface air temperature', 'precipitation flux', 'net primary productivity','land carbon'] #long names of variables 

############ Specify name of directory of data which you want to plot #################
region = 'Amazon'
  
if region == 'World':
    lonmin = -180-0.5
    lonmax = 180+0.5
    latmin = -90-0.5
    latmax = 90+0.5
    lonstep = 60
    latstep = 20
    figsizes = (20, 13)
    
elif region == 'Amazon':
    lonmin = -83-0.5
    lonmax = -35+0.5
    latmin = -20-0.5
    latmax = 13+0.5
    lonstep = 10
    latstep = 10
    figsizes = (17, 13)
 
### initialize figure ###
proj = ccrs.PlateCarree() #select map projection
fig, axes = plt.subplots(2,2, subplot_kw=dict(projection=proj), figsize=figsizes, sharex=True, sharey=True)
axes = axes.flatten()

##### Load in land use data to mask out areas with prescribed land use changes #####
################## Specify path to processed data directory #####################
path = 'C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/land_use/Processed_data/'

# Longitude, latitude min, max, step and figure size for region of interest
for i in range (len(experiments)):
    land_experiment = experiments[i]
    
    # File name of interpolated data
    fname = path+region+'/'+land_var+'_'+land_experiment+'_states_landuse.nc'
    fname = Path(fname)
    
    fname2 = path+region+'/'+land_var2+'_'+land_experiment+'_states_landuse.nc'
    fname2 = Path(fname2)
    
    # Load in data
    f1 = nc.Dataset(fname,'r')
    f2 = nc.Dataset(fname2,'r')
    
    # Extract longitude and latitude
    longitude = f1.variables['longitude'][:]
    latitude = f1.variables['latitude'][:]
    
    x1 = f1.variables[land_var][:,:,:]
    x2 = f2.variables[land_var2][:,:,:]
    
    # Close dataset
    f1.close()
    f2.close()
    
    x = x1+x2 # add two land use variables together 
    lands[i,:,:,:] = x ## save land use data in ensemble array

nt = int(len(x1[:,0,0])) #find length of data 
land1 = np.nanmean(lands[0][nt-10:nt,:,:],axis=0) #find mean of last 10 years of ssp585
land2 = np.nanmean(lands[1][nt-10:nt,:,:],axis=0) #find mean of last 10 years of ssp245

land_change = (land2 - land1) #finding difference between ssp585 and ss245 prescribed land use changes

land_mask_ind = np.where(land_change>=0.1) #find where difference in prescribed land use is greater than 0.1
land_use_mask =ma.masked_where(land_change >=0.1, land_change) #mask areas found 
land_use_mask = land_use_mask.mask*1


for k in range (len(varss)): #loop through variables 
    var = varss[k] #select variable 
    var2 = vars2[k]
    
    ############## Specify units of variable ##################
    if var == 'cVeg':
        units = 'kgC m$^{-2}$'
    elif var == 'tas':
        units = '$^\circ$C'
    elif var == 'pr':
        units = 'mm day$^{-1}$'
    elif var == 'NPP':
        units = 'kgC m$^{-2}$ yr$^{-1}$'
    elif var == 'cSoil' or var == 'cLand':
        units = 'kgC m$^{-2}$'
    
    for m in range(np.size(models)): #loop through model        
        model=models[m] #select model 
        
        ##### select variant ID and date ranges based on experiment ####
        if model == 'CESM2-WACCM':
            date_range1 = '2015-2100'  ##ssp245
            date_range2 = '2020-2100' ## geo
            variant_id1 = 'r1i1p1f1'
            variant_id2 = 'r1i1p1f2'
            
        elif model == 'CNRM-ESM2-1':
            date_range1 = '2015-2100'
            date_range2 = '2015-2100'
            variant_id1 = 'r1i1p1f2'
            variant_id2 = 'r1i1p1f2'
    
        elif model == 'IPSL-CM6A-LR':
            if var =='cLand':
                date_range1 = '2015-2100' ##ssp585
            else:
                date_range1 = '2015-2100'
            date_range2 = '2020-2100' ##Geo
            variant_id1 = 'r1i1p1f1'
            variant_id2 = 'r1i1p1f1'
    
        elif model == 'MPI-ESM1-2-LR':
            date_range1 = '2015-2100' ##ssp585 /ssp245
            date_range2 = '2015-2099' ##Geo
            variant_id1 = 'r1i1p1f1'
            variant_id2 = 'r1i1p1f1'
            
        elif model == 'MPI-ESM1-2-HR':
            date_range1 = '2015-2099' ##ssp585 /ssp245
            date_range2 = '2015-2099' ##Geo
            variant_id1 = 'r1i1p1f1'
            variant_id2 = 'r1i1p1f1'
    
        elif model =='UKESM1-0-LL':
            date_range1 = '2015-2100' ## ssp585 / ssp245
            date_range2 = '2020-2100' ## Geo
            variant_id1 = 'r1i1p1f2' 
            variant_id2 = 'r1i1p1f2'
        
        ################## Specify path to processed data directory #####################
        path1 = 'C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/'+var+'/historical/concatenated_data/'+region+'/'

        # File name of interpolated data
        #1 ssp585 or ssp245
        fname1 = path1+var+'_'+model+'_hist_'+experiment1+'_'+variant_id1+'_'+region+'.nc'
        fname1 = Path(fname1)
        #2 G6
        fname2 = path1+var+'_'+model+'_hist_'+experiment2+'_'+variant_id2+'_'+region+'.nc'
        fname2 = Path(fname2)
    
        # Load in data
        try:
            f1 = nc.Dataset(fname1,'r')
            f2 = nc.Dataset(fname2,'r')
            
        except FileNotFoundError as e: #exception if data fails to load 
            print(model+' failed: {}'.format(e))
    
        else:
            # Extract longitude and latitude 
            longitude = f1.variables['lon'][:]
            latitude = f1.variables['lat'][:]
            
            #extract data for variable 
            if model == 'IPSL-CM6A-LR' and var=='tas': #cut off tas array for IPSL which is 451 long
                x1 = f1.variables[var2][:251,:,:]
            else:
                x1 = f1.variables[var2][:]
            x2 = f2.variables[var2][:]
           
            # Close dataset
            f1.close()
            f2.close()            
            
            #create arrays of latitude and longitude 
            my_lats=np.arange(-20,14,step=1)
            my_lons=np.arange(-83,-34,step=1)
        
            #find lengths of lat and lon arrays 
            ny =  np.size(my_lats) 
            nx = np.size(my_lons)
            
            #create new iris cube coordinates 
            coord1 = iris.coords.DimCoord(my_lats,bounds=np.array([my_lats-0.5,my_lats+0.5]).T, standard_name='latitude', units='degrees', var_name='lat', attributes={'title': 'Latitude', 'type': 'double', 'valid_max': 90.0, 'valid_min': -90.0}, circular=True)
            coord2 = iris.coords.DimCoord(my_lons,bounds=np.array([my_lons-0.5,my_lons+0.5]).T, standard_name='longitude', units='degrees', var_name='lon', attributes={'title': 'Longitude', 'type': 'double', 'valid_max': 180.0, 'valid_min': -180.0}, circular=True)
            cube = iris.cube.Cube(np.zeros((ny,nx),np.float32),dim_coords_and_dims=[(coord1,0),(coord2,1)])
            
            #global grid of areas for weightings 
            areas = iris.analysis.cartography.area_weights(cube, normalize=False)
        
            latitude = coord1.points
            longitude = coord2.points
            
            #### convert uints of variables ####
            if var == 'NPP':
                x1 = x1* 3.154e7 #multiplied to make per year
                x2 = x2* 3.154e7 #multiplied to make per year
            elif var == 'pr':
                x1 = x1*86400 #converting units 
                x2 = x2*86400
             
            # Obtain dimension sizes
            nt = int(len(x1[:,0,0]))
            ny = int(len(x1[0,:,0]))
            nx = int(len(x1[0,0,:]))
               
            if var == 'cLand' or var == 'NPP': #land use map required only fot cland and NPP
            ############ Specify path to land mask ######################
                fname2 = 'C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/Masks/'+region+'.nc'
                fname2 = Path(fname2)
                
                # Load in land mask and remove any ocean values
                f2 = nc.Dataset(fname2,'r')
                sftlf = f2.variables['sftlf'][:]
                
                sftlf_mask = np.broadcast_to(sftlf, (nt, ny,nx))
                x1[np.where(sftlf_mask <= 0)] = np.nan  
                x1[np.where(~np.isfinite(sftlf_mask))] = np.nan #mask out ocean values 
                
                sftlf_mask = np.broadcast_to(sftlf, (int(len(x2[:,0,0])), ny,nx))
                x2[np.where(sftlf_mask <= 0)] = np.nan  
                x2[np.where(~np.isfinite(sftlf_mask))] = np.nan #mask out ocean values 
                
                land_area = np.nansum(areas[np.where((sftlf>0)&(np.isfinite(sftlf)))]) #find total land area in region
                
                f2.close() #close mask dataset
            
            x_futr1 = np.nanmean(x1[-11:-1,:,:],axis=0)
            if model == 'MPI-ESM1-2-LR':
                x_futr2 = np.nanmean(x2[-10:-1,:,:], axis=0)
            else:
                x_futr2 = np.nanmean(x2[-11:-1,:,:], axis=0)

            x_change = x_futr2 - x_futr1
        
            ## mask out areas with changes in prescribed land use change ###
            x_change[land_mask_ind] = np.nan
            x_futr1[land_mask_ind] = np.nan
            x_futr2[land_mask_ind] = np.nan
            
            ### save into ensemble arrays ###
            abss[m,:,:] = x_change
            futrs1[m,:,:] = x_futr1
            futrs2[m,:,:] = x_futr2
            
            ### Create colour maps  ###           
            if var == 'cLand' or var=='cVeg' or var =='cSoil':
                vc = np.concatenate([np.arange(-5, 0, 1.0), np.arange(1, 6, 1.0)])
                levs = range(12)
                cmap = mcolors.LinearSegmentedColormap.from_list(name='red_white_green', 
                                                                 colors =[(0.9, 0.25, 0.21), 
                                                                          (1, 1., 1), 
                                                                          (0.25, 0.67, 0.35)],
                                                                 N=len(levs)-1)
                colors1 = list(cmap(np.arange(len(vc)+1)))
                cmap = colors.ListedColormap(colors1[1:-1], "")
                cmap.set_under(colors1[0])     # set under-color to last color of list 
                cmap.set_over(colors1[-1])     # set over-color to firstst color of list 
                norm = colors.BoundaryNorm(vc, cmap.N)                
                
                # Define colour intervals, and colour maps of 30 year mean plots
                vr = np.arange(0, 30, step=5)
                cmap_name2 = 'Greens'#'RdYlBu_r'#'YlGnBu'#'Greens'#'PRGn'#'RdBu_r'
                cmap2 = plt.cm.get_cmap(cmap_name2,len(vr))
                if var == 'cVeg':
                    colors2 = list(cmap2(np.arange(len(vr))))
                    cmap2 = colors.ListedColormap(colors2[:-1], "")
                    cmap2.set_over(colors2[-1])     # set over-color to last color of list 
                norm2 = colors.BoundaryNorm(vr, cmap2.N)
                
                # color map for relative change 
                vrc = np.arange(-50, 60, 10)
                cmap_name3 = 'RdYlGn'
                cmap3 = plt.cm.get_cmap(cmap_name3,len(vrc))
                norm3 = colors.BoundaryNorm(vrc, cmap3.N)
                
            elif var == 'NPP':
                vc = np.concatenate([np.arange(-0.5, 0, 0.1), np.arange(0.1, 0.6, 0.1)])
                levs = range(12)
                cmap = mcolors.LinearSegmentedColormap.from_list(name='red_white_green', 
                                                                 colors =[(0.9, 0.25, 0.21), 
                                                                          (1, 1., 1), 
                                                                          (0.25, 0.67, 0.35)],
                                                                 N=len(levs)-1)
                colors1 = list(cmap(np.arange(len(vc)+1)))
                cmap = colors.ListedColormap(colors1[1:-1], "")
                cmap.set_under(colors1[0])     # set under-color to last color of list 
                cmap.set_over(colors1[-1])     # set over-color to firstst color of list 
                norm = colors.BoundaryNorm(vc, cmap.N)
    
                # Define colour intervals, and colour maps of 30 year mean plots
                vr = np.arange(0, 2.2, step=0.2)#(0, 120, step=20)#(275, 320, step=5)#(500, 4000, step=500)#np.arange(0, 120, step=20)#(0, 30, step=5)
                cmap_name2 = 'Greens'#'RdYlBu_r'#'YlGnBu'#'Greens'#'PRGn'#'RdBu_r'
                cmap2 = plt.cm.get_cmap(cmap_name2,len(vr))
                if var == 'cVeg':
                    colors2 = list(cmap2(np.arange(len(vr))))
                    cmap2 = colors.ListedColormap(colors2[:-1], "")
                    cmap2.set_over(colors2[-1])     # set over-color to last color of list 
                norm2 = colors.BoundaryNorm(vr, cmap2.N)
                
                # color map for relative change 
                vrc = np.arange(-50, 60, 10)
                cmap_name3 = 'RdYlGn'
                cmap3 = plt.cm.get_cmap(cmap_name3,len(vrc))
                norm3 = colors.BoundaryNorm(vrc, cmap3.N)
            
            elif var == 'tas':
                vc = np.concatenate([np.arange(-5, 0, 1), [-0.5], [0.5], np.arange(1, 6, 1)])
                cmap_name = 'RdBu_r'#'RdBu'#'RdYlGn'#'PRGn'#'RdBu_r'
                cmap = plt.get_cmap(cmap_name,len(vc)+1)
                colors1 = list(cmap(np.arange(len(vc)+1)))
                cmap = colors.ListedColormap(colors1[1:-1], "")
                cmap.set_under(colors1[0])     # set under-color to last color of list 
                cmap.set_over(colors1[-1])     # set over-color to firstst color of list 
                norm = colors.BoundaryNorm(vc, cmap.N)
    
                # Define colour intervals, and colour maps of 30 year mean plots
                vr = np.arange(0, 2.2, step=0.2)#(0, 120, step=20)#(275, 320, step=5)#(500, 4000, step=500)#np.arange(0, 120, step=20)#(0, 30, step=5)
                cmap_name2 = 'RdBu_r'#'RdYlBu_r'#'YlGnBu'#'Greens'#'PRGn'#'RdBu_r'
                cmap2 = plt.cm.get_cmap(cmap_name2,len(vr))
                if var == 'cVeg':
                    colors2 = list(cmap2(np.arange(len(vr))))
                    cmap2 = colors.ListedColormap(colors2[:-1], "")
                    cmap2.set_over(colors2[-1])     # set over-color to last color of list 
                norm2 = colors.BoundaryNorm(vr, cmap2.N)
                
                # color map for relative change 
                vrc = np.arange(-50, 60, 10)
                cmap_name3 = 'RdBu_r'
                cmap3 = plt.cm.get_cmap(cmap_name3,len(vrc))
                norm3 = colors.BoundaryNorm(vrc, cmap3.N)
            
            elif var == 'pr':
                vc = np.concatenate([np.arange(-1, 0, 0.2), np.arange(0.2, 1.2, 0.2)])
                cmap_name = 'BrBG'#'RdBu'#'RdYlGn'#'PRGn'#'RdBu_r'
                cmap = plt.get_cmap(cmap_name,len(vc)+1)
                colors1 = list(cmap(np.arange(len(vc)+1)))
                cmap = colors.ListedColormap(colors1[1:-1], "")
                cmap.set_under(colors1[0])     # set under-color to last color of list 
                cmap.set_over(colors1[-1])     # set over-color to firstst color of list 
                norm = colors.BoundaryNorm(vc, cmap.N)
    
                # Define colour intervals, and colour maps of 30 year mean plots
                vr = np.arange(0, 2.2, step=0.2)#(0, 120, step=20)#(275, 320, step=5)#(500, 4000, step=500)#np.arange(0, 120, step=20)#(0, 30, step=5)
                cmap_name2 = 'PRGn'#'RdYlBu_r'#'YlGnBu'#'Greens'#'PRGn'#'RdBu_r'
                cmap2 = plt.cm.get_cmap(cmap_name2,len(vr))
                if var == 'cVeg':
                    colors2 = list(cmap2(np.arange(len(vr))))
                    cmap2 = colors.ListedColormap(colors2[:-1], "")
                    cmap2.set_over(colors2[-1])     # set over-color to last color of list 
                norm2 = colors.BoundaryNorm(vr, cmap2.N)
                
                # color map for relative change 
                vrc = np.arange(-50, 60, 10)
                cmap_name3 = 'PRGn'
                cmap3 = plt.cm.get_cmap(cmap_name3,len(vrc))
                norm3 = colors.BoundaryNorm(vrc, cmap3.N)
                
            #### calculate the mean change in carbon ####
            weighted_data = x_change*areas # area weight the variable data
            avr_change = np.nanmean(weighted_data) #average across the globe
            
    ### calculate sum across globe ###
    abs_mean = np.nanmean(abss, axis=0) #average across models where abss is the change between first and last decade of run
    #initialize zero arrays 
    ens_change = np.zeros(len(models))
    end1 = np.zeros(len(models))
    end2 = np.zeros(len(models))
    #remove filler values from futrs arrays 
    futrs2[futrs2>1e36]  = np.nan
    futrs1[futrs1>1e36]  = np.nan
    
    ### calculate average change in each variable to be put in panel title ###
    if var == 'NPP' or var =='cLand':
        for i in range (len(models)):
            end1[i] = np.nansum((futrs1[i,:,:]*areas))/land_area
            end1_mean = np.mean(end1)
            end2[i] = np.nansum((futrs2[i,:,:]*areas))/land_area
            end2_mean = np.mean(end2)

    else: 
        for i in range (len(models)):            
            end1[i] = np.nansum((futrs1[i,:,:]*areas))/np.sum(areas)
            end1_mean = np.mean(end1)
            end2[i] = np.nansum((futrs2[i,:,:]*areas))/np.sum(areas)
            end2_mean = np.mean(end2)
             
    ens_avr_change2 = end2_mean-end1_mean #average change to be put in titles    

    # ens_rel_mean = (ens_avr_change2/end1_mean)*100
    # print(str(ens_rel_mean)+' '+str(var)) 

    ##### Create stippling #####
    ens_std = np.std(abss, axis=0)
    stps = abs(ens_std/abs_mean)
    
    ### plot maps ###   
    axes[k].text(.02, 1.07, letters[k], ha='left', va='top', transform=axes[k].transAxes, fontsize='large') #label panels 
    
    axes[k].set_extent([lonmin ,lonmax, latmin, latmax], crs=ccrs.PlateCarree())   
    axes[k].add_feature(cfeature.COASTLINE, zorder=10)
    
    axes[k].set_xticks(np.arange(-80,-30,lonstep),crs=ccrs.PlateCarree())
    axes[k].set_xticklabels(np.arange(-80,-30,lonstep))
    axes[k].set_yticks(np.arange(-20,20,latstep),crs=ccrs.PlateCarree())
    axes[k].set_yticklabels(np.arange(-20,20,latstep))
    
    lon_formatter = cticker.LongitudeFormatter()
    lat_formatter = cticker.LatitudeFormatter()
    
    axes[k].xaxis.set_major_formatter(lon_formatter)
    axes[k].yaxis.set_major_formatter(lat_formatter)
    axes[k].grid(linewidth=0.5, color='black', linestyle='--') #create grid  
    
    #plot image of data
    im = axes[k].imshow(abs_mean,cmap=cmap,norm=norm,transform=ccrs.PlateCarree(),extent=[lonmin ,lonmax, latmax, latmin])
    
    #plot stippling 
    stip2 = axes[k].contourf(longitude, latitude,land_use_mask, tranform=ccrs.PlateCarree(), levels=[0, 0.99, 2], 
                              colors='none', hatches=[None,'xxx'])           
    stip = axes[k].contourf(longitude, latitude, stps, transform=ccrs.PlateCarree(), levels = [0, 1, 100000],
                            colors='none', hatches=[None, '..'])
    
    #set title 
    axes[k].set_title('Change in '+var_names[k]+' ('+f'{ens_avr_change2:.3}'+' '+units+')') # plot figure title
    
    #set colourbar
    if var=='tas':
        cbar = fig.colorbar(im, ticks=[-5,-4,-3,-2,-0.5,0.5,2,3,4,5] ,spacing="proportional", pad=0.08, orientation='horizontal', extend='both')
    else:
        cbar = fig.colorbar(im, spacing="proportional", pad=0.08, orientation='horizontal', extend='both')
    
    cbar.set_label('Change in '+var_names[k]+' ('+units+')')

#save figure 
fig.savefig('C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/Figures/Geoeng_paper/'+experiment1+'_vsG6sulfur_map_Amazon_V1.png', bbox_inches='tight') 

  

       
