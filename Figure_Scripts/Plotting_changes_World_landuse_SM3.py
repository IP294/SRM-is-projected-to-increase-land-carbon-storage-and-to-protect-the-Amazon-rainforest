# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 10:14:22 2020

@author: pdlr201

Plotting land use changes from LUH2 comparing ssp585 and ssp245
"""
### load in libraries ###
from pathlib import Path
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.mpl.ticker as cticker
import matplotlib.colors as mcolors

experiments = ['ssp585', 'ssp245'] #experiments for comparison

################# Specify variables of interest ######################
var = 'primf' 
var2 = 'secdf'

#initialize ensemble zero array 
lands = np.zeros((2,86,181,361))

############## Specify units of variable ###########################
units = 'kgC m$^{-2}$'

############ Specify name of directory of data you want to plot #################
region = 'World'

if region == 'World':
    lonmin = -180-0.5
    lonmax = 180+0.5
    latmin = -90-0.5
    latmax = 90+0.5
    lonstep = 60
    latstep = 20
    figsizes = (22, 20)
elif region == 'Amazon':
    lonmin = -83-0.5
    lonmax = -35+0.5
    latmin = -20-0.5
    latmax = 13+0.5
    lonstep = 10
    latstep = 10
    figsizes = (6, 5.75)
    
############# Plotting the figures, specify where you want the figures saved ############
proj = ccrs.PlateCarree()
 
#### Initialize figure ####
fig1, ax1 = plt.subplots(1, subplot_kw=dict(projection=proj), figsize=(11,10), sharex=True, sharey=True)

################## Specify path to processed data directory #####################
path = 'C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/land_use/Processed_data/'

for i in range (len(experiments)): #loop through experimentd 
    experiment = experiments[i] #select experiment 
    
    # File name of interpolated data for each variable 
    fname = path+region+'/'+var+'_'+experiment+'_states_landuse.nc'
    fname = Path(fname)
    
    fname2 = path+region+'/'+var2+'_'+experiment+'_states_landuse.nc'
    fname2 = Path(fname2)
    
    # Load in both datasets 
    f1 = nc.Dataset(fname,'r')
    f2 = nc.Dataset(fname2,'r')
    
    # Extract longitude and latitude
    longitude = f1.variables['longitude'][:]
    latitude = f1.variables['latitude'][:]
    
    x1 = f1.variables[var][:,:,:]
    x2 = f2.variables[var2][:,:,:]
    
    # Close datasets
    f1.close()
    f2.close()
    
    x = x1+x2 #add two datasets together to get total prescribed changes 
    lands[i,:,:,:] = x #save into ensemble array

# Obtain dimension sizes of data 
nt = int(len(x1[:,0,0]))
ny = int(len(x1[0,:,0]))
nx = int(len(x1[0,0,:]))

########### Specify path to land mask ######################
fname2 = 'C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/Masks/'+region+'.nc'
fname2 = Path(fname2)

# Load in land mask and remove any ocean values
f2 = nc.Dataset(fname2,'r')
sftlf = f2.variables['sftlf'][:]
sftlf_mask = np.broadcast_to(sftlf, (nt, ny,nx))
x[np.where(sftlf_mask <= 0)] = np.nan  
x[np.where(~np.isfinite(sftlf_mask))] = np.nan #removing ocean values 

f2.close() #close land mask dataset


#find mean of last ten years for each eperiment 
x_futr1 = np.nanmean(lands[0][nt-10:nt,:,:],axis=0) #ssp585
x_futr2 = np.nanmean(lands[1][nt-10:nt,:,:],axis=0) #ssp245

x_tot_change = (x_futr2 - x_futr1) #find difference in prescribed change between experiments

### define colour map ###
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



### plotting map ###
ax1.set_extent([lonmin ,lonmax, latmin, latmax], crs=ccrs.PlateCarree())   
ax1.add_feature(cfeature.COASTLINE, zorder=10) #add map image 

#set latitude and longitude ticks 
ax1.set_xticks(np.arange(-180,210,lonstep),crs=ccrs.PlateCarree())
ax1.set_xticklabels(np.arange(-180,210,lonstep))
ax1.set_yticks(np.arange(-80,80,latstep),crs=ccrs.PlateCarree())
ax1.set_yticklabels(np.arange(-80,80,latstep))

lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()

ax1.xaxis.set_major_formatter(lon_formatter)
ax1.yaxis.set_major_formatter(lat_formatter)

ax1.grid(linewidth=0.5, color='black', linestyle='--') #add gridlines

#plot data as image
im = ax1.imshow(x_tot_change,cmap=cmap,norm=norm,transform=ccrs.PlateCarree(),extent=[lonmin ,lonmax, latmax, latmin])

#plot colourbar   
cbar = fig1.colorbar(im, spacing="proportional", pad=0.08, orientation='horizontal', extend='both')
cbar.set_label('Prescribed forest fraction change')
