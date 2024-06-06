    #- coding: utf-8 -*-
"""
Created on Mon Dec 16 09:01:27 2019

@author: pdlr201

Script to interpolate annual 1pctCO2 CMIP6 data onto 1 degree world grid
"""

import iris
import numpy as np
import numpy.ma as ma
from iris.cube import Cube

def shiftlon(lon,lon_0):
    """returns original sequence of longitudes (in degrees) recentered
    in the interval [lon_0-180,lon_0+180]"""
    lon_shift = np.asarray(lon)
    lon_shift = np.where(lon_shift > lon_0+180, lon_shift-360 ,lon_shift)
    lon_shift = np.where(lon_shift < lon_0-180, lon_shift+360 ,lon_shift)
    itemindex = len(lon)-np.where(lon_shift[0:-1]-lon_shift[1:] >= 180)[0]
    return (np.roll(lon_shift,itemindex-1), itemindex-1)


################# Specify variable of interest ######################
var = 'cLand' 
############# Specify temporal resolution of dataset #################
temp_res = 'Emon'
############# Specify experiment of interest ##############
experiment = 'G6sulfur'

################# Specify models here ################################## 
models = ['CESM2-WACCM', 'CNRM-ESM2-1','IPSL-CM6A-LR', 'MPI-ESM1-2-LR', 'UKESM1-0-LL']


for m in range (np.size(models)):
    model = models[m]    
   
    # if model == 'NorESM2-LM':
    #     grid_label = '_gn'
    #     date_ranges = ['014101-014912']
    #     for k in np.arange(25):
    #         date_ranges.extend([(str(150+k*10)+'01').zfill(6)+'-'+(str(159+(k)*10)+'12').zfill(6)])
    #     total_years = 260
    #     date_range2 = '141-399'
    #     variant_id = 'r1i1p1f1'

    if model == 'UKESM1-0-LL':
        grid_label = '_gn'
        date_ranges = ['202001-204912', '205001-210012'] #list of date ranges in each file
        # date_ranges = ['20200101-20491230', '20500101-21001230']
        total_years = 81 #total number of years in model dataset
        date_range2 = '2020-2100' #string of overall date range
        variant_id = 'r1i1p1f2'  #varaint ID of dataset
    
    elif model == 'CNRM-ESM2-1':
        grid_label = '_gr'
        date_ranges = ['201501-210012']
        # date_ranges = ['20150101-21001231']
        total_years = 86
        date_range2 = '2015-2100'
        variant_id = 'r1i1p1f2'
    
    elif model == 'CESM2-WACCM':
        if experiment == 'G6solar':
            date_ranges = ['201901-206812', '206901-210012']
            total_years = 82
            variant_id = 'r1i1p1f1'
        elif experiment == 'G6sulfur':
            date_ranges = ['202001-206912', '207001-210012']
            total_years = 81
            variant_id = 'r1i1p1f2'
        grid_label = '_gn'
        
        date_range2 = '2020-2100'
        
    
    elif model == 'IPSL-CM6A-LR':
        grid_label = '_gr'
        date_ranges = ['202001-210012']
        total_years = 81
        date_range2 = '2020-2100'
        variant_id = 'r1i1p1f1'
    
    elif model == 'MPI-ESM1-2-LR':
        grid_label = '_gn'
        date_ranges = ['201501-203412', '203501-205412', '205501-207412', '207501-209412', '209501-209912']
        total_years = 85
        date_range2 = '2015-2099'
        variant_id = 'r1i1p1f1'
       
    elif model == 'MPI-ESM1-2-HR':
        grid_label = '_gn'
        date_ranges = []
        for k in np.arange(16): # don't forget to change for the different date ranges
            # date_ranges.extend([(str(2015+k*5)+'01').zfill(6)+'-'+(str(2019+(k)*5)+'12').zfill(6)])
            date_ranges.extend([(str(2020+k*5)+'01').zfill(6)+'-'+(str(2024+(k)*5)+'12').zfill(6)])
        total_years = 80 
        date_range2 = '2020-2099' 
        variant_id = 'r1i1p1f1'
        
    ######################## Specify path to data ##########################
    # (I have 2 directories a raw data directory and a processed data directory)
    path = 'C:/Users/ip294/OneDrive - University of Exeter/PhD_Tipping_Points/'+var+'/'+experiment
    
    ######################## Specify directory of raw data ##################
    raw_data = '/Original_data/'
    
    ######################## Specify directory of processed data ##################
    processed_data = '/Processed_data/'
    
    ############ Specify directory within processed data directory for the region interpolated #################
    region = 'World'
    
    # Define latitude and longitude (uses 1 degree by 1 degree grid) 
    my_lats=['latitude',np.arange(-90,91,step=1)]
    my_lons=['longitude',np.arange(-180,181,step=1)]
    
    # Months in the year (or days in the year)
    freq = 12
    
    # Initialise counters
    count = 0
    sub_total_years = 0
    
    # Loop through data files
    for i in range(len(date_ranges)):
        
        # File name of data
        fname = path+raw_data+var+'_'+temp_res+'_'+model+'_'+experiment+'_'+variant_id+grid_label+'_'+date_ranges[i]+'.nc'
    
        # Load file contents into an iris cube
        x = iris.load_cube(fname)
        
        # Extract data values
        data = x.data
    
        # Determine lengths of data dimensions
        nt = int(len(data[:,0,0]))
        ny = int(len(data[0,:,0]))
        nx = int(len(data[0,0,:]))
        
        years = int((nt+1)/freq)
        
        # On first loop through create a new empty cube with original coordinates
        if count == 0:
        
            new_coord = iris.coords.DimCoord(range(total_years), long_name='time', units=1)
            coord1 = x.coord('latitude')
            
            if model == 'CESM2' or model == 'CESM2-WACCM':
                coord2 = iris.coords.DimCoord(np.arange(0,360,step=1.25), standard_name='longitude', units='degrees', var_name='lon', attributes={'title': 'Longitude', 'type': 'double', 'valid_max': 180.0, 'valid_min': -180.0}, circular=True)
            else:
                coord2 = x.coord('longitude')
            
            cube = Cube(ma.zeros((total_years,ny,nx),np.float32),dim_coords_and_dims=[(new_coord,0),(coord1,1),(coord2,2)])
        
        count2 = 0
        
        # Stack all data into 1 cube 
        for j in range(sub_total_years,sub_total_years+years):
    
            cube.data[j,:,:] = np.mean(data[count2*freq:(count2+1)*freq,:,:], axis=0)
            
            count2 = count2 + 1
      
        sub_total_years = sub_total_years + years
        count = count+1
        
    # Centre coordinates about longitude 0
    cube.coord('longitude').points, rollindex = shiftlon(cube.coord('longitude').points, 0)
    
    # Roll data to be consistent with longitude 0
    cube.data = np.roll(cube.data, rollindex, axis=2)
    
    # Interpolate data onto new coordinate system
    region_cube = cube.interpolate([my_lats, my_lons],iris.analysis.Linear(extrapolation_mode='nan'))
    
    # Ensure new cube has correct name
    region_cube.rename(var)
    
    # New file name for data
    outname = path+processed_data+region+'/'+var+'_'+model+'_'+experiment+'_'+variant_id+'_'+date_range2+'_'+region+'.nc'
    
    # Save data to file name
    iris.save(region_cube, outname, unlimited_dimensions=[])
