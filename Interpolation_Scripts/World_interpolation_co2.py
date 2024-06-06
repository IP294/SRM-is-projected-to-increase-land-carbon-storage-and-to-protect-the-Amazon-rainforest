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

# You shouldn't need to change the below entries at least for the 3 models you
# are currently using

experiment = 'historical'

################# Specify variable of interest ######################
var = 'co2'

################# Specify models here ################################## 
models = ['CESM2-WACCM']

for m in range (np.size(models)):
    model = models[m]  
    
    print(model + ' start')
    if model == 'NorESM2-LM':
        grid_label = '_gn'
        date_ranges = ['014101-014912']
        for k in np.arange(25):
            date_ranges.extend([(str(150+k*10)+'01').zfill(6)+'-'+(str(159+(k)*10)+'12').zfill(6)])
        total_years = 260
        date_range2 = '141-399'
        variant_id = 'r1i1p1f1'
        
    elif model == 'CESM2-WACCM':
        grid_label = '_gn'
        date_ranges = ['185001-201412']
        # date_ranges = ['201501-210012', '210101-215012', '215101-220012', '220101-225012', '225101-229912']
        # for k in np.arange(5):
        #     date_ranges.extend([(str(2050+k*100)+'01').zfill(6)+'-'+(str(2149+(k)*100)+'12').zfill(6)])
        # date_ranges.extend(['255001-263912'])
        total_years = 250
        date_range2 = '1850-2014'
        variant_id = 'r1i1p1f1'

    elif model == 'UKESM1-0-LL':
        grid_label = '_gn'
        date_ranges = ['199001-204912']
        for k in np.arange(5):
            date_ranges.extend([(str(2050+k*100)+'01').zfill(6)+'-'+(str(2149+(k)*100)+'12').zfill(6)])
        date_ranges.extend(['255001-263912'])
        total_years = 650
        date_range2 = '1999-2639'
        variant_id = 'r1i1p1f2'
    
    elif model == 'CanESM5':
        grid_label = '_gn'
        date_ranges1 = ['199101-219012', '219101-229012']
        date_ranges2 = ['185001-200012']
        total_years1 = 300
        total_years2 = 150
        date_range1 = '1991-2290'
        date_range2 = '1850-2000'
        variant_id = 'r1i1p2f1'
    
    elif model == 'ACCESS-ESM1-5':
        grid_label = '_gn'
        date_ranges = ['024101-038012', '038101-088012', '088101-100012']
        total_years = 760
        date_range2 = '241-1000'
        variant_id = 'r1i1p1f1'
    
    elif model == 'MIROC-ES2L':
        grid_label = '_gn'
        date_ranges = ['198901-218812', '218901-238812', '238901-249012']
        total_years = 502
        date_range2 = '1989-2388'
        variant_id = 'r1i1p1f2'
    
    elif model == 'CESM2':
        grid_label = '_gn'
        date_ranges = ['000101-004912']
        for k in range(3):
            date_ranges.extend([(str(50+k*50)).zfill(4)+'01-'+str(99+k*50).zfill(4)+'12'])
        date_ranges.extend(['020001-020012'])
        total_years = 250
        date_range2 = '1-199'
        variant_id = 'r1i1p1f1'
    
    elif model == 'GFDL-ESM4':
        grid_label = '_gr1'
        date_ranges = ['014101-024012', '024101-034012']
        # date_ranges2 = ['000101-010012', '010101-015012']
        total_years = 200
        # total_years2 = 150
        date_range2 = '140-340'
        # date_range2 = '01-150'
        variant_id = 'r1i1p1f1'
    
    elif model == 'CNRM-ESM2-1':
        grid_label = '_gr'
        # date_ranges = ['201501-210012']
        date_ranges = ['185001-199912']
        total_years = 150
        # date_range2 = '2015-2100'
        date_range2 = '1850-1999'
        variant_id = 'r1i1p1f2'
    
    elif model == 'BCC-CSM2-MR':
        grid_label = '_gn'
        date_ranges = ['201501-205412', '205501-209412', '209501-210012']
        total_years = 86
        date_range2 = '2015-2100'
        # date_range2 = '1850-1999'
        variant_id = 'r1i1p1f1'
   
    if model == 'NorESM2-LM':
        temp_res = 'AERmon'
    else:
        temp_res = 'Amon'
    
    # if experiment =='1pctCO2':
    #     date_ranges = date_ranges2
    #     total_years = total_years2
    #     date_range = date_range2
    # else:
    #     date_ranges = date_ranges1
    #     total_years = total_years1
    #     date_range = date_range1
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
    
    # Months in the year
    freq = 12
    
    # Initialise counters
    count = 0
    sub_total_years = 0
    
    # Loop through data files
    for i in range(len(date_ranges)):
        
        # File name of data
        fname = path+raw_data+var+'_'+temp_res+'_'+model+'_'+experiment+'_'+variant_id+grid_label+'_'+date_ranges[i]+'.nc'
    
        # Load file contents into an iris cube
        if model =='NorESM2-LM':
            x = iris.load_raw(fname)
            x = x[0]
        else:
            x = iris.load_cube(fname)
        
        # Extract data values
        data = x.data
    #%%
        # Determine lengths of data dimensions
        nt = int(len(data[:,0,0,0]))
        na = int(len(data[0,:,0,0]))
        ny = int(len(data[0,0,:,0]))
        nx = int(len(data[0,0,0,:]))
        
        years = int((nt+1)/freq)
        
        # On first loop through create a new empty cube with original coordinates
        if count == 0:
        
            new_coord = iris.coords.DimCoord(range(total_years), long_name='time', units=1)
            coord1 = x.coord('latitude')
            if model == 'NorESM2-LM':
                coord3 = x.coord('atmosphere_hybrid_sigma_pressure_coordinate')
            else:
                coord3 = x.coord('air_pressure')
            if model == 'CESM2' or model == 'CESM2-WACCM':
                coord2 = iris.coords.DimCoord(np.arange(0,360,step=1.25), standard_name='longitude', units='degrees', var_name='lon', attributes={'title': 'Longitude', 'type': 'double', 'valid_max': 180.0, 'valid_min': -180.0}, circular=True)
            else:
                coord2 = x.coord('longitude')
            
            cube = Cube(ma.zeros((total_years,na,ny,nx),np.float32),dim_coords_and_dims=[(new_coord,0),(coord3,1),(coord1,2),(coord2,3)])
        
        count2 = 0
        
        # Stack all data into 1 cube 
        for j in range(sub_total_years,sub_total_years+years):
    
            cube.data[j,:,:,:] = np.mean(data[count2*freq:(count2+1)*freq,:,:], axis=0)
            
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
