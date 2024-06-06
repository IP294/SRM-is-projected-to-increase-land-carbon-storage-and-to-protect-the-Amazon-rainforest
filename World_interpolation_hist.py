#- coding: utf-8 -*-
"""
Created on Mon Dec 16 09:01:27 2019

@author: pdlr201

Script to interpolate annual historical CMIP6 data onto 1 degree world grid
"""

### import libraries ###
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

#### Specify temoral resolution ###
temp_res = 'Emon'

#### Specify experiment of interest ###
experiment = 'historical'

####### Specify variable of interest ######
var = 'cLand' 

################# Specify model here ################################## 
models = ['CESM2-WACCM', 'CNRM-ESM2-1', 'IPSL-CM6A-LR', 'MPI-ESM1-2-LR', 'UKESM1-0-LL']

for m in range (np.size(models)): #loop through models 
    model = models[m] #select model 
       
    ##### define details associated with each model #####
    if model == 'EC-Earth3-Veg':
        grid_label = '_gr' #specify grid labell
        date_ranges = [] #create list of date ranges for each file in dataset
        for k in np.arange(151):
            date_ranges.extend([str(1850+k)+'01-'+str(1850+k)+'12'])
        total_years = 151 #total number of model years in dataset
        date_range2 = '1850-2000' #total range of years in dataset
        variant_id = 'r1i1p1f1' #specify variant ID
    
    elif model == 'CNRM-ESM2-1':
        grid_label = '_gr'
        date_ranges = ['185001-201412']
        total_years = 165
        date_range2 = '1850-2014'
        variant_id = 'r1i1p1f2'    
        
    elif model == 'MPI-ESM1-2-LR':
        grid_label = '_gn'
        date_ranges = []
        for k in np.arange(8):
            date_ranges.extend([str(1850+k*20)+'01-'+str(1869+(k)*20)+'12'])
        date_ranges.extend(['201001-201412'])
        total_years = 165
        date_range2 = '1850-2014'
        variant_id = 'r1i1p1f1'
    
    elif model == 'MPI-ESM1-2-HR':
        grid_label = '_gn'
        date_ranges = []
        for k in np.arange(33):
            date_ranges.extend([(str(1850+k*5)+'01').zfill(6)+'-'+(str(1854+(k)*5)+'12').zfill(6)])
        total_years = 165
        date_range2 = '1850-2014'
        variant_id = 'r1i1p1f1'
            
    elif model == 'UKESM1-0-LL':
        grid_label = '_gn'
        date_ranges = ['185001-194912', '195001-201412']
        total_years = 165
        date_range2 = '1850-2014'
        variant_id = 'r1i1p1f2'
        
    elif model == 'CESM2-WACCM':
        grid_label = '_gn'
        date_ranges = ['185001-201412']
        total_years = 165
        date_range2 = '1850-2014'
        variant_id = 'r1i1p1f1'
    
    elif model == 'IPSL-CM6A-LR':
        grid_label = '_gr'
        date_ranges = ['185001-201412']
        total_years = 165
        date_range2 = '1850-2014'
        variant_id = 'r1i1p1f1'
        
    elif model == 'GFDL-ESM4':
        
        grid_label = '_gr1'
        date_ranges = ['000101-010012', '010101-015012']
        total_years = 150
        date_range2 = '1850-1999'
        variant_id = 'r1i1p1f1'
    
    elif model == 'NorCPM1':
        grid_label = '_gn'
        date_ranges = ['000102-016412']
        total_years = 164
        date_range2 = '1850-2013'
        variant_id = 'r1i1p1f1'
    
    elif model == 'TaiESM1':
        grid_label = '_gn'
        date_ranges = ['000102-015012']
        total_years = 150
        date_range2 = '1850-2000'
        variant_id = 'r1i1p1f1'
        
    elif model == 'SAM0-UNICON':
        grid_label = '_gn'
        date_ranges = []
        for k in np.arange(15):
            date_ranges.extend([str(1850+(k*10))+'01-'+str(1859+(k*10))+'12'])
        total_years = 150
        date_range2 = '1850-1999'
        variant_id = 'r1i1p1f1'
      
    elif model == 'CanESM5':
        grid_label = '_gn'
        date_ranges = []
        for k in np.arange(5):
            date_ranges.extend([str(5201+(k*200))+'01-'+str(5400+(k*200))+'12'])
        total_years = 1999
        date_range2 = '0000-0000'
        
    elif model == 'BCC-ESM1':
        grid_label = '_gn'
        date_ranges = ['185001-230012']
        total_years = 451
        date_range2 = '1850-2300'
        variant_id = 'r1i1p1f1'
        
    elif model == 'CESM2':
        grid_label = '_gn'
        date_ranges = ['000101-009912']
        for k in np.arange(9):
            date_ranges.extend(['0'+str(100+(k*100))+'01-0'+str(199+(k*100))+'12'])
        total_years = 1200
        date_range2 = '0000-0000'
        variant_id = 'r1i1p1f1'
    
    elif model == 'AWI-ESM-1-1-LR':
        grid_label = '_gn'
        date_ranges = []
        for k in np.arange(100):
            date_ranges.extend([str(1855+(k*1))+'01-'+str(1855+(k*1))+'12'])
        total_years = 100
        date_range2 = '1855-1954'
        variant_id = 'r1i1p1f1'
        
    elif model == 'ACCESS-ESM1-5':
        grid_label = '_gn'
        date_ranges = ['010101-060012', '060101-100012','100101-110012']
        total_years = 1000
        date_range2 = '100-1100'
        variant_id = 'r1i1p1f1'
    
    elif model == 'BCC-CSM2-MR':
        grid_label = '_gn'
        date_ranges = ['185001-244912']
        total_years = 600
        date_range2 = '1850-2449'
        variant_id = 'r1i1p1f1'
    
    elif model == 'CMCC-CM2-SR5':
        grid_label = '_gn'
        date_ranges = ['185001-209912','210001-234912']
        total_years = 500
        date_range2 = '1850-2349'
        variant_id = 'r1i1p1f1'
    
    elif model == 'CMCC-ESM2':
        grid_label = '_gn'
        date_ranges = ['185001-209912','210001-234912']
        total_years = 500
        date_range2 = '1850-2349'
        variant_id = 'r1i1p1f1'
    
    elif model == 'MPI-ESM-1-2-HAM':
        grid_label = '_gn'
        date_ranges = []
        for k in np.arange(50):
            date_ranges.extend([str(1850+(k*20))+'01-'+str(1869+(k*20))+'12'])
        total_years = 1000
        date_range2 = '1850-2849'
        variant_id = 'r1i1p1f1'
    
    elif model == 'NorESM2-LM':
        grid_label = '_gn'
        date_ranges = []
        for k in np.arange(20):
            date_ranges.extend([str(1600+(k*10))+'01-'+str(1609+(k*10))+'12'])
        date_ranges.extend(['180001-180012'])
        for k in np.arange(30):
            date_ranges.extend([str(1801+(k*10))+'01-'+str(1810+(k*10))+'12'])
        total_years = 501
        date_range2 = '1600-2100'
        variant_id = 'r1i1p1f1'
    
    elif model == 'MIROC-ES2L':
        grid_label = '_gn'
        date_ranges = ['185001-204912', '205001-224912', '225001-234912']
        total_years = 500
        date_range2 = '1850-2350'
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
