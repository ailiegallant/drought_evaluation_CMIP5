#-----------------------------------------------------------------------------
'''
NAME
    NC_SSTBOUNDARY
PURPOSE
    Generates an SST boundary climatology file using a series of input years that 
    are user-defined. SST boundary generated using the HadISST data set.
PROGRAMMER(S)
    Ailie Gallant
REVISION HISTORY
    20160629 -- Script written by Ailie Gallant
REFERENCES
    netcdf4-python -- http://code.google.com/p/netcdf4-python/
'''

def nc_sstboundary():
	from netCDF4 import date2num
	from netCDF4 import num2date
	import numpy as np
	from netcdf_tools import ncwrite_climgrid
	from netcdf_tools import ncextractall

	#Input file (this file goes from January 1870 to April 2016 as is)
	file = '/Users/ailieg/Data/HadISST_sst.nc'
	
	nc = ncextractall(file)
	sst = nc['sst']
	lon = nc['longitude']
	nlon = lon.size
	lat = nc['latitude']
	nlat = lat.size
	time = nc['time']
	miss = nc['sst_missing_value']
	units = nc['sst_units']
	
	#convert time units to actual date
	time_u = nc['time_units']
	if 'time_calendar' in nc.keys(): 
	     cal = nc['time_calendar']
	     time = num2date(time,units = time_u, calendar=cal)
	else: time = num2date(time,units = time_u)
	
	#extract years and months from the datetime array
	ntime = time.size
	year = np.zeros(ntime)
	month = np.zeros(ntime)
	i = 0
	while (i < ntime):
		year[i] = time[i].year
		month[i] = time[i].month
		i=i+1
		
	#Extract data from 1950 to 2015 only
	sst = sst[(year > 1949) & (year < 2016), :,:]
	time = time[(year > 1949) & (year < 2016)]
	ntime = time.size
	
	#Reshape the array to determine climatology over the whole period
	nyear = ntime/12
	sst = sst.reshape([nyear,12,nlat,nlon])
	time = time.reshape([nyear, 12])
	time = time[32,:]
	
	sst = np.mean(sst, axis=0)
	
	#Write the output file
	descrip = 'Monthly climatology of HadISST SSTs computed from 1950-2015'
	long_name = 'sst'
	missing = miss
	climunits = units
	time = date2num(time,units = time_u, calendar=cal)
	print(sst.shape,time.size,lon.size,lat.size)
	filename = '/Users/ailieg/Data/IPO_ACCESS/HadISST_climatol_1950_2015.nc'
	
	print("Writing netCDF file...")
	ncw = ncwrite_climgrid(filename, sst, 'sst', descrip, long_name, missing, climunits, time, lon, lat, time_u, cal)
	print("NetCDF file written")
	
	return sst, lat, lon

#-----------------------------------------------------------------------------
'''
NAME
    NC_IPOPHASE
PURPOSE
    Generates an SST boundary climatology file using a series of input years that are
    defined as being above or below a threshold of a particular IPO value
PROGRAMMER(S)
    Ailie Gallant
REVISION HISTORY
    20160629 -- Script written by Ailie Gallant
REFERENCES
    netcdf4-python -- http://code.google.com/p/netcdf4-python/
'''

def nc_ipophase():
	from netCDF4 import date2num
	from netCDF4 import num2date
	import numpy as np
	from netcdf_tools import ncwrite_climgrid
	from netcdf_tools import ncextractall

	#Input file (this file goes from January 1870 to April 2016 as is)
	file = '/Users/ailieg/Data/HadISST_sst.nc'
	
	nc = ncextractall(file)
	sst = nc['sst']
	lon = nc['longitude']
	nlon = lon.size
	lat = nc['latitude']
	nlat = lat.size
	time = nc['time']
	miss = nc['sst_missing_value']
	units = nc['sst_units']
	
	#convert time units to actual date
	time_u = nc['time_units']
	if 'time_calendar' in nc.keys(): 
	     cal = nc['time_calendar']
	     time = num2date(time,units = time_u, calendar=cal)
	else: time = num2date(time,units = time_u)
	
	#extract years and months from the datetime array
	ntime = time.size
	year = np.zeros(ntime)
	month = np.zeros(ntime)
	i = 0
	while (i < ntime):
		year[i] = time[i].year
		month[i] = time[i].month
		i=i+1
		
	#Extract data from 1950 to 2015 only
	sst = sst[(year > 1976) & (year < 2000), :,:]
	time = time[(year > 1976) & (year < 2000)]
	ntime = time.size
	
	#Reshape the array to determine climatology over the whole period
	nyear = ntime/12
	sst = sst.reshape([nyear,12,nlat,nlon])
	time = time.reshape([nyear, 12])
	mid = int(nyear/2)
	time = time[mid,:]
	
	sst = np.mean(sst, axis=0)
	
	#Write the output file
	descrip = 'Monthly climatology of HadISST SSTs from IPO positive years computed as 1977-1999'
	long_name = 'sst'
	missing = miss
	climunits = units
	time = date2num(time,units = time_u, calendar=cal)
	print(sst.shape,time.size,lon.size,lat.size)
	filename = '/Users/ailieg/Data/IPO_ACCESS/HadISST_IPOpos_1977_1999.nc'
	
	print("Writing netCDF file...")
	ncw = ncwrite_climgrid(filename, sst, 'sst', descrip, long_name, missing, climunits, time, lon, lat, time_u, cal)
	print("NetCDF file written")
	
	return sst, lat, lon

#-----------------------------------------------------------------------------
'''
NAME
    NC_IPOPHASE
PURPOSE
    Generates an SST boundary climatology file using a series of input years that are
    defined as being above or below a threshold of a particular IPO value
PROGRAMMER(S)
    Ailie Gallant
REVISION HISTORY
    20160629 -- Script written by Ailie Gallant
REFERENCES
    netcdf4-python -- http://code.google.com/p/netcdf4-python/
'''

	fill_missing(z, y, x)
