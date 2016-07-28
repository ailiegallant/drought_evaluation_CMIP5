#THIS MODULE CONTAINS TOOLS FOR READING AND WRITING NETCDF FILES
#EACH SUBROUTINE IS SEPARATED BY A LINE --------
#-----------------------------------------------------------------------------

'''
NAME
    NCEXTRACTALL
PURPOSE
    Read all variables in a NetCDF file 
PROGRAMMER(S)
    Ailie Gallant
REVISION HISTORY
    20160413 -- Generation of the script
REFERENCES
    netcdf4-python -- http://code.google.com/p/netcdf4-python/
'''
#import datetime as dt  # Python standard library datetime  module
import numpy as np
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/

def ncextractall(filename):
    '''
    ncdump outputs variables.
    ncdump requires a valid instance of Dataset.

    Parameters
    ----------
    nc_fid : netCDF4.Dataset
        A netCDF4 dateset object

    Returns
    -------
    nc_vars : list
        A Python list of the NetCDF file variables
    '''
    
    #Read the netcdf as a read-only
    ncid = Dataset(filename,mode='r')
    
    #Read the variables and variable names
    #Variable information will be read as an OrderedDict and will
    #require individual items to be named
    varnc = ncid.variables
    
    #Determine how many items in the dictionary
    varitems = list(varnc.items())
    nvar = len(varitems)
    out = {}

    #Loop over names and make into array
    for i in range(0,nvar):
        thisvar = (varitems[i])[0]
        
        dictstring = thisvar
        attlist = ncid.variables[thisvar]
        out[dictstring] = ncid.variables[thisvar][:]
        
        if hasattr(attlist, 'units'): \
           out[dictstring+'_units'] = ncid.variables[thisvar].units
            
        if hasattr(attlist, 'long_name'): \
            out[dictstring+'_name'] = ncid.variables[thisvar].long_name
            
        if hasattr(attlist, 'missing_value'): \
            out[dictstring+'_missing_value'] = ncid.variables[thisvar].missing_value
        
        if hasattr(attlist, 'add_offset'): \
            out[dictstring+'_add_offset'] = ncid.variables[thisvar].add_offset

        if hasattr(attlist, 'scale_factor'): \
            out[dictstring+'_scale_factor'] = ncid.variables[thisvar].scale_factor
            
        if hasattr(attlist, 'calendar'): \
           out[dictstring+'_calendar'] = ncid.variables[thisvar].calendar
        

    ncid.close()
    return out
        

#-----------------------------------------------------------------------------

'''
NAME
    NCWRITE_CLIMGRID
PURPOSE
    To write a NetCDF file of 
PROGRAMMER(S)
    Ailie Gallant
REVISION HISTORY
    20160414 -- Script written by Ailie Gallant
REFERENCES
    netcdf4-python -- http://code.google.com/p/netcdf4-python/
'''

def ncwrite_climgrid(filename, climdata, climname, descrip, long_name, missing, climunits, 
                       time, lon, lat, time_units, time_cal):

     '''
     ncwrite_climgrid(filename, climdata, descrip, long_name):

     Must be an input array of climate data of the form climdata(time, lat, lon).
     Input vectors of time, lon and lat must be provided. Time must be in the format to
     write (not datetime format).
   
     time_units - must be a string in the format of <time units> since <reference time>. 
     For example, "days since 1800-1-1 00:00:0.0"
     '''

     import numpy as np
     from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/

     #Create NetCDF file to write
     w_nc = Dataset(filename, 'w', format='NETCDF4')
     
     #File description
     w_nc.description = "%s" % (descrip)
     
     #File dimensions for TIME
     w_nc.createDimension('time', len(time))
     w_nc_time = w_nc.createVariable('time', time.dtype, ('time',))

     w_nc_time.setncatts({'long_name': 'time',\
                    'units': time_units, 'calendar': time_cal})
                    
     # Assign the dimension data to the new NetCDF file.
     w_nc.variables['time'][:] = time
     
     #File dimensions for LAT
     w_nc.createDimension('lat', len(lat))
     w_nc_lat = w_nc.createVariable('lat', lat.dtype,('lat',))

     w_nc_lat.setncatts({'long_name': 'Latitude',\
                    'units': 'Degrees North'})
                    
     # Assign the dimension data to the new NetCDF file.
     w_nc.variables['lat'][:] = lat
     
     #File dimensions for LON
     w_nc.createDimension('lon', len(lon))
     w_nc_lon = w_nc.createVariable('lon', lon.dtype,('lon',))

     w_nc_lon.setncatts({'long_name': 'Longitude',\
                    'units': 'Degrees East'})
                    
     # Assign the dimension data to the new NetCDF file.
     w_nc.variables['lon'][:] = lon

     # Assign the climate variable
     w_nc_var = w_nc.createVariable(climname, 'f', ('time','lat','lon'))
     w_nc_var.setncatts({'long_name': long_name,\
                    'units': climunits,\
                    'missing_value': missing})
     w_nc.variables[climname][:] = climdata
     
     w_nc.close()




