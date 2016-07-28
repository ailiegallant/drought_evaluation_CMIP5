#THIS MODULE ALLOWS YOU TO COMPUTE THE SPI FROM A TIME SERIES USING EITHER SINGLE 
#TIMESERIES OR GRIDDED DATA.
#
#EACH SUBROUTINE IS SEPARATED BY A LINE --------
#-----------------------------------------------------------------------------

#SPI computes the SPI for a given data input

'''
Given an input time series of monthly precipitation data (1D) this code 
computes the Standardised Precipitation Index (SPI) for the time period 
provided in the input to the function. Typical time periods are 3, 6, 12 &
24 months.

INPUT

 SERIES - an input time series of monthly data of length, n.
 X - the period of time over which to calculate the SPI

#Created by Ailie Gallant 15/02/2016
'''

def spi(series, x):
   import numpy as np
   from scipy.stats import gamma
   from timeseries_tools import moving_average
   from scipy.interpolate import interp1d

   # Take monthly data and compute into moving averages of length, X
   # Note that series will be length n-(x-1) following the application of
   # the moving average 
   series = moving_average(series, x)

   #The data must be stratified into months so that the SPI may be computed
   #for each month so data is deseasonalised. 
   #Reshape the array to stratify into months.
   lens = len(series)
   lens = lens/12
   lens = int(lens)
   series = series.reshape([lens,12]) 
      
   #Create dummy arrays to store the SPI data
   z = np.zeros([lens,12],float)
   
   #Compute the SPI, one month at a time
   for i in range(0,12):    
       tmp = series[:,i]
       tmpz = spicalc(tmp)   
       z[:,i] = tmpz

   
   #Reshape the array into it's original time series
   return z


#-----------------------------------------------------------------------------
#SPICALC

'''
Does the transform calculation to compute the SPI. The data
provided for the calculation has been pre-computed into 
appropriate running means etc.

INPUT

 DATA- input data from which the SPI is to be computed

Created by Ailie Gallant 07/04/2016
'''
def spicalc(data):
   import numpy as np
   from scipy.stats import gamma

#remove any NaNs (i.e. missing numbers) from data so only real numbers exist
   tmp = data[~np.isnan(data)]
   
#if there are less than 10 real datum with which to fit the distribution, 
#then return an array of NaN, otherwise, do the calculations
 
   spireturn = np.zeros(len(data))+np.nan

   if len(tmp) > 10:
   
       #compute the shape and scale parameters using more than one non-zero data point
       #otherwise computation of the log will fail
       tmpnonz = tmp[np.where(tmp > 0.0)]
       if len(tmpnonz) > 1:
          A = np.log(np.mean(tmpnonz)) - (np.sum(np.log(tmpnonz))/len(tmpnonz))
          shp = (1.0/(4*A)) * (1 + ((1 + ((4*A)/3) )**0.5))
          scl = np.mean(tmpnonz)/shp 
          gam = gamma.cdf(tmpnonz,shp,scale=scl)
       else:
          #if there are no or one non-zero number, then the probability of non-zero numbers
          #is set as 0 or 1/len(tmp) (depending on len(tmpnonz))
          gam = len(tmpnonz)/len(tmp)
     
       #fit the gamma distribution, G(x), already calculated as gam if there is more than
       #one non-zero number in the time series
       #if there are zero values, the cdf becomes H(x) = q + (1-q)G(x) where q is the
       #probability of a zero value in the time series
       numzero = len(tmp[np.where(tmp == 0.0)])
       if numzero > 0:
          q = numzero/len(tmp)
          gam = q + (1-q)*gam
          gcdf = np.zeros(len(tmp))
          i = np.where(tmp >0.0)
          j = np.where(tmp == 0.0)
          gcdf[i] = gam
          gcdf[j] = q
       else: gcdf = gam

       #define the constants for the approximation    
       c0 = 2.515517
       c1 = 0.802853
       c2 = 0.010328
       d1 = 1.432788
       d2 = 0.189269
       d3 = 0.001308
   
       #compute the SPI values when the gamma cdf is non-uniform
       if len(gcdf[np.where(gcdf == 1.0)]) == 0: 
           t = np.where(gcdf<=0.5,(np.log(1/(gcdf**2)))**0.5,(np.log(1/((1.0-gcdf)**2))**0.5))
           ztmp = (t - ((c0 + c1*t + c2*(t**2))/(1 + d1*t + d2*(t**2) + d3*(t**3)))) 
           s = np.where(gcdf<=0.5, -1*ztmp, ztmp)
           
       #if the grid cell is always dry (i.e. precip of zero, then SPI returns 0s as dry
       #is always "normal"
       else: s = np.zeros(len(gcdf))
       

       spireturn[~np.isnan(data)]=s
       

   
   return spireturn
   
#-----------------------------------------------------------------------------
#SPI_NETCDF_GRID

'''
This script computes the SPI for a set of gridded data that is input in
NetCDF format

Assumes that the dimensions will be time, lat, lon (in that order).

The resultant SPI grid is saved in the same format as the input grid so the
two can be used interchangeably.

INPUT:
  file - the filename with path (as appropriate) of the netcdf file to be read 
         and the SPI to be computed. Must contain precipitation as a variable.
  vname - the variable name to be extracted for calculation of the SPI
 

Created by Ailie Gallant 13/04/2016
'''

def spi_netcdf_grid(rfile, wfile, vname, nspi):
	from netcdf_tools import ncextractall
	from netcdf_tools import ncwrite_climgrid
	from netCDF4 import num2date
	from netCDF4 import date2num
	import numpy as np
	
	#Import and extract the netcdf file, returns a dictionary
	datdict = ncextractall(rfile)
	
	#Read each variable
	data = datdict[vname]
	lon = datdict['lon']
	lat = datdict['lat']
	time = datdict['time']
	
	#convert missing numbers to NaN
	miss = datdict[vname+"_missing_value"] 
	data = np.where(data == miss,np.nan,data)
	
	#convert time units to actual date
	time_u = datdict['time_units']
	if 'time_calendar' in datdict.keys(): 
	     cal = datdict['time_calendar']
	     time = num2date(time,units = time_u, calendar=cal)
	else: time = num2date(time,units = time_u)
	
	trimmed = grid_tools.trim_time_jandec(data, time)
	data = trimmed[0]
	time = trimmed[1]
	
	#Create an empty array to store the SPI data
	spigrid = np.zeros(data.shape)
	
	#Compute the SPI at all locations
	for i in range(0,len(lat)):
	    for j in range(0,len(lon)):
	          tmpdata = data[:,i,j]
	          tmpdata = tmpdata.flatten()
	          tmpspi = spi(tmpdata,nspi)
	          tmpspi = tmpspi.flatten()
	          spigrid[:,i,j] = tmpspi
	
	#convert missing numbers back to a float
	spigrid = np.where(spigrid == np.nan, miss, spigrid)

	#convert time back to original units
	if 'time_calendar' in datdict.keys(): 
	     cal = datdict['time_calendar']
	     time = date2num(time,units = time_u, calendar=cal)
	else: time = date2num(time,units = time_u)
	
	spidescrip = "The Standardised Precipitation Index (SPI) computed as per McKee et al. (1993)"
	spilong_name = str(nspi)+"-month Standardised Precipitation Index"
	spiname = "SPI"+str(nspi)
	
	write = ncwrite_climgrid(wfile, spigrid, spiname, spidescrip, spilong_name, \
	                          miss, "standardised units", time, lon, lat, time_u)
	                          
	return spigrid,lon,lat,time
	
#-----------------------------------------------------------------------------

