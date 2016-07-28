#THIS MODULE CONTAINS TOOLS AND CODE FOR TIME SERIES ANALYSIS
#EACH SUBROUTINE IS SEPARATED BY A LINE --------
#-----------------------------------------------------------------------------



#------------------------------------------------------------------------------

'''
interp_landgrid
Takes an input gridded time series and interpolates to a new grid using a cubic spline.

INPUT
------
oldx - the lon of the grid from which you are interpolating
oldy - the lat of the grid from which you are interpolating
newx - the lon of the grid to which you are interpolating
newy - the lat of the grid to which you are interpolating
data - the data on the grid which you are interpolating

OUTPUT
--------
an interpolated grid where DATA is interpolated on the surface NEWX, NEWY 

HISTORY
-------
20160524 - Created by Ailie Gallant
'''

def interp_landgrid(oldx, oldy, newx, newy, data): 
	from scipy import interpolate
	import numpy as np
	
	gnew_x, gnew_y = np.meshgrid(newx, newy)
	
	gold_x, gold_y = np.meshgrid(oldx, oldy)
	gold_x = gold_x.flatten()
	gold_y = gold_y.flatten()
	oldpoints = np.vstack((gold_x, gold_y)).T
	
	values = data.flatten()
	
	int = interpolate.griddata(oldpoints, values, (gnew_x, gnew_y), method='linear')
	
	#Compute ambiguous land values (e.g. int = 0.25, 0.5 or 0.75). Make land only where
	#int > 0.7
	
	int[np.where(int < 0.7,0,1)]
	
	return int




	
#------------------------------------------------------------------------------

'''
Takes grids from two netcdf files and creates a land mask for interpolation.

INPUT
------
FILE - the file from which the data will be interpolated onto

OUTPUT
--------
an interpolated land mask

HISTORY
-------
20160524 - Created by Ailie Gallant
'''

def make_land_mask(file): 
	from netcdf_tools import ncextractall
	import numpy as np
	
	#Read in file from which new mask will be generated
	maskdat = ncextractall('/Users/ailieg/Data/drought_model_eval_data/lsmask.nc')
	mlon = maskdat['lon']
	mlat = maskdat['lat']
	mask = maskdat['mask']
	mask = mask[0,:,:]
	mask = mask[::-1,:]
	
	#For some reason land in this mask is zero and ocean is one - flip this
	mask[np.where(mask == 1)] = 3
	mask[np.where(mask == 0)] = 1
	mask[np.where(mask == 3)] = 0
	
	
	#Read in file containing the longitude and latitude to which the mask will be interpolated
	dat = ncextractall(file)
	lon = dat['lon']
	lat = dat['lat']
	
	newmask = interp_landgrid(mlon, mlat, lon, lat, mask)
	
	return newmask
	
#------------------------------------------------------------------------------

'''
convert_grid
Takes an input grid and computes it to a new, predefined grid

INPUT
------
oldx - the lon of the grid from which you are interpolating
oldy - the lat of the grid from which you are interpolating
newx - the lon of the grid to which you are interpolating
newy - the lat of the grid to which you are interpolating
data - the data on the grid which you are interpolating in y, x co-ordinates (e.g lat,lon)

OUTPUT
--------
an interpolated grid where DATA is interpolated on the surface NEWY, NEWX 

HISTORY
-------
20160524 - Created by Ailie Gallant
'''

def convert_grid(oldy, oldx, newy, newx, data): 
	from scipy import interpolate
	import numpy as np
	
	intf = interpolate.interp2d(oldy,oldx,data.T, kind='linear')
	
	int = intf(newy,newx)
	
	return int.T


#--------------------------------------------------------------------------

'''
fill_missing(data, lon, lat)

Given an input data grid, longitude and latitude, this function will 
infill any missing data using cubic spline interpolation.

Input:
----------
z - an input array of the form (lat,lon)
y - the latitudes of the data array
x - the longitudes of the data array

Output:
----------
A complete grid
'''

def fill_missing(z, y, x):
	from scipy import interpolate
	import numpy as np
	
	nx = x.size
	ny = y.size
	vals = z.flatten()
	
	g_y, g_x = np.meshgrid(y, x)
	g_x = g_x.flatten()
	g_y = g_y.flatten()
	
	pts = np.vstack((g_y, g_x)).T
	
	if (z[np.isnan(z)]).size > 0:
		pts = pts[~np.isnan(vals),:]
		vals = vals[~np.isnan(vals)]
		
		int = interpolate.griddata(pts, vals, (g_y, g_x), method='linear')
		int = int.reshape(ny,nx)
	else: 
		int = z
	
	return int
#--------------------------------------------------------------------------



'''
grid_to_grid(data, lon, lat, newlon, newlat)

Given an input data grid, longitude and latitude, this function will 
first infill any missing data through linear interpolation, and will then
interpolate onto the new grid.

Input:
----------
data - an input array of the form (lat,lon)
lon - the longitudes of the data array
lat - the latitudes of the data array
newlon - the longitudes of the data array on which the data will be interpolated
newlat - the latitudes of the data array on which the data will be interpolated

Output:
----------
The data array interpolated onto the newlon and newlat points
'''

def grid_to_grid(data, lat, lon, newlat, newlon):
	import numpy as np
	
	#Infill any missing data
	data = fill_missing(data, lat, lon)
	
	#get dimensions
	nlat = lat.size
	nlon = lon.size
	
	#Append values before and after to ensure no missing values
	lon = np.concatenate(([lon[0]-(lon[1]-lon[0])], lon))
	lon = np.concatenate((lon, [lon[nlon]+(lon[1]-lon[0])]))
	
	lat = np.concatenate(([lat[0]-(lat[1]-lat[0])], lat))
	lat = np.concatenate((lat, [lat[nlat]+(lat[1]-lat[0])]))

	if newlat[(newlat.size)-1] < newlat[0]: newlat = newlat[::-1]

	gex = np.zeros((nlat+2,nlon+2))+ np.nan
	gex[1:nlat+1,1:nlon+1] = data

	gex = interp_nearest(gex, lat, lon)

	grid = convert_grid(lat, lon, newlat, newlon, gex)
	
	return grid

#-------------------------------------------------------------------------
'''
interp_nearest(z,y,x)

Given an input data grid, longitude and latitude, this function will 
infill any missing data using nearest neighbour interpolation.

Input:
----------
z - an input array of the form (lat,lon)
y - the latitudes of the data array
x - the longitudes of the data array

Output:
----------
A complete grid
'''

def interp_nearest(z, y, x):
	from scipy import interpolate
	import numpy as np
	
	nx = x.size
	ny = y.size
	vals = z.flatten()
	
	g_y, g_x = np.meshgrid(y, x)
	g_x = g_x.flatten()
	g_y = g_y.flatten()
	
	pts = np.vstack((g_y, g_x)).T
	
	if (z[np.isnan(z)]).size > 0:
		pts = pts[~np.isnan(vals),:]
		vals = vals[~np.isnan(vals)]
		
		int = interpolate.griddata(pts, vals, (g_y, g_x), method='nearest')
		int = int.reshape(ny,nx)
	else: 
		int = z
	
	return int
#--------------------------------------------------------------------------
'''
trim_time_jandec

Function trims a continuous 3D data array so that the time begins in January and 
ends in December

Input:
data - 3D z[t, y, x]
time - [t1,t2,...,tn]

'''	
def trim_time_jandec(data, time):
	

	#check that the data array begins in January 
	i = 0
	smon = time[i].month
	while (smon != 1):
	    i = i + 1
	    smon = time[i].month 
	    
	#clip the array at the start   
	data = data[i:,:,:]
	time = time[i:]
	
	#check that the data array ends in December 
	i = len(time) - 1
	emon = time[i].month
	while (emon != 12):
	    i = i - 1
	    emon = time[i].month
	#clip the array at the end    
	data = data[:i+1,:,:] #remember that Python does not include the final value, so it has to be +1
	time = time[:i+1]
	
	return data, time
	
#--------------------------------------------------------------------------
'''
INTERP_TOCOARSE

Given two sets of input lons, lats and temporal gridded data, this function 
interpolates the finer resolution grid to the coarser resolution. The temporal
portion of the grids need not be equal.

INPUT
a - the first input array of dimensions [time, lat, lon]
alat - the values for the latitudinal dimension of a
alon - the values for the longitudinal dimension of a
b - the second input array of dimensions [time, lat, lon]
blat - the values for the latitudinal dimension of b
blon - the values for the longitudinal dimension of b

Returns a and b on the same grid spacing

'''

def interp_tocoarse(a, alat, alon, b, blat, blon):
	import numpy as np
	
	#Interpolate data to the coarsest grid spacing
	
	#Determine which of the model or observations has the coarsest grid spacing
	alatdim = alat.shape[0]
	alondim = alon.shape[0]
	adimmean = (alatdim + alondim)/2.0
	
	blatdim = blat.shape[0]
	blondim = blon.shape[0]
	bdimmean = (blatdim + blondim)/2.0
	
	#Interpolate to the coarser spacing
	#A larger value for adimmean/bdimmean means that there are more elements in the lat
	#and lon arrays, which means a/b has finer resolution than b/a 
	
	
	#If A is coarser than B (adimmean smaller value than bdimmean)
	if adimmean < bdimmean:
	
		#Define time dimension if applicable and create and loop the array
		if len(b.shape) > 2:
			btime = b.shape[0]
			newb = np.zeros([btime, alatdim, alondim])
		
			#Loop over time dimension and interpolate
			for i in range(0,btime):
				tmpb = b[i,:,:]
				newb[i,:,:] = grid_to_grid(tmpb, blat, blon, alat, alon)
		
			b = newb
			lon = alon
			lat = alat
	
		else:
			b = grid_to_grid(b, blat, blon, alat, alon)
			lon = alon
			lat = alat
			
	#If B is coarser than A (adimmean larger value than bdimmean)
	if adimmean > bdimmean:
	
		#Define time dimension if applicable and create and loop the array
		if len(a.shape) > 2:
			atime = a.shape[0]
			newa = np.zeros([atime, blatdim, blondim])
		
			#Loop over time dimension and interpolate
			for i in range(0,atime):
				tmpa = a[i,:,:]
				newa[i,:,:] = grid_to_grid(tmpa, alat, alon, blat, blon)
		
			a = newa
			lon = blon
			lat = blat
		else: 
			a = grid_to_grid(a, alat, alon, blat, blon)
			lon = blon
			lat = blat
		
	return a, b, lon, lat
	
#--------------------------------------------------------------------------
'''
INTERP_TOGRID

Given two sets of input lons, lats and temporal gridded data, this function 
interpolates the resolution of grid B onto the spacing of grid A. The temporal
portion of the grids need not be equal.

INPUT
a - the first input array of dimensions [time, lat, lon]
alat - the values for the latitudinal dimension of a
alon - the values for the longitudinal dimension of a
b - the second input array of dimensions [time, lat, lon]
blat - the values for the latitudinal dimension of b
blon - the values for the longitudinal dimension of b

Returns a and b on the same grid spacing

'''

def interp_togrid(a, alat, alon, b, blat, blon):
	import numpy as np

	
	#Define time dimension if applicable and create and loop the array
	if len(b.shape) > 2:
		btime = b.shape[0]
		newb = np.zeros([btime, alatdim, alondim])
	
		#Loop over time dimension and interpolate
		for i in range(0,btime):
			tmpb = b[i,:,:]
			newb[i,:,:] = grid_to_grid(tmpb, blat, blon, alat, alon)
	
		b = newb
		lon = alon
		lat = alat
	
	else:
		b = grid_to_grid(b, blat, blon, alat, alon)
		lon = alon
		lat = alat
		
	
	return a, b, lon, lat
	
		
#--------------------------------------------------------------------------	
'''
FOLD_GRID

Appends data onto the end of a lon/lat gridded file to make circular for neat plotting.

INPUT
data - in the form of [time, lat, lon]
lat - lat points for data
lon - lon points for data


'''

def fold_grid(data, lat, lon):
	import numpy as np
	
	d = np.zeros((lat.size, lon.size+1))
	d[:,0:-1] = data
	d[:,-1] = data[:,0]

	lon = np.append(lon,lon[0]+360)
	
	return d, lat, lon
	
	
	