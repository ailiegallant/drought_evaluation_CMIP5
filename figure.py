#THIS MODULE CREATES FIGURES
#
#EACH SUBROUTINE IS SEPARATED BY A LINE --------
#-----------------------------------------------------------------------------
'''
plot_gpcp_corr

Plots correlations between seasonal time series of GPCP data

'''

def plot_gpcp_corr(sone, stwo):
	from matplotlib import pyplot as plt
	from cartopy import config
	import cartopy.crs as ccrs
	from clim_stats import gpcp_corr
	import numpy as np
	
	gp = gpcp_corr(sone, stwo)
	corr = gp[0]
	lat = gp[1]
	lon = gp[2]
	
	ax = plt.axes(projection=ccrs.PlateCarree())
	
	lev = np.arange(-1.,1.1,0.1)
	pc = plt.contourf(lon, lat, corr, levels=lev, transform=ccrs.PlateCarree())
	cb = plt.colorbar(pc)

	ax.coastlines()
	plt.show()
	

#-----------------------------------------------------------------------------
'''
plot_telecorr

Plots correlations between seasonal time series of GPCP data

'''

def plot_tele_corr():
	from matplotlib import pyplot as plt
	from cartopy import config
	import cartopy.crs as ccrs
	from netcdf_tools import ncextractall
	import numpy as np
	from scipy import stats
	from netCDF4 import num2date

	import numpy as np
	
	with open('sam_1979_2010.dat') as file:
		sam = [[float(digit) for digit in line.split()] for line in file]
	s = 5
	e = 11
	sam = np.array(sam)
	sam = sam[:,1:]
	sam = sam[:,s:e]
	sam = np.mean(sam,axis=1)
	
	with open('nino34_1979_2010.dat') as file:
		enso = [[float(digit) for digit in line.split()] for line in file]
	
	enso = np.array(enso)
	enso = enso[:,1:]
	enso = enso[:,s:e]
	enso = np.mean(enso, axis=1)
	
	with open('iod_1979_2010.dat') as file:
		iod = [[float(digit) for digit in line.split()] for line in file]
	
	iod = np.array(iod)
	iod = iod[:,2]
	iod = iod.reshape(31,12)
	iod = iod[:,s:e]
	iod = np.mean(iod, axis=1)
	print("IOD:",iod.size)
	print("ENSO:",enso.size)
	print("SAM:",sam.size)
	
	obfile = '/Users/ailieg/Data/drought_model_eval_data/data/obs/GPCP/precip.mon.mean.nc'
	obsnc = ncextractall(obfile)
	obsp = obsnc['precip']
	lon = obsnc['lon']
	nlon = lon.size
	lat = obsnc['lat']
	nlat = lat.size
	time = obsnc['time']
	obsmiss = obsnc['precip_missing_value']
	obsp[np.where(obsp == obsmiss)] = np.nan
	
		#convert time units to actual date
	time_u = obsnc['time_units']
	if 'time_calendar' in obsnc.keys(): 
	     cal = obsnc['time_calendar']
	     time = num2date(time,units = time_u, calendar=cal)
	else: time = num2date(time,units = time_u)
	
		#check that the data array begins in January 
	i = 0
	smon = time[i].month
	while (smon != 1):
	    i = i + 1
	    smon = time[i].month 
	    
	#clip the array at the start   
	obsp = obsp[i:,:,:]
	time = time[i:]
	
	#check that the data array ends in December 
	i = len(time) - 1
	emon = time[i].month
	while (emon != 12):
	    i = i - 1
	    emon = time[i].month
	#clip the array at the end    
	obsp = obsp[:i+1,:,:] #remember that Python does not include the final value, so it has to be +1
	time = time[:i+1]
	ntime = time.size
	
	corr = np.zeros((nlat,nlon))
	
	
	for i in range(0,nlon):
		for j in range(0,nlat):
		
			series = obsp[:,j,i]
			lens = series.size
			lens = lens/12
			lens = int(lens)
			series = series.reshape([lens,12])
			srs = series[0:31,s:e]
			series = np.mean(srs, axis=1)
			corriod = stats.pearsonr(iod, series)[0]
			correnso = stats.pearsonr(enso, series)[0]
			corrsam = stats.pearsonr(sam, series)[0]
			
			corr[j,i] = ((corriod**2)+(correnso**2)+(corrsam**2))*100
	
	ax = plt.axes(projection=ccrs.PlateCarree())
	
	lev = np.arange(-100.,105,5)
	pc = plt.contourf(lon, lat, corr, levels=lev, transform=ccrs.PlateCarree())
	cb = plt.colorbar(pc)

	ax.coastlines()
	plt.show()

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
'''
plot_globrainvar

Plots the variance of rainfall compared to mean rainfall

'''

def plot_globrainvar():
	from matplotlib import pyplot as plt
	from cartopy import config
	import cartopy.crs as ccrs
	from netcdf_tools import ncextractall
	import numpy as np
	from scipy import stats
	from netCDF4 import num2date

	import numpy as np
	
	obfile = '/Users/ailieg/Data/drought_model_eval_data/data/obs/GPCP/precip.mon.mean.nc'
	obsnc = ncextractall(obfile)
	obsp = obsnc['precip']
	lon = obsnc['lon']
	nlon = lon.size
	lat = obsnc['lat']
	nlat = lat.size
	time = obsnc['time']
	obsmiss = obsnc['precip_missing_value']
	obsp[np.where(obsp == obsmiss)] = np.nan
	
		#convert time units to actual date
	time_u = obsnc['time_units']
	if 'time_calendar' in obsnc.keys(): 
	     cal = obsnc['time_calendar']
	     time = num2date(time,units = time_u, calendar=cal)
	else: time = num2date(time,units = time_u)
	
		#check that the data array begins in January 
	i = 0
	smon = time[i].month
	while (smon != 1):
	    i = i + 1
	    smon = time[i].month 
	    
	#clip the array at the start   
	obsp = obsp[i:,:,:]
	time = time[i:]
	
	#check that the data array ends in December 
	i = len(time) - 1
	emon = time[i].month
	while (emon != 12):
	    i = i - 1
	    emon = time[i].month
	#clip the array at the end    
	obsp = obsp[:i+1,:,:] #remember that Python does not include the final value, so it has to be +1
	time = time[:i+1]
	ntime = time.size
	
	var = np.zeros((nlat,nlon))
	
	
	for i in range(1,nlon):
		for j in range(1,nlat):
		
			series = obsp[:,j,i]
			lens = series.size
			lens = lens/12
			lens = int(lens)
			series = series.reshape([lens,12])
			series = np.sum(series,axis=1)
			
			sm = np.mean(series)
			sd = np.std(series)
			
			frac = (sd/sm)*100
			frac = np.where(frac > 99.9,99.9,frac)
			
			var[j,i] = frac
			

	ax = plt.axes(projection=ccrs.PlateCarree())
	
	lev = np.arange(-100.,105,5)
	pc = plt.contourf(lon, lat, var, levels=lev, transform=ccrs.PlateCarree())
	cb = plt.colorbar(pc)

	ax.coastlines()
	plt.show()
	
#------------------------------------------------------------------------------	
'''
plot_telecorr

Plots correlations between seasonal time series of GPCP data

'''

def plot_model():
	from matplotlib import pyplot as plt
	import numpy as np
	import pylab
	
	with open('/Users/ailieg/Documents/Papers/In_Prep/Melb_water_catchment/gfdl_2_0_prcp.dat') as file:
		gfdl = [[float(digit) for digit in line.split()] for line in file]
	gfdl = np.array(gfdl)
	gfdlpcp = gfdl[:,2]
	gfdlpcp = np.roll(gfdlpcp,-4)
	gfdlpcp = gfdlpcp.reshape((gfdlpcp.size)/12,12)
	gfdlpcp = np.sum(gfdlpcp,axis=1)
	gfdlpcp = gfdlpcp - np.mean(gfdlpcp)
	
	gfdlqqp = gfdl[:,3]
	gfdlqqp = np.roll(gfdlqqp,-4)
	gfdlqqp = gfdlqqp.reshape((gfdlqqp.size)/12,12)
	gfdlqqp = np.sum(gfdlqqp,axis=1)
	gfdlqqp = gfdlqqp - np.mean(gfdlqqp)

	with open('/Users/ailieg/Documents/Papers/In_Prep/Melb_water_catchment/mirochi_prcp.dat') as file:
		miroc = [[float(digit) for digit in line.split()] for line in file]
	miroc = np.array(miroc)
	mirocpcp = miroc[:,2]
	mirocpcp = np.roll(mirocpcp,-4)
	mirocpcp = mirocpcp.reshape((mirocpcp.size)/12,12)
	mirocpcp = np.sum(mirocpcp,axis=1)
	mirocpcp = mirocpcp - np.mean(mirocpcp)
	
	mirocqqp = miroc[:,3]
	mirocqqp = np.roll(mirocqqp,-4)
	mirocqqp = mirocqqp.reshape((mirocqqp.size)/12,12)
	mirocqqp = np.sum(mirocqqp,axis=1)
	mirocqqp = mirocqqp - np.mean(mirocqqp)
	
	with open('/Users/ailieg/Documents/Papers/In_Prep/Melb_water_catchment/mpi1_prcp.dat') as file:
		mpi = [[float(digit) for digit in line.split()] for line in file]
	mpi = np.array(mpi)
	mpipcp = mpi[:,2]
	mpipcp = np.roll(mpipcp,-4)
	mpipcp = mpipcp.reshape((mpipcp.size)/12,12)
	mpipcp = np.sum(mpipcp,axis=1)
	mpipcp = mpipcp - np.mean(mpipcp)
	
	mpiqqp = mpi[:,3]
	mpiqqp = np.roll(mpiqqp,-4)
	mpiqqp = mpiqqp.reshape((mpiqqp.size)/12,12)
	mpiqqp = np.sum(mpiqqp,axis=1)
	mpiqqp = mpiqqp - np.mean(mpiqqp)

	bins = np.arange(-600,900,100)
	
	hgfdl = np.histogram(gfdlpcp,bins=bins)
	hgfdl = hgfdl[0]/np.sum(hgfdl[0])
	hgfdl = np.concatenate(([0],hgfdl))
	
	hmiroc = np.histogram(mirocpcp,bins=bins)
	hmiroc = hmiroc[0]/np.sum(hmiroc[0])
	hmiroc = np.concatenate(([0],hmiroc))
	
	hmpi = np.histogram(mpipcp,bins=bins)
	hmpi = hmpi[0]/np.sum(hmpi[0])
	hmpi = np.concatenate(([0],hmpi))
	
	pylab.ylim([0,0.7])
	plt.plot(bins,hgfdl,'red')
	plt.plot(bins,hmiroc,'blue')
	plt.plot(bins,hmpi,'green')
	
	plt.show()