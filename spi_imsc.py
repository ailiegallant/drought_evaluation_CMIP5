def plot_spi_metrics(spiname, metname, state, modfile, savename):
	from matplotlib import pyplot as plt
	import numpy as np
	from cartopy import config
	import cartopy.crs as ccrs
	from scipy import interpolate
	
	a = spi_obs_model_compare(spiname, metname, state, modfile)
	dif = a[0]
	sigdif = a[1]
	lat = a[2]
	lon = a[3]
	
	#Make circular
	d = np.zeros((lat.size, lon.size+1))
	d[:,0:-1] = dif
	d[:,-1] = dif[:,0]
	
	s = np.zeros((lat.size, lon.size+1))
	s[:,0:-1] = sigdif
	s[:,-1] = sigdif[:,0]

	lon = np.append(lon,lon[0]+360)

	levels = np.arange(-10,11)
	#d[d >5] = 5
	#d[d<-5] = -5
	
	ax = plt.axes(projection=ccrs.PlateCarree())
	p = plt.contourf(lon, lat, d, levels=levels, transform=ccrs.PlateCarree(), cmap=plt.cm.Spectral)
	
	if np.sum(np.where(s == np.isnan,1,0)) > 0:
		plt.contourf(lon, lat, s, 2, transform=ccrs.PlateCarree(), 
	             colors='none',hatches = '.')
	cbar = plt.colorbar(p)
	cbar.ax.set_ylabel(r'$\mu$')
	ax.coastlines()

	outfile = '/Users/ailieg/Data/drought_model_eval_data/analysis/'+spiname+'.'+metname+'.'+state+'.'+savename
	plt.savefig(outfile, dpi=400, format='png',bbox_inches='tight')
	plt.close()

	return

#----------------------------------------------------------------------------------------
def plot_spi_metrics_coupledmmm(spiname, metname, state, simtype, levs, unit):
	from matplotlib import pyplot as plt
	import numpy as np
	from cartopy import config
	import cartopy.crs as ccrs
	from scipy import interpolate
	
	a = spi_obs_model_compare_mmm(spiname, metname, state, simtype)
	dif = a[0]
	sigdif = a[1]
	lat = a[2]
	lon = a[3]
	
	#Make circular
	d = np.zeros((lat.size, lon.size+1))
	d[:,0:-1] = dif
	d[:,-1] = dif[:,0]
	
	s = np.zeros((lat.size, lon.size+1))
	s[:,0:-1] = sigdif
	s[:,-1] = sigdif[:,0]

	lon = np.append(lon,lon[0]+360)

	levels = levs
	sd = np.std(d)*5
	d[(d > sd)]=np.nan
	d[(d < -1*sd)]=np.nan
	
	flag = np.sum(np.where(np.isnan(s),0,1))
	print("Flag signif:", flag)
	
	plt.ioff()
	
	ax = plt.axes(projection=ccrs.Robinson())
	p = plt.contourf(lon, lat, d, levels=levels, transform=ccrs.Robinson(), cmap=plt.cm.Spectral)
	
	if flag != 0:
		plt.contourf(lon, lat, s, 2, transform=ccrs.Robinson(), 
	             colors='none',hatches = '.')
	cbar = plt.colorbar(p)
	cbar.ax.set_ylabel(unit)
	ax.coastlines()
	#	plt.show()
	outfile = '/Users/ailieg/Data/drought_model_eval_data/analysis/'+spiname+'.'+metname+'.'+state+'.'+simtype+'_MMM.png'
	plt.savefig(outfile, dpi=400, format='png',bbox_inches='tight')
	plt.close()

	return d,s


#-----------------------------------------------------------------------------------
#SPI_OBS_MODEL_COMPARE
'''
spi_obs_model_compare

Purpose:
---------
Compares the metrics from the CMIP5 or AMIP simulations to the synthetic 95th percentile
confidence interval.

*Comparison metrics:
-Perkins skill score (Sum(1,n) min(Z1, Z2)) where Zm, Zo are the relative frequencies from
datasets 1 and 2.
- Gamma distribution shape and scale parameters
- Mean and median of the data

Inputs:
---------

Outputs:
---------
Saved NetCDF files of the metrics for each of SPI3, 6, 12 and 24

History:
---------
20160524 - Created by Ailie Gallant, Monash University

'''
def spi_obs_model_compare_mmm(spiname, metname, state, simtype):
	from netcdf_tools import ncextractall
	from netCDF4 import Dataset
	from netCDF4 import date2num
	from netCDF4 import num2date
	import time as comptime
	import numpy as np
	from grid_tools import modgrid_to_obsgrid
	from grid_tools import fill_missing
	from matplotlib import pyplot as plt
	from grid_tools import make_land_mask
	
	if state == 'mod': sind = 0
	if state == 'sev': sind = 1
	
	#Metrics for the surrogate data set
	print("Reading metrics from bootstrapped data...")
	surrmetname = '/Users/ailieg/Data/drought_model_eval_data/analysis/'+spiname+'.mon.GPCPsynthetic_landonly.nc'
	surrnc = ncextractall(surrmetname)
	surrlat = surrnc['lat']
	surrlon = surrnc['lon']
	surrmet = surrnc[metname]
	surrmetone = surrmet[sind,:,:,:]

	print("Reading metrics from bootstrapped data...")
	surrmetname = '/Users/ailieg/Data/drought_model_eval_data/analysis/'+spiname+'.mon.GPCPsynthetic_landonly_2.nc'
	surrnc = ncextractall(surrmetname)
	surrmet = surrnc[metname]
	surrmettwo = surrmet[sind,:,:,:]
	
	surrmet = np.zeros((300,surrlat.size,surrlon.size))
	surrmet[0:100,:,:] = surrmetone
	surrmet[100:,:,:] = surrmettwo
	
	if spiname == 'SPI3': spiind = 0
	if spiname == 'SPI6': spiind = 1
	if spiname == 'SPI12': spiind = 2
	if spiname == 'SPI24': spiind = 3

	#Metrics for the observed metrics
	print("Reading metrics from observed data...")
	obmetname = '/Users/ailieg/Data/drought_model_eval_data/analysis/spi_obs_metrics_nomask.nc'
	obnc = ncextractall(obmetname)
	oblat = obnc['lat']
	oblat = oblat[::-1]
	oblon = obnc['lon']
	obmet = obnc[metname]
	obmet = (obmet[spiind,sind,:,:]).T
	
	#Read in model data and take the multi-model mean
	if simtype == 'cmip':	
		modfiles = [spiname+'_ACCESS1-0_historical_r1i1p1_stats_.nc',\
				spiname+'_GISS-E2-R_historical_r6i1p1_stats_.nc',\
				spiname+'_CCSM4_historical_r1i1p1_stats_.nc',\
				spiname+'_HadGEM2-CC_historical_r1i1p1_stats_.nc',\
				spiname+'_CanESM2_historical_r1i1p1_stats_.nc',\
				spiname+'_FGOALS-s2_historical_r1i1p1_stats_.nc',\
				spiname+'_MPI-ESM-P_historical_r1i1p1_stats_.nc',\
				spiname+'_GFDL-CM3_historical_r1i1p1_stats_.nc']
	#'SPI3_IPSL-CM5B-LR_historical_r1i1p1_stats_.nc',\
	
	if simtype == 'cmipamip':	
		modfiles = [spiname+'_HadGEM2-CC_historical_r1i1p1_stats_.nc',\
				spiname+'_FGOALS-s2_historical_r1i1p1_stats_.nc',\
				spiname+'_GFDL-CM3_historical_r1i1p1_stats_.nc']
	
	if simtype == 'amip':	
		modfiles = [spiname+'_FGOALS-s2_amip_r1i1p1_stats_.nc',\
				spiname+'_HadGEM2-A_amip_r1i1p1_stats_.nc',\
				spiname+'_GFDL-CM3_amip_r1i1p1_stats_.nc',\
				spiname+'_NorESM1-M_amip_r1i1p1_stats_.nc']
	
	if simtype == 'amipcmip':	
		modfiles = [spiname+'_FGOALS-s2_amip_r1i1p1_stats_.nc',\
				spiname+'_HadGEM2-A_amip_r1i1p1_stats_.nc',\
				spiname+'_GFDL-CM3_amip_r1i1p1_stats_.nc']
				
	
	modfiles = np.array(modfiles)
	nmod = modfiles.size

	
	modmmm = np.zeros((nmod, oblat.size, oblon.size))
	
	for i in range(0,nmod):
		file = '/Users/ailieg/Data/drought_model_eval_data/analysis/'+modfiles[i]
		ncmod = ncextractall(file)
		modmet = ncmod[metname]
		modlon = ncmod['lon']
		modlat = ncmod['lat']
		modmet = modmet[sind,:,:]
		modmet = modgrid_to_obsgrid(modmet, modlat, modlon, oblat, oblon)
		modmmm[i,:,:] = modmet
	
	m = np.mean(modmmm, axis=0)
	
	miss = -9e35
	sigdif = m - obmet
	sigdif[obmet <= miss] = miss
	sigdif[modmet <= miss] = miss
	
	surrmin = np.min(surrmet, axis=0)
	surrmax = np.max(surrmet, axis=0)
	
	sigdif[modmet < surrmin] = miss
	sigdif[modmet > surrmax] = miss
	
		#Make a land mask
	obfile = '/Users/ailieg/Data/drought_model_eval_data/data/obs/GPCP/precip.mon.mean.nc'
	mask = make_land_mask(obfile)
	sigdif[(sigdif != miss)] = 1.0
	sigdif[(mask < 1)] = miss
	sigdif[(sigdif == miss)] = np.nan

	obmet[np.where(obmet < miss)] = np.nan
	dif = m - obmet

	return dif, sigdif, oblat, oblon
	
	
	
#------------------------------------------------------------------------------------

def plot_multi_metrics(spiname, metname, state):
	import numpy as np

	modfile = ['SPI3_ACCESS1-0_historical_r1i1p1_stats_.nc',\
			'SPI3_GISS-E2-R_historical_r6i1p1_stats_.nc',\
			'SPI3_CCSM4_historical_r1i1p1_stats_.nc',\
			'SPI3_HadGEM2-CC_historical_r1i1p1_stats_.nc',\
			'SPI3_CanESM2_historical_r1i1p1_stats_.nc',\
			'SPI3_IPSL-CM5B-LR_historical_r1i1p1_stats_.nc',\
			'SPI3_FGOALS-s2_historical_r1i1p1_stats_.nc',\
			'SPI3_MPI-ESM-P_historical_r1i1p1_stats_.nc',\
			'SPI3_GFDL-CM3_historical_r1i1p1_stats_.nc']
			
	savename = ['ACCESS1-0_historical_r1i1p1.png',\
			'GISS-E2-R_historical_r6i1p1.png',\
			'CCSM4_historical_r1i1p1.png',\
			'HadGEM2-CC_historical_r1i1p1.png',\
			'CanESM2_historical_r1i1p1.png',\
			'IPSL-CM5B-LR_historical_r1i1p1.png',\
			'FGOALS-s2_historical_r1i1p1.png',\
			'MPI-ESM-P_historical_r1i1p1.png',\
			'GFDL-CM3_historical_r1i1p1.png']
					
	modfile = np.array(modfile)
	nmod = modfile.size

	for i in range(0,nmod): 
		plot = plot_spi_metrics(spiname, metname, state, modfile[i], savename[i])
		
#------------------------------------------------------------------------------------		

#----------------------------------

def plot_mmm_metrics(sim, state):
	import numpy as np


	lev = np.arange(-5,6)
	u = r'$\xi$' 
	plot = plot_spi_metrics_coupledmmm('SPI3', 'shape', state, sim, lev, u)
	#plot = plot_spi_metrics_coupledmmm('SPI12', 'shape', state, sim, lev, u)
	
	lev = np.arange(-5,6)
	u = r'$\sigma$' 
	plot = plot_spi_metrics_coupledmmm('SPI3', 'scale', state, sim, lev, u)
	#plot = plot_spi_metrics_coupledmmm('SPI12', 'scale', state, sim, lev, u)
	
	lev = np.arange(-5,6)
	u = r'$\mu$' 
	plot = plot_spi_metrics_coupledmmm('SPI3', 'loc', state, sim, lev, u)
	#plot = plot_spi_metrics_coupledmmm('SPI12', 'loc', state, sim, lev, u)
	
	
		
		
		
		
		
		
		
	return	
	
#------------------------------------------------------------------------------------		

#----------------------------------

def plot_gev_dist():
	from netcdf_tools import ncextractall
	from matplotlib import pyplot as plt
	from scipy import stats
	import numpy as np
	from grid_tools import make_land_mask
	
	#plt.ion()
	obfile = '/Users/ailieg/Data/drought_model_eval_data/data/obs/GPCP/precip.mon.mean.nc'
	mask = make_land_mask(obfile)
	
	obmetname = '/Users/ailieg/Data/drought_model_eval_data/analysis/spi_obs_metrics_nomask.nc'
	obnc = ncextractall(obmetname)
	obshape = obnc['shape']
	obshape = (obshape[0,0,:,:]).T
	obshape = obshape[(mask == 1)]
	obshape = obshape.flatten()
	shapemean = np.mean(obshape) 
	shapelo = np.percentile(obshape,10)
	shapehi = np.percentile(obshape,90)
	
	obscale = obnc['scale']
	obscale = (obscale[0,0,:,:]).T
	obscale = obscale[(mask == 1)]
	obscale = obscale.flatten()
	scalemean = np.mean(obscale) 
	scalelo = np.percentile(obscale,10)
	scalehi = np.percentile(obscale,90)
	
	obloc = obnc['loc']
	obloc = (obloc[0,0,:,:]).T
	obloc = obloc[(mask == 1)]
	obloc = obloc.flatten()
	locmean = np.mean(obloc) 
	loclo = np.percentile(obloc,10)
	lochi = np.percentile(obloc,90)
	
	#shape plot
	sz = 1000000
	normbins = np.arange(-0.5,5.6,0.1)#normbins[0:-1]
	
	randnorm = stats.genextreme.rvs(shapemean, loc = locmean, scale = scalemean, size=sz)
	h = np.histogram(randnorm, bins=normbins)
	hnorm = h[0]/sz
	normbins = h[1]
	
	randlo = stats.genextreme.rvs(shapemean, loc = locmean, scale = scalelo, size=sz)
	h = np.histogram(randlo, bins=normbins)
	hlo = h[0]/sz
	randhi = stats.genextreme.rvs(shapemean, loc = locmean, scale = scalehi, size=sz)
	h = np.histogram(randhi, bins=normbins)
	hhi = h[0]/sz
	
	normbins = normbins[0:-1]
	#print(normbins.shape, hnorm.shape)
	p=plt.plot(normbins,hnorm, color='black')
	plt.setp(p, linewidth=2)
	#q=plt.plot(normbins, hlo, color='red')
	#plt.setp(q, linewidth=2)
	#r=plt.plot(normbins, hhi, color='blue')
	#plt.setp(r, linewidth=2)
	#plt.title('Scale ('+r'$\sigma$'+')')

	plt.show()
	
	return
#----------------------------------------------------------------------------------

def plot_location_gev(ilat, jlon, b):
	from netcdf_tools import ncextractall
	from matplotlib import pyplot as plt
	from scipy import stats
	import numpy as np
	from grid_tools import make_land_mask
	from timeseries_tools import average_run
	
	obfile = '/Users/ailieg/Data/drought_model_eval_data/data/obs/GPCP/precip.mon.mean.nc'
	mask = make_land_mask(obfile)
	print("Mask value:", mask[ilat,jlon])
	if mask[ilat, jlon] == 1: print("Mask says location is on land")
	if mask[ilat, jlon] == 0:print("Mask says location is in the ocean")
	
	obfile = '/Users/ailieg/Data/drought_model_eval_data/data/obs/GPCP/precip.mon.mean.nc'
	obsnc = ncextractall(obfile)
	oblon = obsnc['lon']
	oblat = obsnc['lat']
	oblat = oblat[::-1]
	
	
	obs = compute_obs_pdist(ilat,jlon)
	h = np.histogram(obs[0], bins=b)
	hob = h[0]/obs[0].size
	obsbins = h[1]
	obidx = obs[1] 
	mrun=average_run(obidx)
	print("Mean run-len from OBS:", mrun)
	cmippath = '/Users/ailieg/Data/drought_model_eval_data/data/CMIP5/'
	
	prfile = ['ACCESS1-0/r1i1p1/pr/pr_Amon_ACCESS1-0_historical_r1i1p1_185001-200512.nc',\
	'CanESM2/r1i1p1/pr/pr_Amon_CanESM2_historical_r1i1p1_185001-200512.nc',\
	'GFDL-CM3/r1i1p1/pr/pr_Amon_GFDL-CM3_historical_r1i1p1_186001-200512.nc',\
	'HadGEM2-CC/r1i1p1/pr/pr_Amon_HadGEM2-CC_historical_r1i1p1_185912-200511.nc',\
	'MPI-ESM-P/r1i1p1/pr/pr_Amon_MPI-ESM-P_historical_r1i1p1_185001-200512.nc',\
	'CCSM4/r1i1p1/pr/pr_Amon_CCSM4_historical_r1i1p1_185001-200512.nc',\
	'FGOALS-s2/r1i1p1/pr/pr_Amon_FGOALS-s2_historical_r1i1p1_185001-200512.nc',\
	'GISS-E2-R/r6i1p1/pr/pr_Amon_GISS-E2-R_historical_r6i1p1_185001-200512.nc',\
	'NorESM1-M/r1i1p1/pr/pr_Amon_NorESM1-M_historical_r1i1p1_185001-200512.nc']
	#'IPSL-CM5B-LR/r1i1p1/pr/pr_Amon_IPSL-CM5B-LR_historical_r1i1p1_185001-200512.nc',\
	
	nmod = 9
	hmod = np.zeros((nmod,len(b)-1))
	io = oblat[ilat]
	jo = oblon[jlon]
	runlen = np.zeros(nmod)

	for i in range(0,9):
		modfile = cmippath+prfile[i]
		ncmod = ncextractall(modfile)
		modlon = ncmod['lon']
		modlat = ncmod['lat']
		idx =(np.abs(modlat-io)).argmin()
		idy = np.abs(modlon-jo).argmin()

		tmp = compute_mod_pdist(idx,idy, modfile)
		tmpmod = np.zeros((5,tmp[0].size))
		tmprun = np.zeros(5)
		tmpmod[0,:] = tmp[0]
		tmprun[0] = average_run(tmp[1])
		
		tmp = compute_mod_pdist(idx-1,idy-1, modfile)
		tmprun[1] = average_run(tmp[1])
		#tmpmod[1,:] = tmp[0] 
		
		tmp = compute_mod_pdist(idx+1,idy+1, modfile)
		tmprun[2] = average_run(tmp[1])
		#tmpmod[2,:] = tmp[0]
		
		tmp = compute_mod_pdist(idx-1,idy+1, modfile)
		tmprun[3] = average_run(tmp[1])
		#tmpmod[3,:] = tmp[0]
		
		tm = compute_mod_pdist(idx+1,idy-1, modfile)
		tmprun[4] = average_run(tmp[1])
		#tmpmod[4,:] = tmp[0]
		
		
		runlen[i]=np.mean(tmprun)
		
		
		#tmpmod = np.mean(tmpmod,axis=0)
		
		#h = np.histogram(tmpmod, bins=b)
		#hmod[i,:] = h[0]/tmpmod.size
	print("Mean run-len from MMM:", np.mean(runlen))
	#hmod = np.mean(hmod,axis=0)
	
	obsbins = obsbins[0:-1]
	p=plt.plot(obsbins,hob, color='black')
	plt.setp(p, linewidth=2)
	q=plt.plot(obsbins, hmod, color='red')
	plt.setp(q, linewidth=2)
	plt.title('Western Canada')

	#plt.show()
	
	return
	
#--------------------------------------------------------------------------------------
def compute_obs_pdist(ilat,jlon):
	import numpy as np
	from netcdf_tools import ncextractall
	from netCDF4 import Dataset
	from netCDF4 import date2num
	from netCDF4 import num2date
	import grid_tools
	from timeseries_tools import moving_average
	
	#Read in GPCP observations
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

	#Extract the observed GPCP data for lat, lon
	savedo = obsp[:,ilat,jlon]

	#Smooth GPCP synthetic data as for the SPI
	tmpobs = moving_average(savedo, 3)
	
	#Extract the observed SPI data for GPCP for lat, lon
	tmpspiob = spi(savedo,3)
	tmpspiob = tmpspiob
	tmpspiob = tmpspiob.flatten()
	
	#Use real data only
	tmpobs = tmpobs[~np.isnan(tmpobs)]
	tmpspiob = tmpspiob[~np.isnan(tmpspiob)]
	
	#Compute shape, scale, location, mean and median metrics
	dmodobs = tmpobs[np.where(tmpspiob < -1.0)]
	idx = np.arange(0,tmpobs.size)
	dobidx = idx[np.where(tmpspiob < -1.0)]
	
	return dmodobs, dobidx
	
#--------------------------------------------------------------------------------------
def compute_mod_pdist(latval,lonval, modfile):
	from scipy import stats
	import numpy as np
	from netcdf_tools import ncextractall
	from netCDF4 import Dataset
	from netCDF4 import date2num
	from netCDF4 import num2date
	from timeseries_tools import moving_average
	from grid_tools import modgrid_to_obsgrid

	
	#Read model precip data file
	datdict = ncextractall(modfile)
	
	#Read each variable
	data = datdict['pr']
	lon = datdict['lon']
	nlon = lon.size
	lat = datdict['lat']
	nlat = lat.size
	time = datdict['time']
	ilat = (np.abs(lat-latval)).argmin()
	jlon = (np.abs(lon-lonval)).argmin()
	
	#convert time units to actual date
	time_u = datdict['time_units']
	if 'time_calendar' in datdict.keys(): 
	     cal = datdict['time_calendar']
	     time = num2date(time,units = time_u, calendar=cal)
	else: time = num2date(time,units = time_u)
	
	#convert missing numbers to NaN
	miss = 1e10
	data = np.where(data > miss,np.nan,data)
	
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
	ntime = time.size
	
	#Convert model data to same units as obs
	data = data * 86400.

	#Extract the observed GPCP data for lat, lon
	tmpdata = data[:,ilat,jlon]
	tmpdata = tmpdata.flatten()
	
	tmpspi = spi(tmpdata,3)
	tmpspi = tmpspi.flatten()
	
	#Smooth GPCP synthetic data as for the SPI
	tmpdata = moving_average(tmpdata, 3)
	
	#Use real data only
	tmpdata = tmpdata[~np.isnan(tmpdata)]
	tmpspi = tmpspi[~np.isnan(tmpspi)]

	
	#Compute shape, scale, location, mean and median metrics
	dmod = tmpdata[np.where(tmpspi < -1.0)]
	idx = np.arange(0,tmpdata.size)
	dmodidx = idx[np.where(tmpspi < -1.0)]
	
	return dmod, dmodidx
	
#--------------------------------------------------------------------------------------
def plot_mean_rain_bias():
	from netcdf_tools import ncextractall
	from matplotlib import pyplot as plt
	from scipy import stats
	from scipy import ndimage
	import numpy as np
	from grid_tools import modgrid_to_obsgrid
	from cartopy import config
	import cartopy.crs as ccrs
	
	obfile = '/Users/ailieg/Data/drought_model_eval_data/data/obs/GPCP/precip.mon.mean.nc'
	obsnc = ncextractall(obfile)
	obs = obsnc['precip']
	lon = obsnc['lon']
	lat = obsnc['lat']
	lat = lat[::-1]
	obs = ndimage.filters.uniform_filter(obs,size=[3,1,1])
	obs = np.mean(obs,axis=0)

	
	cmippath = '/Users/ailieg/Data/drought_model_eval_data/data/CMIP5/'
	
	prfile = ['ACCESS1-0/r1i1p1/pr/pr_Amon_ACCESS1-0_historical_r1i1p1_185001-200512.nc',\
	'CanESM2/r1i1p1/pr/pr_Amon_CanESM2_historical_r1i1p1_185001-200512.nc',\
	'GFDL-CM3/r1i1p1/pr/pr_Amon_GFDL-CM3_historical_r1i1p1_186001-200512.nc',\
	'HadGEM2-CC/r1i1p1/pr/pr_Amon_HadGEM2-CC_historical_r1i1p1_185912-200511.nc',\
	'MPI-ESM-P/r1i1p1/pr/pr_Amon_MPI-ESM-P_historical_r1i1p1_185001-200512.nc',\
	'CCSM4/r1i1p1/pr/pr_Amon_CCSM4_historical_r1i1p1_185001-200512.nc',\
	'FGOALS-s2/r1i1p1/pr/pr_Amon_FGOALS-s2_historical_r1i1p1_185001-200512.nc',\
	'GISS-E2-R/r6i1p1/pr/pr_Amon_GISS-E2-R_historical_r6i1p1_185001-200512.nc',\
	'NorESM1-M/r1i1p1/pr/pr_Amon_NorESM1-M_historical_r1i1p1_185001-200512.nc']
	#'IPSL-CM5B-LR/r1i1p1/pr/pr_Amon_IPSL-CM5B-LR_historical_r1i1p1_185001-200512.nc',\
	
	nmod = 9

	mod = np.zeros((nmod,72,144))
	
	for k in range(0,9):
		modfile = cmippath+prfile[k]
		ncmod = ncextractall(modfile)
		modlon = ncmod['lon']
		modlat = ncmod['lat']
		moddata = ncmod['pr']
		moddata = moddata*86400.
		moddata = ndimage.filters.uniform_filter(moddata,size=[3,1,1])
		moddata = np.mean(moddata, axis=0)
		mod[k,:,:] = modgrid_to_obsgrid(moddata, modlat, modlon, lat, lon)
	mod = np.mean(mod, axis=0)
	bias = mod - obs
	
	ax = plt.axes(projection=ccrs.Robinson())
	levels = np.arange(-10,11)
	p = plt.contourf(lon, lat, bias, levels=levels, transform=ccrs.Robinson(), cmap=plt.cm.Spectral)

	cbar = plt.colorbar(p)
	unit = 'mm/day'
	cbar.ax.set_ylabel(unit)
	ax.coastlines()

	plt.show()
	
	return
	
#--------------------------------------------------------------------------------------
def plot_bias_zonal_land():
	from netcdf_tools import ncextractall
	from matplotlib import pyplot as plt
	from scipy import stats
	from scipy import ndimage
	import numpy as np
	from grid_tools import modgrid_to_obsgrid
	from grid_tools import make_land_mask

	
	obfile = '/Users/ailieg/Data/drought_model_eval_data/data/obs/GPCP/precip.mon.mean.nc'
	obsnc = ncextractall(obfile)
	obs = obsnc['precip']
	lon = obsnc['lon']
	lat = obsnc['lat']
	lat = lat[::-1]
	obs = ndimage.filters.uniform_filter(obs,size=[3,1,1])
	obs = np.mean(obs,axis=0)

	
	cmippath = '/Users/ailieg/Data/drought_model_eval_data/data/CMIP5/'
	
	prfile = ['ACCESS1-0/r1i1p1/pr/pr_Amon_ACCESS1-0_historical_r1i1p1_185001-200512.nc',\
	'CanESM2/r1i1p1/pr/pr_Amon_CanESM2_historical_r1i1p1_185001-200512.nc',\
	'GFDL-CM3/r1i1p1/pr/pr_Amon_GFDL-CM3_historical_r1i1p1_186001-200512.nc',\
	'HadGEM2-CC/r1i1p1/pr/pr_Amon_HadGEM2-CC_historical_r1i1p1_185912-200511.nc',\
	'MPI-ESM-P/r1i1p1/pr/pr_Amon_MPI-ESM-P_historical_r1i1p1_185001-200512.nc',\
	'CCSM4/r1i1p1/pr/pr_Amon_CCSM4_historical_r1i1p1_185001-200512.nc',\
	'FGOALS-s2/r1i1p1/pr/pr_Amon_FGOALS-s2_historical_r1i1p1_185001-200512.nc',\
	'GISS-E2-R/r6i1p1/pr/pr_Amon_GISS-E2-R_historical_r6i1p1_185001-200512.nc',\
	'NorESM1-M/r1i1p1/pr/pr_Amon_NorESM1-M_historical_r1i1p1_185001-200512.nc']
	#'IPSL-CM5B-LR/r1i1p1/pr/pr_Amon_IPSL-CM5B-LR_historical_r1i1p1_185001-200512.nc',\
	
	nmod = 9

	mod = np.zeros((nmod,72,144))
	
	for k in range(0,9):
		modfile = cmippath+prfile[k]
		ncmod = ncextractall(modfile)
		modlon = ncmod['lon']
		modlat = ncmod['lat']
		moddata = ncmod['pr']
		moddata = moddata*86400.
		moddata = ndimage.filters.uniform_filter(moddata,size=[3,1,1])
		moddata = np.mean(moddata, axis=0)
		mod[k,:,:] = modgrid_to_obsgrid(moddata, modlat, modlon, lat, lon)
	mod = np.mean(mod, axis=0)
	bias = mod - obs
	
	
	metname='loc'
	spiind=0
	sind=0
	spiname='SPI3'
		#Metrics for the observed metrics
	print("Reading metrics from observed data...")
	obmetname = '/Users/ailieg/Data/drought_model_eval_data/analysis/spi_obs_metrics_nomask.nc'
	obnc = ncextractall(obmetname)
	oblat = obnc['lat']
	oblat = oblat[::-1]
	oblon = obnc['lon']
	obmet = obnc[metname]
	obmet = (obmet[spiind,sind,:,:]).T
	

	modfiles = [spiname+'_ACCESS1-0_historical_r1i1p1_stats_.nc',\
				spiname+'_GISS-E2-R_historical_r6i1p1_stats_.nc',\
				spiname+'_CCSM4_historical_r1i1p1_stats_.nc',\
				spiname+'_HadGEM2-CC_historical_r1i1p1_stats_.nc',\
				spiname+'_CanESM2_historical_r1i1p1_stats_.nc',\
				spiname+'_FGOALS-s2_historical_r1i1p1_stats_.nc',\
				spiname+'_MPI-ESM-P_historical_r1i1p1_stats_.nc',\
				spiname+'_GFDL-CM3_historical_r1i1p1_stats_.nc']

				
	
	modfiles = np.array(modfiles)
	nmod = modfiles.size
	
	modmmm = np.zeros((nmod, oblat.size, oblon.size))
	
	for i in range(0,nmod):
		file = '/Users/ailieg/Data/drought_model_eval_data/analysis/'+modfiles[i]
		ncmod = ncextractall(file)
		modmet = ncmod[metname]
		modlon = ncmod['lon']
		modlat = ncmod['lat']
		modmet = modmet[sind,:,:]
		modmet = modgrid_to_obsgrid(modmet, modlat, modlon, oblat, oblon)
		modmmm[i,:,:] = modmet
	
	m = np.mean(modmmm, axis=0)
	
	drought_bias = m-obmet
	
	print(bias.shape)
	print(drought_bias.shape)
	
	obfile = '/Users/ailieg/Data/drought_model_eval_data/data/obs/GPCP/precip.mon.mean.nc'
	mask = make_land_mask(obfile)
	drought_bias[drought_bias > 50] = np.nan
	drought_bias[mask < 0.7] = np.nan
	bias[mask < 0.7] = np.nan
	
	
	zmb = np.zeros(72)
	zmdb = np.zeros(72)
	
	for i in range(0,72):
		tb = bias[i,:]
		zmb[i] = np.mean(tb[~np.isnan(tb)])
		tb = drought_bias[i,:]
		zmdb[i] = np.mean(tb[~np.isnan(tb)])
	
	p=plt.plot(oblat, zmb, 'black')
	q=plt.plot(oblat, zmdb, 'red')
	plt.setp(p, linewidth=2)
	plt.setp(q, linewidth=2)
	plt.show()
	return
'''
synth_spi_confint

Purpose:
---------
Computes the two-tailed 95% confidence intervals of the comparison metrics* from the synthetic
precipitation time series. 

*Comparison metrics:
-Perkins skill score (Sum(1,n) min(Z1, Z2)) where Zm, Zo are the relative frequencies from
datasets 1 and 2.
- Gamma distribution shape and scale parameters
- Mean and median of the data

Inputs:
---------

Outputs:
---------
Saved NetCDF files of the metrics for each of SPI3, 6, 12 and 24

History:
---------
20160524 - Created by Ailie Gallant, Monash University

'''
def obs_spi_metrics():
	from scipy import stats
	import numpy as np
	from netcdf_tools import ncextractall
	from netCDF4 import Dataset
	from netCDF4 import date2num
	from netCDF4 import num2date
	import grid_tools
	import sys
	from timeseries_tools import moving_average
	from matplotlib import pyplot as plt
	import time as comptime
	
	print("Starting calcs...")
	
	#Read in GPCP observations
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
	print("Observed data has been extracted")
	
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
	
	print("Observed data clipped to complete year as necessary. New time dimension is",time.size)
	
	spival = np.array([3,6,12,24])
	nspistr = np.array(['SPI3','SPI6','SPI12','SPI24'])
	nspi = spival.size
	
	mask = grid_tools.make_land_mask(obfile)
	
	#Create empty arrays to save metrics for the obs and synthetic data
	#Arrays are of the form (state,lon,lat), where state indicates moderate
	#or severe drought.
	oshape = np.zeros((nspi,2,nlon,nlat))+ obsmiss
	oloc = np.zeros((nspi,2,nlon,nlat))+ obsmiss
	oscale = np.zeros((nspi,2,nlon,nlat))+ obsmiss
	omean = np.zeros((nspi,2,nlon,nlat))+ obsmiss
	omedian = np.zeros((nspi,2,nlon,nlat))+ obsmiss
	
	#Loop through LON and LAT (compute land only)
	for j in range(0,nlon):
	
		for k in range(0,nlat):

				#Extract the observed GPCP data for lat, lon
				savedo = obsp[:,k,j]

				for i in range(0,nspi):
					#Smooth GPCP synthetic data as for the SPI
					tmpobs = moving_average(savedo, spival[i])
			
					#Extract the observed SPI data for GPCP for lat, lon
					tmpspiob = spi(savedo,spival[i])
					tmpspiob = tmpspiob
					tmpspiob = tmpspiob.flatten()
				
					#Use real data only
					tmpobs = tmpobs[~np.isnan(tmpobs)]
					tmpspiob = tmpspiob[~np.isnan(tmpspiob)]
				
					if tmpobs.size != tmpspiob.size:
						print(tmpobs.size, tmpspiob.size)
						print("Error at lon:", j, "lat:", k)
						sys.exit("Error: Observed PRCP & SPI not equal length!")
					
					#Compute shape, scale, location, mean and median metrics
					dmodobs = tmpobs[np.where(tmpspiob < -1.0)]
					nmod = dmodobs.size
					dsevobs = tmpobs[np.where(tmpspiob < -1.5)]
					nsev = dsevobs.size
					if nmod < 10: print("Obs mod droughts:",nmod) 
					if nsev < 10: print("Obs sev droughts:",nsev)
					
					if nmod > 10:
						try:
							#gevmodpar = stats.genextreme.fit(dmodobs)
							gevmodpar = stats.genpareto.fit(dmodobs)
						except ValueError:
							gevmodpar = np.zeros(3) + obsmiss
							print("ValueError encountered, missing values applied ")
						
						oshape[i,0,j,k] = gevmodpar[0]
						oloc[i,0,j,k] = gevmodpar[1]
						oscale[i,0,j,k] = gevmodpar[2]
						omean[i,0,j,k] = np.mean(dmodobs)
						omedian[i,0,j,k] = np.median(dmodobs)
					else:
						continue
						
					if nsev > 10:
						try:
							#gevsevpar = stats.genextreme.fit(dsevobs)
							gevsevpar = stats.genpareto.fit(dsevobs)
						except ValueError:
							gevsevpar = np.zeros(3) + obsmiss
							print("ValueError encountered, missing values applied")
						
						oshape[i,1,j,k] = gevsevpar[0]
						oloc[i,1,j,k] = gevsevpar[1]
						oscale[i,1,j,k] = gevsevpar[2]
						omean[i,1,j,k] = np.mean(dsevobs)
						omedian[i,1,j,k] = np.median(dsevobs)
					else:
						continue

		print("Longitude ",j+1, "of", nlon, "complete")
		#your code here    

#Create a NetCDF file

	#Create NetCDF file to write
	outfile = '/Users/ailieg/Data/drought_model_eval_data/analysis/spi_obs_metrics_pareto.nc'
	print(outfile)
	w_nc = Dataset(outfile, 'w', format='NETCDF4')

	#File description
	long_name = "The statistical metrics of GEV shape, location and scale parameters, mean and median",\
	             " computed from the observed monthly rainfall data."
	w_nc.description = long_name
	
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
	
	#File dimensions for STATE 
	w_nc.createDimension('state', 2)
	state = np.array(['Moderate', 'Severe'],dtype='str')
	w_nc_nsyn = w_nc.createVariable('state', state.dtype,('state',))
	
	w_nc_nsyn.setncatts({'long_name': 'Moderate drought (SPI < -1) or Severe drought (SPI < -1.5)'})
	                    
	# Assign the dimension data to the new NetCDF file.
	w_nc.variables['state'][:] = state
	
	#File dimensions for SPI 
	w_nc.createDimension('spivar', nspi)
	spivar = spival
	w_nc_nsyn = w_nc.createVariable('spivar', spivar.dtype,('spivar',))
	
	w_nc_nsyn.setncatts({'long_name': 'SPI3, SPI6, SPI12, SPI24'})
	                    
	# Assign the dimension data to the new NetCDF file.
	w_nc.variables['spivar'][:] = spivar
	
	# Assign the statistics
	#OBS SHAPE
	w_nc_var = w_nc.createVariable('shape', 'f', ('spivar','state', 'lon', 'lat'))
	
	long_name = "GEV shape parameter"
	         
	w_nc_var.setncatts({'long_name': long_name, 'shape_missing_value': obsmiss})
	w_nc.variables['shape'][:] = oshape
	
	#OBS LOC
	w_nc_var = w_nc.createVariable('loc', 'f', ('spivar','state', 'lon', 'lat'))
	
	long_name = "GEV location parameter"
	         
	w_nc_var.setncatts({'long_name': long_name, 'loc_missing_value': obsmiss})
	w_nc.variables['loc'][:] = oloc
	
	#OBS SCALE
	w_nc_var = w_nc.createVariable('scale', 'f', ('spivar','state', 'lon', 'lat'))
	
	long_name = "GEV scale parameter"
	         
	w_nc_var.setncatts({'long_name': long_name, 'scale_missing_value': obsmiss})
	w_nc.variables['scale'][:] = oscale
	
	#OBS MEAN
	w_nc_var = w_nc.createVariable('mean', 'f', ('spivar','state', 'lon', 'lat'))
	
	long_name = "mean"
	         
	w_nc_var.setncatts({'long_name': long_name, 'mean_missing_value': obsmiss})
	w_nc.variables['mean'][:] = omean
	
	#OBS MEDIAN
	w_nc_var = w_nc.createVariable('median', 'f', ('spivar','state', 'lon', 'lat'))
	
	long_name = "median"
	         
	w_nc_var.setncatts({'long_name': long_name, 'median_missing_value': obsmiss})
	w_nc.variables['median'][:] = omedian

	w_nc.close()
	
	return oshape, oloc, oscale, omean, omedian
	
	
#-----------------------------------------------------------------------------------

#SYNTH_SPI_CONFINT
'''
synth_spi_confint

Purpose:
---------
Computes the two-tailed 95% confidence intervals of the comparison metrics* from the synthetic
precipitation time series. 

*Comparison metrics:
-Perkins skill score (Sum(1,n) min(Z1, Z2)) where Zm, Zo are the relative frequencies from
datasets 1 and 2.
- Gamma distribution shape and scale parameters
- Mean and median of the data

Inputs:
---------

Outputs:
---------
Saved NetCDF files of the metrics for each of SPI3, 6, 12 and 24

History:
---------
20160524 - Created by Ailie Gallant, Monash University

'''
def synth_spi_confint(spival, nspistr):
	from scipy import stats
	import numpy as np
	from netcdf_tools import ncextractall
	from netCDF4 import Dataset
	from netCDF4 import date2num
	from netCDF4 import num2date
	import grid_tools
	import sys
	from timeseries_tools import moving_average
	from timeseries_tools import block_bootstrap
	from timeseries_tools import average_run
	from matplotlib import pyplot as plt
	from stat_metrics import perk_skillscore
	import time as comptime
	
	print("Starting calcs...")
	#Define synthetic data file for reading in later
	synfile = '/Users/ailieg/Data/drought_model_eval_data/analysis/'+nspistr+'.mon.GPCPsynthetic_landonly_2.nc'
	print(synfile)
	#Read in GPCP observations
	obfile = '/Users/ailieg/Data/drought_model_eval_data/data/obs/GPCP/precip.mon.mean.nc'
	obsnc = ncextractall(obfile)
	obsp = obsnc['precip']
	lon = obsnc['lon']
	nlon = lon.size
	lat = obsnc['lat']
	nlat = lat.size
	nsyn = 200
	time = obsnc['time']
	obsmiss = obsnc['precip_missing_value']
	obsp[np.where(obsp == obsmiss)] = np.nan
	print("Observed data has been extracted")
	
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
	
	#Make a land mask
	mask = grid_tools.make_land_mask(obfile)

	#Create arrays for saving the data. These are in the form
	#(state, 95CI, lon, lat) where 95CI is the bounds of the 
	#two-tailed 95% confidence interval from the synthetic data
	sshape = np.zeros((2,nsyn,nlat,nlon))+ obsmiss
	sloc = np.zeros((2,nsyn,nlat,nlon))+ obsmiss
	sscale = np.zeros((2,nsyn,nlat,nlon))+ obsmiss
	smean = np.zeros((2,nsyn,nlat,nlon))+ obsmiss
	smedian = np.zeros((2,nsyn,nlat,nlon))+ obsmiss
	spss = np.zeros((2,nsyn,nlat,nlon))+ obsmiss
	
	#Loop through LON and LAT (compute land only)
	for j in range(0,nlon):
		starttime = comptime.clock()
	
		for k in range(0,nlat):
			
			
			if mask[k,j] > 0.5:	
				#Extract the observed GPCP data for lat, lon
				savedo = obsp[:,k,j]
				
				#Smooth GPCP synthetic data as for the SPI
				tmpobs = moving_average(savedo, spival)
			
				#Extract the observed SPI data for GPCP for lat, lon
				tmpspiob = spi(savedo,spival)
				tmpspiob = tmpspiob
				tmpspiob = tmpspiob.flatten()
				
				#Use real data only
				tmpobs = tmpobs[~np.isnan(tmpobs)]
				tmpspiob = tmpspiob[~np.isnan(tmpspiob)]
				
				idx = np.arange(0,tmpobs.size)
				
				if tmpobs.size != tmpspiob.size:
					print(tmpobs.size, tmpspiob.size)
					print("Error at lon:", j, "lat:", k)
					sys.exit("Error: Observed PRCP & SPI not equal length!")
				
				#Compute shape, scale, location, mean and median metrics
				dmodobs = tmpobs[np.where(tmpspiob < -1.0)]
				modidx = idx[np.where(tmpspiob < -1.0)]
				nmod = dmodobs.size
				if nmod > 10: modblk = average_run(modidx)

				dsevobs = tmpobs[np.where(tmpspiob < -1.5)]
				sevidx = idx[np.where(tmpspiob < -1.5)]
				nsev = dsevobs.size
				if nsev > 10: sevblk = average_run(sevidx)
				
				#Turn on a flag if sample sizes are too small. This will return missing data.
				if nmod < 10: 
					modflag = 1
					print("Flag set for moderate drought, will return missing values")
				else:
					modflag = 0
					
				if nsev < 10: 
					sevflag = 1
					print("Flag set for severe drought, will return missing values")
				else:
					sevflag = 0
					
				
				#Loop through SYN series
				for m in range(0,nsyn):
					
					#Create synthetic series
					if modflag == 0:
						dmodsyn = block_bootstrap(dmodobs, modblk)
					else: dmodsyn = np.zeros(3)+obsmiss
					
					if sevflag == 0: 
						dsevsyn = block_bootstrap(dsevobs, sevblk)
					else: dsevsyn = np.zeros(3)+obsmiss
					
					
					#Do moderate droughts first, if samples sizes are < 20 then skip
					#iteration - the non-calculation will result in missing data
					if dmodsyn.size < 10:
						continue
						print("Mod. synth drought calculations skipped")
					else:
						try:
							gevmodpar = stats.genextreme.fit(dmodsyn)
						except ValueError:
							gevmodpar = np.zeros(3) + obsmiss
							print("ValueError encountered, missing values applied")
					
						sshape[0,m,k,j] = gevmodpar[0]
						sloc[0,m,k,j] = gevmodpar[1]
						sscale[0,m,k,j] = gevmodpar[2]
						smean[0,m,k,j] = np.mean(dmodsyn)
						smedian[0,m,k,j] = np.median(dmodsyn)
						
						#print("Saved:",gevmodpar[0], gevmodpar[1],gevmodpar[2], np.mean(dmodsyn), np.median(dmodsyn))
						#Perkins skill score
						allmoddata = np.concatenate((dmodobs, dmodsyn))
						nfrac = np.int((dmodobs.size)*0.1)+1
						if nfrac < 4: nfrac = 5
						minbin = np.min(allmoddata)
						maxbin = np.max(allmoddata)
						binw = (maxbin-minbin)/nfrac
						bins=np.linspace(minbin-binw, maxbin+binw, nfrac)
						spss[0,m,k,j] = perk_skillscore(dmodobs, dmodsyn, bins)
					
						
						
					#Do severe droughts second, if samples sizes are < 20 then skip
					#iteration - the non-calculation will result in missing data
					if dsevsyn.size < 10:
						continue
						print("Sev. synth drought calculations skipped")
					else:
						try:
							gevsevpar = stats.genextreme.fit(dsevsyn)
						except ValueError:
							gevsevpar = np.zeros(3) + obsmiss
							print("ValueError encountered, missing values applied")
				
						sshape[1,m,k,j] = gevsevpar[0]
						sloc[1,m,k,j] = gevsevpar[1]
						sscale[1,m,k,j] = gevsevpar[2]
						smean[1,m,k,j] = np.mean(dsevsyn)
						smedian[1,m,k,j] = np.median(dsevsyn)
				
						#Perkins Skill Score
						allsevdata = np.concatenate((dsevobs, dsevsyn))
						nfrac = np.int((dsevobs.size)*0.1)+1
						if nfrac < 4: nfrac = 5
						minbin = np.min(allsevdata)
						maxbin = np.max(allsevdata)
						binw = (maxbin-minbin)/nfrac
						bins=np.linspace(minbin-binw, maxbin+binw, nfrac)
						spss[1,m,k,j] = perk_skillscore(dsevobs, dsevsyn, bins)
						
			else: print("Ocean - calculations skipped")
			print("Calculating ", k+1, "of ", nlat)
		print("Time taken for iteration:",(comptime.clock() - starttime)/60.,"mins")
		print("Longitude ",j+1, "of", nlon, "complete")
		#your code here    
	
#Create a NetCDF file

	#Create NetCDF file to write
	outfile = synfile
	w_nc = Dataset(outfile, 'w', format='NETCDF4')

	#File description
	long_name = "Surrogates of the statistical metrics of GEV shape, location and scale ",\
	"parameters, mean, median and Perkins skill score."
	
	w_nc.description = long_name
	
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
	
	#File dimensions for STATE 
	w_nc.createDimension('state', 2)
	state = np.array(['Moderate', 'Severe'],dtype='str')
	w_nc_state = w_nc.createVariable('state', state.dtype,('state',))
	
	w_nc_state.setncatts({'long_name': 'Moderate drought (SPI < -1) or Severe drought (SPI < -1.5)'})
	                    
	# Assign the dimension data to the new NetCDF file.
	w_nc.variables['state'][:] = state
	
	#File dimensions for CI_BOUNDS 
	w_nc.createDimension('nsyn', nsyn)
	syn = np.arange(1,nsyn+1)
	w_nc_nsyn = w_nc.createVariable('nsyn', syn.dtype,('nsyn',))
	
	w_nc_nsyn.setncatts({'long_name': 'Number of surrogates'})
	                    
	# Assign the dimension data to the new NetCDF file.
	w_nc.variables['nsyn'][:] = syn


	#SURROGATE SHAPE
	w_nc_var = w_nc.createVariable('shape', 'f', ('state','nsyn', 'lat', 'lon'))
	
	long_name = "Shape parameters of surrogate series"
	         
	w_nc_var.setncatts({'long_name': long_name, 'shape_missing_value': obsmiss})
	w_nc.variables['shape'][:] = sshape
	
	#SURROGATE LOC
	w_nc_var = w_nc.createVariable('loc', 'f', ('state','nsyn', 'lat', 'lon'))
	
	long_name = "Location parameters of surrogate series"
	         
	w_nc_var.setncatts({'long_name': long_name, 'loc_missing_value': obsmiss})
	w_nc.variables['loc'][:] = sloc
	
	#SURROGATE SCALE
	w_nc_var = w_nc.createVariable('scale', 'f', ('state','nsyn', 'lat', 'lon'))
	
	long_name = "Scale parameters of surrogate series"
	         
	w_nc_var.setncatts({'long_name': long_name, 'scale_missing_value': obsmiss})
	w_nc.variables['scale'][:] = sscale
	
	#SURROGATE MEAN
	w_nc_var = w_nc.createVariable('mean', 'f', ('state','nsyn', 'lat', 'lon'))
	
	long_name = "Mean parameters of surrogate series"
	         
	w_nc_var.setncatts({'long_name': long_name, 'mean_missing_value': obsmiss})
	w_nc.variables['mean'][:] = smean
	
	#SURROGATE MEDIAN
	w_nc_var = w_nc.createVariable('median', 'f', ('state','nsyn', 'lat', 'lon'))
	
	long_name = "Median parameters of surrogate series"
	         
	w_nc_var.setncatts({'long_name': long_name, 'median_missing_value': obsmiss})
	w_nc.variables['median'][:] = smedian
	
	#SURROGATE Perkins Skill Score (PSS)
	w_nc_var = w_nc.createVariable('pss', 'f', ('state','nsyn', 'lat', 'lon'))
	
	long_name = "Perkins Skill Score of surrogate series"
	         
	w_nc_var.setncatts({'long_name': long_name, 'pss_missing_value': obsmiss})
	w_nc.variables['pss'][:] = spss

	w_nc.close()
	
	return 
	
	
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
#SPI_OBS_MODEL_COMPARE
'''
spi_obs_model_compare

Purpose:
---------
Compares the metrics from the CMIP5 or AMIP simulations to the synthetic 95th percentile
confidence interval.

*Comparison metrics:
-Perkins skill score (Sum(1,n) min(Z1, Z2)) where Zm, Zo are the relative frequencies from
datasets 1 and 2.
- Gamma distribution shape and scale parameters
- Mean and median of the data

Inputs:
---------

Outputs:
---------
Saved NetCDF files of the metrics for each of SPI3, 6, 12 and 24

History:
---------
20160524 - Created by Ailie Gallant, Monash University

'''
def spi_obs_model_compare(spiname, metname, state, modfile):
	from netcdf_tools import ncextractall
	from netCDF4 import Dataset
	from netCDF4 import date2num
	from netCDF4 import num2date
	import time as comptime
	import numpy as np
	from grid_tools import modgrid_to_obsgrid
	from matplotlib import pyplot as plt
	from grid_tools import make_land_mask
	
	if state == 'mod': sind = 0
	if state == 'sev': sind = 1
	
	#Metrics for the surrogate data set
	print("Reading metrics from bootstrapped data...")
	surrmetname = '/Users/ailieg/Data/drought_model_eval_data/analysis/'+spiname+'.mon.GPCPsynthetic_landonly.nc'
	surrnc = ncextractall(surrmetname)
	surrlat = surrnc['lat']
	surrlon = surrnc['lon']
	surrmet = surrnc[metname]
	surrmetone = surrmet[sind,:,:,:]

	print("Reading metrics from bootstrapped data...")
	surrmetname = '/Users/ailieg/Data/drought_model_eval_data/analysis/'+spiname+'.mon.GPCPsynthetic_landonly_2.nc'
	surrnc = ncextractall(surrmetname)
	surrmet = surrnc[metname]
	surrmettwo = surrmet[sind,:,:,:]
	
	surrmet = np.zeros((300,surrlat.size,surrlon.size))
	surrmet[0:100,:,:] = surrmetone
	surrmet[100:,:,:] = surrmettwo
	
	if spiname == 'SPI3': spiind = 0
	if spiname == 'SPI6': spiind = 1
	if spiname == 'SPI12': spiind = 2
	if spiname == 'SPI24': spiind = 3

	#Metrics for the observed metrics
	print("Reading metrics from observed data...")
	obmetname = '/Users/ailieg/Data/drought_model_eval_data/analysis/spi_obs_metrics_nomask.nc'
	obnc = ncextractall(obmetname)
	oblat = obnc['lat']
	oblat = oblat[::-1]
	oblon = obnc['lon']
	obmet = obnc[metname]
	obmet = (obmet[spiind,sind,:,:]).T
	
	#Read in model data
	file = '/Users/ailieg/Data/drought_model_eval_data/analysis/'+modfile
	ncmod = ncextractall(modfile)
	modmet = ncmod[metname]
	modlon = ncmod['lon']
	modlat = ncmod['lat']
	modmet = modmet[sind,:,:]
	
	modmet = modgrid_to_obsgrid(modmet, modlat, modlon, oblat, oblon)
	
	miss = -9e35
	sigdif = modmet - obmet
	sigdif[obmet <= miss] = miss
	sigdif[modmet <= miss] = miss
	
	surrmin = np.min(surrmet, axis=0)
	surrmax = np.max(surrmet, axis=0)
	
	sigdif[modmet < surrmin] = miss
	sigdif[modmet > surrmax] = miss
	
		#Make a land mask
	obfile = '/Users/ailieg/Data/drought_model_eval_data/data/obs/GPCP/precip.mon.mean.nc'
	mask = make_land_mask(obfile)
	sigdif[(sigdif != miss)] = 1.0
	sigdif[(mask < 1)] = miss
	sigdif[(sigdif == miss)] = np.nan

	obmet[np.where(obmet < miss)] = np.nan
	dif = modmet - obmet
	

	return dif, sigdif, oblat, oblon

#--------------------------------------------------------------------------------------
#READ_SPI_SURR_METRICS
'''
read_spi_surr_metrics

Purpose:
---------
Reads in the metrics

*Comparison metrics:
-Perkins skill score (Sum(1,n) min(Z1, Z2)) where Zm, Zo are the relative frequencies from
datasets 1 and 2.
- Gamma distribution shape and scale parameters
- Mean and median of the data

Outputs:
---------
Saved NetCDF files of the metrics for each of SPI3, 6, 12 and 24

History:
---------
20160524 - Created by Ailie Gallant, Monash University

'''

def read_spi_surr_metrics(metricname):
	import numpy as np
	from netCDF4 import Dataset
	from netcdf_tools import ncextractall
	
	
	#Read in metrics
	path = '/Users/ailieg/Data/drought_model_eval_data/analysis/'
	startfile = 'SPI3_SPI12_statbounds_lon0_1.nc'
	
	#Read in GPCP observations
	file = path+startfile
	nc = ncextractall(file)
	state = nc['state']
	nstate = state.size
	spival = nc['spivar']
	nspi = 2
	ci = nc['ci_bounds']
	nci = ci.size

	obfile = '/Users/ailieg/Data/drought_model_eval_data/data/obs/GPCP/precip.mon.mean.nc'
	obsnc = ncextractall(obfile)
	lon = obsnc['lon']
	nlon = lon.size
	lat = obsnc['lat']
	nlat = lat.size
	
	if metricname[0] == 'o':
		metric = np.zeros((nspi, nstate, nlon, nlat))
	else:
		metric = np.zeros((nspi, nstate, nci, nlon, nlat))
	
	for i in range(0,nlon-2,2):
		
		if i < 140:
			
			if i < 87:
				tmpfile = 'SPI3_SPI12_statbounds_lon'+str(i)+'_'+str(i+1)+'.nc'
				nc = ncextractall(tmpfile)
		
				if metricname[0] == 'o':
					metric[:,:,i:i+2,:] = nc[metricname]
				else:
					metric[:,:,:,i:i+2,:] = nc[metricname]
			else:
				tmpfile = 'SPI3_SPI12_statbounds_lon'+str(i+1)+'_'+str(i+2)+'.nc'
				print(tmpfile)
				nc = ncextractall(tmpfile)
		
				if metricname[0] == 'o':
					metric[:,:,i+1:i+3,:] = nc[metricname]
				else:
					metric[:,:,:,i+1:i+3,:] = nc[metricname]
		else:
			
			tmpfile = 'SPI3_SPI12_statbounds_lon'+str(i+1)+'_'+str(i+3)+'.nc'
			nc = ncextractall(tmpfile)

			if metricname[0] == 'o':
				metric[:,:,i+1:i+3,:] = nc[metricname]
			else:
				metric[:,:,:,i+1:i+3,:] = nc[metricname]
		
	return metric, lat, lon

#------------------------------------------------------------------------------------
'''
compute_spi_model_metrics

Purpose:
---------
Computes the metrics for the CMIP5 and AMIP models

*Comparison metrics:
-Perkins skill score (Sum(1,n) min(Z1, Z2)) where Zm, Zo are the relative frequencies from
datasets 1 and 2.
- Gamma distribution shape and scale parameters
- Mean and median of the data

Outputs:
---------
Saved NetCDF files of the metrics for each of SPI3, 6, 12 and 24

History:
---------
20160524 - Created by Ailie Gallant, Monash University

'''
def compute_spi_model_metrics(modfile, modspifile, modname, spival,spiname):

	from scipy import stats
	import numpy as np
	from netcdf_tools import ncextractall
	from netCDF4 import Dataset
	from netCDF4 import date2num
	from netCDF4 import num2date
	import grid_tools
	import sys
	from timeseries_tools import moving_average
	from stat_metrics import perk_skillscore
	
	print("Reading model data...")
	
	#Read model precip data file
	datdict = ncextractall(modfile)
	
	#Read each variable
	data = datdict['pr']
	lon = datdict['lon']
	nlon = lon.size
	lat = datdict['lat']
	nlat = lat.size
	time = datdict['time']
	
	#convert time units to actual date
	time_u = datdict['time_units']
	if 'time_calendar' in datdict.keys(): 
	     cal = datdict['time_calendar']
	     time = num2date(time,units = time_u, calendar=cal)
	else: time = num2date(time,units = time_u)
	
	#convert missing numbers to NaN
	miss = datdict["pr_missing_value"] 
	data = np.where(data == miss,np.nan,data)
	
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
	ntime = time.size
	
	#Convert model data to same units as obs
	data = data * 86400.
	
	#Create empty arrays to save metrics for the obs and synthetic data
	#Arrays are of the form (state,lon,lat), where state indicates moderate
	#or severe drought.
	mshape = np.zeros((2,nlat,nlon))+ np.nan
	mloc = np.zeros((2,nlat,nlon))+ np.nan
	mscale = np.zeros((2,nlat,nlon))+ np.nan
	mmean = np.zeros((2,nlat,nlon))+ np.nan
	mmedian = np.zeros((2,nlat,nlon))+ np.nan
	
	#Loop through LON and LAT (compute land only)
	for j in range(0,nlon):
		print("Now running calc for longitude ",j, "of",nlon)
		for k in range(0,nlat):
			
			#Extract the observed GPCP data for lat, lon
			tmpdata = data[:,k,j]
			tmpdata = tmpdata.flatten()
			
			tmpspi = spi(tmpdata,spival)
			tmpspi = tmpspi.flatten()
			
			#Smooth GPCP synthetic data as for the SPI
			tmpdata = moving_average(tmpdata, spival)
			
			#Use real data only
			tmpdata = tmpdata[~np.isnan(tmpdata)]
			tmpspi = tmpspi[~np.isnan(tmpspi)]
			
			if tmpdata.size != tmpspi.size:
				print("Error at lon:", k, "lat:", j)
				sys.exit("Error: Model PRCP & SPI not equal length!")
			
			#Compute shape, scale, location, mean and median metrics
			dmod = tmpdata[np.where(tmpspi < -1.0)]
			nmod = dmod.size
			dsev = tmpdata[np.where(tmpspi < -1.5)]
			nsev = dsev.size
			
			if nmod > 10:
				try:
					gevmodpar = stats.genextreme.fit(dmod)
				except ValueError:
					gevmodpar = np.zeros(3) + np.nan
					print("ValueError encountered, missing values applied ")
				
				mshape[0,k,j] = gevmodpar[0]
				mloc[0,k,j] = gevmodpar[1]
				mscale[0,k,j] = gevmodpar[2]
				mmean[0,k,j] = np.mean(dmod)
				mmedian[0,k,j] = np.median(dmod)
			else:
				continue
				
			if nsev > 10:
				try:
					gevsevpar = stats.genextreme.fit(dsev)
				except ValueError:
					gevsevpar = np.zeros(3) + np.nan
					print("ValueError encountered, missing values applied")
				
				mshape[1,k,j] = gevsevpar[0]
				mloc[1,k,j] = gevsevpar[1]
				mscale[1,k,j] = gevsevpar[2]
				mmean[1,k,j] = np.mean(dsev)
				mmedian[1,k,j] = np.median(dsev)
	
	#Convert missing values back to the miss data
	mshape[np.isnan(mshape)] = miss
	mloc[np.isnan(mloc)] = miss
	mscale[np.isnan(mscale)] = miss
	mmedian[np.isnan(mmedian)] = miss
	mmean[np.isnan(mmean)] = miss
	
	#Create a NetCDF file
	state = [0,1]
	
	#Create NetCDF file to write
	outpath = '/Users/ailieg/Data/drought_model_eval_data/analysis/'
	outfile = spiname+'_'+modname+'_stats_.nc'
	w_nc = Dataset(outfile, 'w', format='NETCDF4')

	#File description
	long_name = "The statistical metrics of GEV shape, location and scale parameters, mean and median",\
	             "Moderate and severe drought was defined as below an SPI of -1 and -1.5 respectively. "
	w_nc.description = long_name
	
	#File dimensions for LAT
	w_nc.createDimension('lat', nlat)
	w_nc_lat = w_nc.createVariable('lat', lat.dtype,('lat',))
	
	w_nc_lat.setncatts({'long_name': 'Latitude',\
	                    'units': 'Degrees North'})
	                    
	# Assign the dimension data to the new NetCDF file.
	w_nc.variables['lat'][:] = lat
	     
	#File dimensions for LON
	w_nc.createDimension('lon', nlon)
	w_nc_lon = w_nc.createVariable('lon', lon.dtype,('lon',))
	
	w_nc_lon.setncatts({'long_name': 'Longitude',\
	                    'units': 'Degrees East'})
	                    
	# Assign the dimension data to the new NetCDF file.
	w_nc.variables['lon'][:] = lon
	
	#File dimensions for STATE 
	w_nc.createDimension('state', 2)
	state = np.array(['Moderate', 'Severe'],dtype='str')
	w_nc_nsyn = w_nc.createVariable('state', state.dtype,('state',))
	
	w_nc_nsyn.setncatts({'long_name': 'Moderate drought (SPI < -1) or Severe drought (SPI < -1.5)'})
	                    
	# Assign the dimension data to the new NetCDF file.
	w_nc.variables['state'][:] = state
	
	# Assign the statistics
	#OBS SHAPE
	w_nc_var = w_nc.createVariable('shape', 'f', ('state', 'lat', 'lon'))
	
	long_name = "GEV shape parameter"
	         
	w_nc_var.setncatts({'long_name': long_name, 'shape_missing_value': miss})
	w_nc.variables['shape'][:] = mshape
	
	#OBS LOC
	w_nc_var = w_nc.createVariable('loc', 'f', ('state', 'lat', 'lon'))
	
	long_name = "GEV location parameter"
	         
	w_nc_var.setncatts({'long_name': long_name, 'loc_missing_value': miss})
	w_nc.variables['loc'][:] = mloc
	
	#OBS SCALE
	w_nc_var = w_nc.createVariable('scale', 'f', ('state', 'lat', 'lon'))
	
	long_name = "GEV scale parameter"
	         
	w_nc_var.setncatts({'long_name': long_name, 'scale_missing_value': miss})
	w_nc.variables['scale'][:] = mscale
	
	#OBS MEAN
	w_nc_var = w_nc.createVariable('mean', 'f', ('state', 'lat', 'lon'))
	
	long_name = "mean"
	         
	w_nc_var.setncatts({'long_name': long_name, 'mean_missing_value': miss})
	w_nc.variables['mean'][:] = mmean
	
	#OBS MEDIAN
	w_nc_var = w_nc.createVariable('median', 'f', ('state', 'lat', 'lon'))
	
	long_name = "median"
	         
	w_nc_var.setncatts({'long_name': long_name, 'median_missing_value': miss})
	w_nc.variables['median'][:] = mmedian

	w_nc.close()
	return 

#---------------------------------------------------------------------------
'''
batch_spi_model_metrics

Purpose:
---------
Computes the metrics for the CMIP5 and AMIP models

*Comparison metrics:
-Perkins skill score (Sum(1,n) min(Z1, Z2)) where Zm, Zo are the relative frequencies from
datasets 1 and 2.
- Gamma distribution shape and scale parameters
- Mean and median of the data

Outputs:
---------
Saved NetCDF files of the metrics for each of SPI3, 6, 12 and 24

History:
---------
20160524 - Created by Ailie Gallant, Monash University

'''
def batch_spi_model_metrics(nspi):

	spi = 'SPI'+str(nspi)
	
	amippath = '/Users/ailieg/Data/drought_model_eval_data/data/AMIP/'
	
	prfile = ['FGOALS-s2/r1i1p1/pr/pr_Amon_FGOALS-s2_amip_r1i1p1_197901-200812.nc',\
	'GFDL-CM3/r1i1p1/pr/pr_Amon_GFDL-CM3_amip_r1i1p1_197901_200812.nc',\
	'HadGEM2-A/r1i1p1/pr/pr_Amon_HadGEM2-A_amip_r1i1p1_197809-200811.nc',\
	'NorESM1-M/r1i1p1/pr/pr_Amon_NorESM1-M_amip_r1i1p1_197901-200512.nc']
	
	spifile = ['FGOALS-s2/r1i1p1/'+spi+'_FGOALS-s2_amip_r1i1p1_197901-200812.nc',\
	'GFDL-CM3/r1i1p1/'+spi+'_GFDL-CM3_amip_r1i1p1_197901_200812.nc',\
	'HadGEM2-A/r1i1p1/'+spi+'_HadGEM2-A_amip_r1i1p1_197809-200811',\
	'NorESM1-M/r1i1p1/'+spi+'_NorESM1-M_amip_r1i1p1_197901-200512.nc']
	
	modname = ['FGOALS-s2_amip_r1i1p1','GFDL-CM3_amip_r1i1p1','HadGEM2-A_amip_r1i1p1','NorESM1-M_amip_r1i1p1']
	
	for i in range(0,4):
		modfile = amippath+prfile[i]
		spiin = amippath+spifile[i]
		mname = modname[i]
		comp = compute_spi_model_metrics(modfile, spiin, mname,nspi,spi)
	
	cmippath = '/Users/ailieg/Data/drought_model_eval_data/data/CMIP5/'
	
	prfile = ['ACCESS1-0/r1i1p1/pr/pr_Amon_ACCESS1-0_historical_r1i1p1_185001-200512.nc',\
	'CanESM2/r1i1p1/pr/pr_Amon_CanESM2_historical_r1i1p1_185001-200512.nc',\
	'GFDL-CM3/r1i1p1/pr/pr_Amon_GFDL-CM3_historical_r1i1p1_186001-200512.nc',\
	'HadGEM2-CC/r1i1p1/pr/pr_Amon_HadGEM2-CC_historical_r1i1p1_185912-200511.nc',\
	'MPI-ESM-P/r1i1p1/pr/pr_Amon_MPI-ESM-P_historical_r1i1p1_185001-200512.nc',\
	'CCSM4/r1i1p1/pr/pr_Amon_CCSM4_historical_r1i1p1_185001-200512.nc',\
	'FGOALS-s2/r1i1p1/pr/pr_Amon_FGOALS-s2_historical_r1i1p1_185001-200512.nc',\
	'GISS-E2-R/r6i1p1/pr/pr_Amon_GISS-E2-R_historical_r6i1p1_185001-200512.nc',\
	'IPSL-CM5B-LR/r1i1p1/pr/pr_Amon_IPSL-CM5B-LR_historical_r1i1p1_185001-200512.nc',\
	'NorESM1-M/r1i1p1/pr/pr_Amon_NorESM1-M_historical_r1i1p1_185001-200512.nc']
	
	spifile = ['ACCESS1-0/r1i1p1/'+spi+'_ACCESS1-0_historical_r1i1p1_185001-200512.nc',\
	'CanESM2/r1i1p1/'+spi+'_CanESM2_historical_r1i1p1_185001-200512.nc',\
	'GFDL-CM3/r1i1p1/'+spi+'_GFDL-CM3_historical_r1i1p1_186001-200512.nc',\
	'HadGEM2-CC/r1i1p1/'+spi+'_HadGEM2-CC_historical_r1i1p1_185912-200511.nc',\
	'MPI-ESM-P/r1i1p1/'+spi+'_MPI-ESM-P_historical_r1i1p1_185001-200512.nc',\
	'CCSM4/r1i1p1/'+spi+'_CCSM4_historical_r1i1p1_185001-200512.nc',\
	'FGOALS-s2/r1i1p1/'+spi+'_FGOALS-s2_historical_r1i1p1_185001-200512.nc',\
	'GISS-E2-R/r6i1p1/'+spi+'_GISS-E2-R_historical_r6i1p1_185001-200512.nc',\
	'IPSL-CM5B-LR/r1i1p1/'+spi+'_IPSL-CM5B-LR_historical_r1i1p1_185001-200512.nc',\
	'NorESM1-M/r1i1p1/'+spi+'_NorESM1-M_historical_r1i1p1_185001-200512.nc']
	
	modname = ['ACCESS1-0_historical_r1i1p1','CanESM2_historical_r1i1p1', 'GFDL-CM3_historical_r1i1p1',\
	'HadGEM2-CC_historical_r1i1p1','MPI-ESM-P_historical_r1i1p1','CCSM4_historical_r1i1p1',\
	'FGOALS-s2_historical_r1i1p1','GISS-E2-R_historical_r6i1p1','IPSL-CM5B-LR_historical_r1i1p1',\
	'NorESM1-M_historical_r1i1p1']
	
	for i in range(0,10):
		modfile = cmippath+prfile[i]
		spiin = cmippath+spifile[i]
		mname = modname[i]
		print(mname)
		comp = compute_spi_model_metrics(modfile, spiin, mname,nspi,spi)

	return
	
#-------------------------------------------------------------------------------------
