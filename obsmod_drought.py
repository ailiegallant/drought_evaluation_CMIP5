#THIS MODULE HAS SCRIPTS WHICH COMPARE OBSERVATIONS AND MODELS FOR THE 
#2016 DROUGHT MODEL EVALUATION PAPER
#
#EACH SUBROUTINE IS SEPARATED BY A LINE --------



#-----------------------------------------------------------------------------
'''
RCS_GPCP

This script reads (R), clips (C) and smooths (S) the GPCP data.

INPUT:
  winlen - the length of the period over which to average (e.g. 3 months, winlen = 3)
  

Created by Ailie Gallant 20/07/2016
'''
def rcs_gpcp(winlen):
	from grid_tools import trim_time_jandec
	from netCDF4 import num2date
	from netCDF4 import date2num
	from scipy import ndimage
	from netcdf_tools import ncextractall
	from convert import mmd_mmm
	
	#Extract the observed data and clip to the required start and end months
	obfile = '/Users/ailieg/Data/drought_model_eval_data/data/obs/GPCP/precip.mon.mean.nc'
	obsnc = ncextractall(obfile)
	odata = obsnc['precip']
	olon = obsnc['lon']
	olat = obsnc['lat']
	olat = olat[::-1]
	odata = odata[:,::-1,:]
	otime = obsnc['time']
	obsmiss = obsnc['precip_missing_value']
	odata[np.where(odata == obsmiss)] = np.nan
	
	time_u = obsnc['time_units']
	if 'time_calendar' in obsnc.keys(): 
	     cal = obsnc['time_calendar']
	     otime = num2date(otime,units = time_u, calendar=cal)
	else: otime = num2date(otime,units = time_u)
	
	odata, otime = trim_time_jandec(odata, otime)

	odata = mmd_mmm(odata)
	odata = ndimage.filters.uniform_filter(odata,size=[winlen,1,1])
	
		#Trim first or last values if required as they are unrepresentative
	trim = int(winlen/2)
	odata = odata[trim:,:,:]
	if winlen % 2 == 0: trim = trim - 1
	odata = odata[:-trim,:,:]

	return, odata, olat, olon

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
'''
RCS_MODEL

This script reads (R), clips (C) and smooths (S) CMIP5 or AMIP data.
 

INPUT:
  winlen - the length of the period over which to average (e.g. 3 months, winlen = 3)
  

Created by Ailie Gallant 20/07/2016
'''
def rcs_model(winlen, modfile):
	from grid_tools import trim_time_jandec
	from netCDF4 import num2date
	from netCDF4 import date2num
	from scipy import ndimage
	from netcdf_tools import ncextractall
	from convert import mmd_mmm

	#Extract the model data and clip to the required start and end months
	modnc = ncextractall(modfile)
	mdata = modnc['pr']
	mdata = mdata*86400. #convert to same units as obs
	mlon = modnc['lon']
	mlat = modnc['lat']
	mtime = modnc['time']
	
	time_u = modnc['time_units']
	if 'time_calendar' in modnc.keys(): 
	     cal = modnc['time_calendar']
	     mtime = num2date(mtime,units = time_u, calendar=cal)
	else: mtime = num2date(mtime,units = time_u)
	
	mdata, mtime = trim_time_jandec(mdata, mtime)
	
	mdata = mmd_mmm(mdata)
	mdata = ndimage.filters.uniform_filter(mdata,size=[winlen,1,1])
	
		#Trim first or last values if required as they are unrepresentative
	trim = int(winlen/2)
	mdata = mdata[trim:,:,:]
	if winlen % 2 == 0: trim = trim - 1
	mdata = mdata[:-trim,:,:]
	
	return, mdata, mlat, mlon

#-----------------------------------------------------------------------------
'''
PERC_COMPARE

This script compares the percentile threshold value from observations with model data.
The input threshold is provided and raw precipitation data are compared as the difference
between the two (model - obs). A positive value indicates the models overestimate the
threshold, a negative value indicates the models underestimate the threshold.

INPUT:
  odata - a gridded data set of the observations to be compared*
  mdata - a gridded data set of the model to be compared*
  pc - the percentile threshold to compare
  winlen - the length of the period over which to average (e.g. 3 months, winlen = 3)
  season - the starting month of the season of length winlen to examine, e.g. if season = 1
           and winlen = 12 this defines a january-december year. 
  
*Inputs should be in the form of [time, lat, lon]

Created by Ailie Gallant 20/07/2016
'''

def perc_compare(pc, winlen, season, obfile, modfile, intitle):
	import numpy as np
	from grid_tools import interp_togrid
	
	#Extract model and obs data
	odata, olat, olon = rcs_gpcp(winlen)
	mdata, mlat, mlon = rcs_gpcp(winlen, modfile)
	
	#Shift data to start of required season
	rlen = -1*(season - 1)
	odata = np.roll(odata, rlen, axis=0)
	odata = odata[0::12, :, :]
	mdata = np.roll(mdata, rlen, axis=0)
	mdata = mdata[0::12, :, :]
	
	#Calculate thresholds
	opc = np.percentile(odata, pc, axis=0)
	opc[opc < 0.0] = 0.0 #ensure percentile is not below zero, if below zero then set to 0
	mpc = np.percentile(mdata, pc, axis=0)
	mpc[mpc < 0.0] = 0.0 #ensure percentile is not below zero, if below zero then set to 0
	
	#Interpolate to the same grid spacing (coarsest) and compute difference as a percentage
	ob, mod, lon, lat = interp_togrid(opc, olat, olon, mpc, mlat, mlon)
	d = ((mod - ob)/ob)*100
	
	plot = plot_perc_compare (d, lat, lon, pc, winlen, season, intitle)

	return d, lat, lon
	
#-----------------------------------------------------------------------------	
'''
PLOT_PERC_COMPARE

This script plots the comparison of the percentile thresholds

INPUT:
  data - a gridded data set of the observations to be compared* 
  lat - an array of latitudes associated with the data
  lon - an array of longitudes associated with the data
  season - 
  
*Inputs should be in the form of [time, lat, lon]

Created by Ailie Gallant 20/07/2016
'''
def plot_perc_compare (d, lat, lon, pc, winlen, season, intitle):
	from matplotlib import pyplot as plt
	import cartopy.crs as ccrs
	from matplotlib.colors import BoundaryNorm
	from matplotlib.ticker import MaxNLocator
	import numpy as np
	from grid_tools import fold_grid
	
	#Set levels and colormap
	levels = np.arange(-200,220,20)
	cmap = plt.get_cmap('Spectral')
	norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
	
	#Make the grid circular
	d, lat, lon = fold_grid(d, lat, lon)
	
	#Set axes and plot
	ax = plt.axes(projection=ccrs.PlateCarree())
	p=plt.pcolormesh(lon, lat, d, cmap=cmap, norm=norm)
	
	#Add a colorbar
	cbar = plt.colorbar(p, extend='both')
	cbar.ax.set_ylabel('%')
	ax.coastlines()
	
	#Create title for saved plot
	seasname = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec',\
	'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec',\
	'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
	title = intitle+'_'+str(winlen)+'mth_'+seasname[season-1]+seasname[season+winlen-2]+'_'+str(pc)+'th%ile'
	plt.title(title, fontsize=10)
	
	#Save the output
	savefile = title+'.png'
	outfile = '/Users/ailieg/Data/drought_model_eval_data/analysis/'+savefile
	plt.savefig(outfile, dpi=400, format='png',bbox_inches='tight')
	plt.close()
	
#----------------------------------------------------------------------------------	
'''
PERC_COMPARE_BSMODLEN

This script compares the percentile threshold value from observations with model data.
The input threshold is provided and raw precipitation data are compared as the difference
between the two (model - obs) as a percentage. A positive value indicates the models overestimate the
threshold, a negative value indicates the models underestimate the threshold. The significance
of the difference between the two is computed as a bootstrap of resampled observations of the same
length as the model time series OR a clipped model time series the same length as the 
observations.

INPUT:
  odata - a gridded data set of the observations to be compared*
  mdata - a gridded data set of the model to be compared*
  pc - the percentile threshold to compare
  winlen - the length of the period over which to average (e.g. 3 months, winlen = 3)
  season - the starting month of the season of length winlen to examine, e.g. if season = 1
           and winlen = 12 this defines a january-december year.
  
*Inputs should be in the form of [time, lat, lon]

Created by Ailie Gallant 26/07/2016
'''
def perc_compare_bsmodlen(pc, winlen, season, obfile, modfile, intitle):
	import numpy as np
	from grid_tools import interp_togrid
	
	#Extract model and obs data
	odata, olat, olon = rcs_gpcp(winlen)
	mdata, mlat, mlon = rcs_gpcp(winlen, modfile)
	
	#Shift data to start of required season and trim
	rlen = -1*(season - 1)
	odata = np.roll(odata, rlen, axis=0)
	odata = odata[0::12, :, :]
	mdata = np.roll(mdata, rlen, axis=0)
	mdata = mdata[0::12, :, :]
	
	#Calculate thresholds
	#Set the length of the bootstrapped data set
	n = mdata.shape[0]  #length is the length of the model data
	print("Commencing bootstrap...")
	opc, bspc = bootstrap_perc_obs(odata, n, pc)
	print("Bootstrap completed...")
	
	mpc = np.percentile(mdata, pc, axis=0)
	mpc[mpc < 0.0] = 0.0 #ensure percentile is not below zero, if below zero then set to 0
	
	#Interpolate to the same grid spacing (coarsest) and compute difference as a percentage
	ob, mod, lon, lat = interp_togrid(opc, olat, olon, mpc, mlat, mlon)
	d = ((mod - ob)/ob)*100
	
	#Calculate significance
	bslow = np.percentile(bspc, 2.5, axis=0)
	bshi = np.percentile(bspc, 97.5, axis=0)
	
	sig = np.zeros((lat.size, lon.size))
	sig[np.logical_and(mod >= bslow, mod <= bshi)] = 1
	
	plot = plot_perc_compare_bs (d, sig, lat, lon, pc, winlen, season, intitle)

	return d, sig, lat, lon

#----------------------------------------------------------------------------------	
'''
PERC_COMPARE_BSOBLEN

This script compares the percentile threshold value from observations with model data.
The input threshold is provided and raw precipitation data are compared as the difference
between the two (model - obs) as a percentage. A positive value indicates the models overestimate the
threshold, a negative value indicates the models underestimate the threshold. The significance
of the difference between the two is computed as a bootstrap of resampled observations of the same
length as the model time series OR a clipped model time series the same length as the 
observations.

INPUT:
  odata - a gridded data set of the observations to be compared*
  mdata - a gridded data set of the model to be compared*
  pc - the percentile threshold to compare
  winlen - the length of the period over which to average (e.g. 3 months, winlen = 3)
  season - the starting month of the season of length winlen to examine, e.g. if season = 1
           and winlen = 12 this defines a january-december year.
  
*Inputs should be in the form of [time, lat, lon]

Created by Ailie Gallant 26/07/2016
'''
def perc_compare_bsoblen(pc, winlen, season, obfile, modfile, intitle):
	import numpy as np
	from grid_tools import interp_togrid
	
	#Extract model and obs data
	odata, olat, olon = rcs_gpcp(winlen)
	mdata, mlat, mlon = rcs_gpcp(winlen, modfile)
	
	#Shift data to start of required season and trim
	rlen = -1*(season - 1)
	odata = np.roll(odata, rlen, axis=0)
	odata = odata[0::12, :, :]
	mdata = np.roll(mdata, rlen, axis=0)
	mdata = mdata[0::12, :, :]
	
	#Calculate thresholds
	#Set the length of the bootstrapped data set
	n = mdata.shape[0]  #length is the length of the model data
	print("Commencing bootstrap...")
	opc, bspc = bootstrap_perc_obs(odata, n, pc)
	print("Bootstrap completed...")
	
	mpc = np.percentile(mdata, pc, axis=0)
	mpc[mpc < 0.0] = 0.0 #ensure percentile is not below zero, if below zero then set to 0
	
	#Interpolate to the same grid spacing (coarsest) and compute difference as a percentage
	ob, mod, lon, lat = interp_togrid(opc, olat, olon, mpc, mlat, mlon)
	d = ((mod - ob)/ob)*100
	
	#Calculate significance
	bslow = np.percentile(bspc, 2.5, axis=0)
	bshi = np.percentile(bspc, 97.5, axis=0)
	
	sig = np.zeros((lat.size, lon.size))-1
	sig[np.logical_and(mod >= bslow, mod <= bshi)] = 1
	
	plot = plot_perc_compare_bs (d, sig, lat, lon, pc, winlen, season, intitle)

	return d, sig, lat, lon




#---------------------------------------------------------------------------------
'''
PLOT_PERC_COMPARE_BS

This script plots the comparison of the percentile thresholds including bootstrapped signif

INPUT:
  data - a gridded data set of the observations to be compared* 
  lat - an array of latitudes associated with the data
  lon - an array of longitudes associated with the data
  season - 
  
*Inputs should be in the form of [time, lat, lon]

Created by Ailie Gallant 26/07/2016
'''

def plot_perc_compare_bs (d, sig, lat, lon, pc, winlen, season, intitle):
	from matplotlib import pyplot as plt
	import cartopy.crs as ccrs
	from matplotlib.colors import BoundaryNorm
	from matplotlib.ticker import MaxNLocator
	import numpy as np
	from grid_tools import fold_grid
	
	#Set levels and colormap
	levels = np.arange(-200,220,20)
	cmap = plt.get_cmap('Spectral')
	norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
	
	#Make the grid circular
	d, lat, lon = fold_grid(d, lat, lon)
	lon = lon[0:-1]
	sig, lat, lon = fold_grid(sig, lat, lon)
	
	#Set axes and plot
	ax = plt.axes(projection=ccrs.PlateCarree())
	p=plt.pcolormesh(lon, lat, d, cmap=cmap, norm=norm)
	cont = plt.contourf(lon, lat, sig, levels=[-1,0,1], transform = ccrs.PlateCarree(),colors='none', hatches=[" ","."])
	
	#Add a colorbar
	cbar = plt.colorbar(p, extend='both')
	cbar.ax.set_ylabel('%')
	ax.coastlines()
	
	#Create title for saved plot
	seasname = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec',\
	'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec',\
	'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
	title = intitle+'_'+str(winlen)+'mth_'+seasname[season-1]+seasname[season+winlen-2]+'_'+str(pc)+'th%ile'
	plt.title(title, fontsize=10)
	
	#Save the output
	savefile = title+'.png'
	outfile = '/Users/ailieg/Data/drought_model_eval_data/analysis/'+savefile
	plt.savefig(outfile, dpi=400, format='png',bbox_inches='tight')
	plt.close()
	


#-----------------------------------------------------------------------------	
'''
PLOT_MULTI_COMPARE_INDIV

Plots multiple plot_perc_compare_indiv 

Created by Ailie Gallant 25/07/2016
'''

def plot_multi_compare_indiv():
	import numpy as np
	
	ofile = '/Users/ailieg/Data/drought_model_eval_data/data/obs/GPCP/precip.mon.mean.nc'
	
	cmiptitle = ['ACCESS1-0_historical_r1i1p1',\
				'CanESM2_historical_r1i1p1',\
				'GFDL-CM3_historical_r1i1p1',\
				'HadGEM2-CC_historical_r1i1p1',\
				'MPI-ESM-P_historical_r1i1p1',\
				'CCSM4_historical_r1i1p1',\
				'FGOALS-s2_historical_r1i1p1',\
				'GISS-E2-R_historical_r6i1p1',\
				'NorESM1-M_historical_r1i1p1',\
				'IPSL-CM5B-LR_historical_r1i1p1']
				
	amiptitle = ['FGOALS-s2_amip_r1i1p1',\
				'GFDL-CM3_amip_r1i1p1',\
				'HadGEM2-A_amip_r1i1p1',\
				'NorESM1-M_amip_r1i1p1']
				
	cmipfile = ['CMIP5/ACCESS1-0/r1i1p1/pr/pr_Amon_ACCESS1-0_historical_r1i1p1_185001-200512.nc',\
				'CMIP5/CanESM2/r1i1p1/pr/pr_Amon_CanESM2_historical_r1i1p1_185001-200512.nc',\
				'CMIP5/GFDL-CM3/r1i1p1/pr/pr_Amon_GFDL-CM3_historical_r1i1p1_186001-200512.nc',\
				'CMIP5/HadGEM2-CC/r1i1p1/pr/pr_Amon_HadGEM2-CC_historical_r1i1p1_185912-200511.nc',\
				'CMIP5/MPI-ESM-P/r1i1p1/pr/pr_Amon_MPI-ESM-P_historical_r1i1p1_185001-200512.nc',\
				'CMIP5/CCSM4/r1i1p1/pr/pr_Amon_CCSM4_historical_r1i1p1_185001-200512.nc',\
				'CMIP5/FGOALS-s2/r1i1p1/pr/pr_Amon_FGOALS-s2_historical_r1i1p1_185001-200512.nc',\
				'CMIP5/GISS-E2-R/r6i1p1/pr/pr_Amon_GISS-E2-R_historical_r6i1p1_185001-200512.nc',\
				'CMIP5/NorESM1-M/r1i1p1/pr/pr_Amon_NorESM1-M_historical_r1i1p1_185001-200512.nc',\
				'CMIP5/IPSL-CM5B-LR/r1i1p1/pr/pr_Amon_IPSL-CM5B-LR_historical_r1i1p1_185001-200512.nc']
				
	amipfile = ['AMIP/FGOALS-s2/r1i1p1/pr/pr_Amon_FGOALS-s2_amip_r1i1p1_197901-200812.nc',\
				'AMIP/GFDL-CM3/r1i1p1/pr/pr_Amon_GFDL-CM3_amip_r1i1p1_197901_200812.nc',\
				'AMIP/HadGEM2-A/r1i1p1/pr/pr_Amon_HadGEM2-A_amip_r1i1p1_197809-200811.nc',\
				'AMIP/NorESM1-M/r1i1p1/pr/pr_Amon_NorESM1-M_amip_r1i1p1_197901-200512.nc']
	
	modfile = amipfile
	intitle = amiptitle
	
	modpath = '/Users/ailieg/Data/drought_model_eval_data/data/'
	
	
	for i in range(0,len(modfile)):
		
		mfile = modpath+modfile[i]
		it = intitle[i]
		print("CALCULATING MODEL:", it)
		
		wl = 3
		print("WINDOW LENGTH:", wl)
		sm = [3,6,9,12]
		
		for j in range(0,4):
			test = perc_compare(10,wl,sm[j],ofile, mfile, it)
			test = perc_compare(5,wl,sm[j],ofile, mfile, it)
		
		wl = 6
		print("WINDOW LENGTH:", wl)
		sm = [11,5]
		for j in range(0,2):
			test = perc_compare(10,wl,sm[j],ofile, mfile, it)
			test = perc_compare(5,wl,sm[j],ofile, mfile, it)
			
		wl = 12
		print("WINDOW LENGTH:", wl)
		sm = [1]
		for j in range(0,2):
			test = perc_compare(10,wl,sm[j],ofile, mfile, it)
			test = perc_compare(5,wl,sm[j],ofile, mfile, it)
		
		wl = 24
		print("WINDOW LENGTH:", wl)
		sm = [1]
		for j in range(0,2):
			test = perc_compare(10,wl,sm[j],ofile, mfile, it)
			test = perc_compare(5,wl,sm[j],ofile, mfile, it)
	

#-----------------------------------------------------------------------------	
'''
PLOT_MULTI_COMPARE_BOOTSTRAP

Plots multiple perc_compare_bsoblen for every model 

Created by Ailie Gallant 25/07/2016
'''

def plot_multi_compare_bootstrap():
	import numpy as np
	
	ofile = '/Users/ailieg/Data/drought_model_eval_data/data/obs/GPCP/precip.mon.mean.nc'
	
	intitle = ['ACCESS1-0_historical_r1i1p1',\
				'CanESM2_historical_r1i1p1',\
				'GFDL-CM3_historical_r1i1p1',\
				'HadGEM2-CC_historical_r1i1p1',\
				'MPI-ESM-P_historical_r1i1p1',\
				'CCSM4_historical_r1i1p1',\
				'FGOALS-s2_historical_r1i1p1',\
				'GISS-E2-R_historical_r6i1p1',\
				'NorESM1-M_historical_r1i1p1',\
				'IPSL-CM5B-LR_historical_r1i1p1',\
				'FGOALS-s2_amip_r1i1p1',\
				'GFDL-CM3_amip_r1i1p1',\
				'HadGEM2-A_amip_r1i1p1',\
				'NorESM1-M_amip_r1i1p1']
				
	modfile = ['CMIP5/ACCESS1-0/r1i1p1/pr/pr_Amon_ACCESS1-0_historical_r1i1p1_185001-200512.nc',\
				'CMIP5/CanESM2/r1i1p1/pr/pr_Amon_CanESM2_historical_r1i1p1_185001-200512.nc',\
				'CMIP5/GFDL-CM3/r1i1p1/pr/pr_Amon_GFDL-CM3_historical_r1i1p1_186001-200512.nc',\
				'CMIP5/HadGEM2-CC/r1i1p1/pr/pr_Amon_HadGEM2-CC_historical_r1i1p1_185912-200511.nc',\
				'CMIP5/MPI-ESM-P/r1i1p1/pr/pr_Amon_MPI-ESM-P_historical_r1i1p1_185001-200512.nc',\
				'CMIP5/CCSM4/r1i1p1/pr/pr_Amon_CCSM4_historical_r1i1p1_185001-200512.nc',\
				'CMIP5/FGOALS-s2/r1i1p1/pr/pr_Amon_FGOALS-s2_historical_r1i1p1_185001-200512.nc',\
				'CMIP5/GISS-E2-R/r6i1p1/pr/pr_Amon_GISS-E2-R_historical_r6i1p1_185001-200512.nc',\
				'CMIP5/NorESM1-M/r1i1p1/pr/pr_Amon_NorESM1-M_historical_r1i1p1_185001-200512.nc',\
				'CMIP5/IPSL-CM5B-LR/r1i1p1/pr/pr_Amon_IPSL-CM5B-LR_historical_r1i1p1_185001-200512.nc',\
				'AMIP/FGOALS-s2/r1i1p1/pr/pr_Amon_FGOALS-s2_amip_r1i1p1_197901-200812.nc',\
				'AMIP/GFDL-CM3/r1i1p1/pr/pr_Amon_GFDL-CM3_amip_r1i1p1_197901_200812.nc',\
				'AMIP/HadGEM2-A/r1i1p1/pr/pr_Amon_HadGEM2-A_amip_r1i1p1_197809-200811.nc',\
				'AMIP/NorESM1-M/r1i1p1/pr/pr_Amon_NorESM1-M_amip_r1i1p1_197901-200512.nc']
	
	modpath = '/Users/ailieg/Data/drought_model_eval_data/data/'
	
	
	for i in range(0,len(modfile)):
		
		mfile = modpath+modfile[i]
		it = intitle[i]+'_bootstrap'
		print("CALCULATING MODEL:", it)
		
		wl = 3
		print("NOW COMPUTING WINDOW LENGTH:", wl)
		sm = [3,6,9,12]
		
		for j in range(0,4):
			test = perc_compare_bsoblen(10,wl,sm[j],ofile, mfile, it)
			test = perc_compare_bsoblen(5,wl,sm[j],ofile, mfile, it)
		
		wl = 6
		print("NOW COMPUTING WINDOW LENGTH:", wl)
		sm = [11,5]
		for j in range(0,2):
			test = perc_compare_bsoblen(10,wl,sm[j],ofile, mfile, it)
			test = perc_compare_bsoblen(5,wl,sm[j],ofile, mfile, it)
			
		wl = 12
		print("NOW COMPUTING WINDOW LENGTH:", wl)
		sm = [1,6]
		for j in range(0,2):
			test = perc_compare_bsoblen(10,wl,sm[j],ofile, mfile, it)
			test = perc_compare_bsoblen(5,wl,sm[j],ofile, mfile, it)
		
		wl = 24
		print("NOW COMPUTING WINDOW LENGTH:", wl)
		sm = [1,6]
		for j in range(0,2):
			test = perc_compare_bsoblen(10,wl,sm[j],ofile, mfile, it)
			test = perc_compare_bsoblen(5,wl,sm[j],ofile, mfile, it)
	

#----------------------------------------------------------------------------------	
'''
BOOTSTRAP_PERC_OBS

Uses bootstrapping to randomly sample a data series 1000 times.

INPUT:
  data -  a set of observations in format of [time, lat, lon]
  n - the number of random samples required


Created by Ailie Gallant 26/07/2016
'''
def bootstrap_perc_obs(data, n, pc):
	import numpy as np
	import random
	
	#Compute percentile from the raw data
	opc = np.percentile(data, pc, axis=0)
	opc[opc < 0.0] = 0.0 #ensure percentile is not below zero, if below zero then set to 0
	
	#Compute percentile from bootstrapped data
	#Define a 3D index array
	idxtime = np.arange(data.shape[0])
	
	#Randomly select data across time index and reshape
	bs = np.random.choice(idxtime, size=n*1000, replace=True)
	
	#Reshape data, extract the bootstrapped time points and calculate percentile
	randdata = data[bs,:,:]
	randdata = randdata.reshape(1000, n, data.shape[1], data.shape[2])
	bspc = np.percentile(randdata, pc, axis=1)
	bspc[bspc < 0.0] = 0.0 #ensure percentile is not below zero, if below zero then set to 0
	
	return opc, bspc

#---------------------------------------------------------------------------------------	
	
'''
MMM_PERC_COMPARE

This script compares the percentile threshold value from observations with model data for
each model and then computes the multi-model mean of these differences. 


INPUT:
  ofile - the netcdf file of the observed data to be compared*
  pc - the percentile threshold to compare
  winlen - the length of the period over which to average (e.g. 3 months, winlen = 3)
  season - the starting month of the season of length winlen to examine, e.g. if season = 1
           and winlen = 12 this defines a january-december year.
  thresh - the number of models that are consistent in sign to show for hatching (e.g. 8 = 8 models)
  
*Inputs should be in the form of [time, lat, lon]

Created by Ailie Gallant 27/07/2016
'''
def mmm_perc_compare(pc, winlen, season, obfile, thresh):
	import numpy as np
	from grid_tools import interp_togrid
	
	#Extract model and obs data
	odata, olat, olon = rcs_gpcp(winlen)
	
	#Begin at correct time 
	rlen = -1*(season - 1)
	odata = np.roll(odata, rlen, axis=0)
	odata = odata[0::12, :, :]

	#Calculate thresholds
	opc = np.percentile(odata, pc, axis=0)
	opc[opc < 0.0] = 0.0 #ensure percentile is not below zero, if below zero then set to 0
	
	cmipfile = ['CMIP5/ACCESS1-0/r1i1p1/pr/pr_Amon_ACCESS1-0_historical_r1i1p1_185001-200512.nc',\
				'CMIP5/CanESM2/r1i1p1/pr/pr_Amon_CanESM2_historical_r1i1p1_185001-200512.nc',\
				'CMIP5/GFDL-CM3/r1i1p1/pr/pr_Amon_GFDL-CM3_historical_r1i1p1_186001-200512.nc',\
				'CMIP5/HadGEM2-CC/r1i1p1/pr/pr_Amon_HadGEM2-CC_historical_r1i1p1_185912-200511.nc',\
				'CMIP5/MPI-ESM-P/r1i1p1/pr/pr_Amon_MPI-ESM-P_historical_r1i1p1_185001-200512.nc',\
				'CMIP5/CCSM4/r1i1p1/pr/pr_Amon_CCSM4_historical_r1i1p1_185001-200512.nc',\
				'CMIP5/FGOALS-s2/r1i1p1/pr/pr_Amon_FGOALS-s2_historical_r1i1p1_185001-200512.nc',\
				'CMIP5/GISS-E2-R/r6i1p1/pr/pr_Amon_GISS-E2-R_historical_r6i1p1_185001-200512.nc',\
				'CMIP5/NorESM1-M/r1i1p1/pr/pr_Amon_NorESM1-M_historical_r1i1p1_185001-200512.nc',\
				'CMIP5/IPSL-CM5B-LR/r1i1p1/pr/pr_Amon_IPSL-CM5B-LR_historical_r1i1p1_185001-200512.nc']
				
	amipfile = ['AMIP/FGOALS-s2/r1i1p1/pr/pr_Amon_FGOALS-s2_amip_r1i1p1_197901-200812.nc',\
				'AMIP/GFDL-CM3/r1i1p1/pr/pr_Amon_GFDL-CM3_amip_r1i1p1_197901_200812.nc',\
				'AMIP/HadGEM2-A/r1i1p1/pr/pr_Amon_HadGEM2-A_amip_r1i1p1_197809-200811.nc',\
				'AMIP/NorESM1-M/r1i1p1/pr/pr_Amon_NorESM1-M_amip_r1i1p1_197901-200512.nc']
	
	modpath = '/Users/ailieg/Data/drought_model_eval_data/data/'
	modfile = cmipfile
	
	d = np.zeros((len(modfile), olat.size, olon.size))
	
	for i in range(0,len(modfile)):
		
		#Extract model data
		mdata, mlat, mlon = rcs_gpcp(winlen, modfile)
	
		#Shift data to start of required season and trim
		mdata = np.roll(mdata, rlen, axis=0)
		mdata = mdata[0::12, :, :]
	
		#Calculate thresholds
		mpc = np.percentile(mdata, pc, axis=0)
		mpc[mpc < 0.0] = 0.0 #ensure percentile is not below zero, if below zero then set to 0
	
		#Interpolate to the same grid spacing (coarsest) and compute difference as a percentage
		ob, mod, lon, lat = interp_togrid(opc, olat, olon, mpc, mlat, mlon)
		d[i,:,:] = ((mod - ob)/ob)*100
		print("Model "+str(i+1)+" of "+str(len(modfile))+" completed.")
		
	#Compute multi-model mean
	mmm = np.median(d, axis=0)
	
	#Compute the agreement in the sign of the models
	sign = d
	sign[d < 0.0] = -1
	sign[d > 0.0] = 1
	
	sign = abs(np.sum(sign, axis=0))
	
	plot = plot_perc_compare_mmm (mmm, sign, lat, lon, pc, winlen, season, thresh)
	
	return mmm, sign, lat, lon


#---------------------------------------------------------------------------------
'''
PLOT_PERC_COMPARE_MMM

This script plots the comparison of the percentile thresholds for the multi model mean.
Hatching shows where N out of M models agree in sign

INPUT:
  data - a gridded data set of the observations to be compared* 
  sign - an integer array containing the number of models where the sign is in agreement
  lat - an array of latitudes associated with the data
  lon - an array of longitudes associated with the data
  pc - the percentile threshold being investigated
  season - an integer of the starting month
  winlen - the window length for averaging
  thresh - the number of models in SIGN to consider significant and to show for hatching
  intitle - the title of the saved file 
  
*Inputs should be in the form of [time, lat, lon]

Created by Ailie Gallant 26/07/2016
'''

def plot_perc_compare_mmm (d, sign, lat, lon, pc, winlen, season, thresh):
	from matplotlib import pyplot as plt
	import cartopy.crs as ccrs
	from matplotlib.colors import BoundaryNorm
	from matplotlib.ticker import MaxNLocator
	import numpy as np
	from grid_tools import fold_grid
	
	#Set levels and colormap
	levels = np.arange(-200,220,20)
	cmap = plt.get_cmap('Spectral')
	norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
	
	#Make the grid circular
	d, lat, lon = fold_grid(d, lat, lon)
	lon = lon[0:-1]
	sign, lat, lon = fold_grid(sign, lat, lon)
	
	sign[sign < thresh] = -1.0
	sign[sign >= thresh] = 1.0
	
	#Set axes and plot
	ax = plt.axes(projection=ccrs.PlateCarree())
	p=plt.pcolormesh(lon, lat, d, cmap=cmap, norm=norm)
	cont = plt.contourf(lon, lat, sign, levels=[-1,0,1], transform = ccrs.PlateCarree(),colors='none', hatches=[" ","."])

	#Add a colorbar
	cbar = plt.colorbar(p, extend='both')
	cbar.ax.set_ylabel('%')
	ax.coastlines()
	
	#Create title for saved plot
	seasname = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec',\
	'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec',\
	'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
	title = 'MMM_CMIP_'+str(winlen)+'mth_'+seasname[season-1]+seasname[season+winlen-2]+'_'+str(pc)+'th%ile_'+str(thresh)+'modhatch'
	plt.title(title, fontsize=10)
	
	#Save the output
	savefile = title+'.png'
	outfile = '/Users/ailieg/Data/drought_model_eval_data/analysis/'+savefile
	plt.savefig(outfile, dpi=400, format='png',bbox_inches='tight')
	plt.close()
	
#-----------------------------------------------------------------------------	
'''
PLOT_MMM_BOOTSTRAP

Plots the bootstrapped individual files as well as the multi model mean

Created by Ailie Gallant 25/07/2016
'''

def plot_mmm_bootstrap(pc, wl, s):
	import numpy as np
	from netcdf_tools import ncextractall
	from matplotlib import pyplot as plt
	import cartopy.crs as ccrs
	from matplotlib.colors import BoundaryNorm
	from matplotlib.ticker import MaxNLocator
	from grid_tools import fold_grid
	
	ofile = '/Users/ailieg/Data/drought_model_eval_data/data/obs/GPCP/precip.mon.mean.nc'
	
	cmiptitle = ['ACCESS1-0_historical_r1i1p1',\
				'CanESM2_historical_r1i1p1',\
				'GFDL-CM3_historical_r1i1p1',\
				'HadGEM2-CC_historical_r1i1p1',\
				'MPI-ESM-P_historical_r1i1p1',\
				'CCSM4_historical_r1i1p1',\
				'FGOALS-s2_historical_r1i1p1',\
				'GISS-E2-R_historical_r6i1p1',\
				'NorESM1-M_historical_r1i1p1',\
				'IPSL-CM5B-LR_historical_r1i1p1']
				
	amiptitle = ['FGOALS-s2_amip_r1i1p1',\
				'GFDL-CM3_amip_r1i1p1',\
				'HadGEM2-A_amip_r1i1p1',\
				'NorESM1-M_amip_r1i1p1']
				
	cmipfile = ['CMIP5/ACCESS1-0/r1i1p1/pr/pr_Amon_ACCESS1-0_historical_r1i1p1_185001-200512.nc',\
				'CMIP5/CanESM2/r1i1p1/pr/pr_Amon_CanESM2_historical_r1i1p1_185001-200512.nc',\
				'CMIP5/GFDL-CM3/r1i1p1/pr/pr_Amon_GFDL-CM3_historical_r1i1p1_186001-200512.nc',\
				'CMIP5/HadGEM2-CC/r1i1p1/pr/pr_Amon_HadGEM2-CC_historical_r1i1p1_185912-200511.nc',\
				'CMIP5/MPI-ESM-P/r1i1p1/pr/pr_Amon_MPI-ESM-P_historical_r1i1p1_185001-200512.nc',\
				'CMIP5/CCSM4/r1i1p1/pr/pr_Amon_CCSM4_historical_r1i1p1_185001-200512.nc',\
				'CMIP5/FGOALS-s2/r1i1p1/pr/pr_Amon_FGOALS-s2_historical_r1i1p1_185001-200512.nc',\
				'CMIP5/GISS-E2-R/r6i1p1/pr/pr_Amon_GISS-E2-R_historical_r6i1p1_185001-200512.nc',\
				'CMIP5/NorESM1-M/r1i1p1/pr/pr_Amon_NorESM1-M_historical_r1i1p1_185001-200512.nc',\
				'CMIP5/IPSL-CM5B-LR/r1i1p1/pr/pr_Amon_IPSL-CM5B-LR_historical_r1i1p1_185001-200512.nc']
				
	amipfile = ['AMIP/FGOALS-s2/r1i1p1/pr/pr_Amon_FGOALS-s2_amip_r1i1p1_197901-200812.nc',\
				'AMIP/GFDL-CM3/r1i1p1/pr/pr_Amon_GFDL-CM3_amip_r1i1p1_197901_200812.nc',\
				'AMIP/HadGEM2-A/r1i1p1/pr/pr_Amon_HadGEM2-A_amip_r1i1p1_197809-200811.nc',\
				'AMIP/NorESM1-M/r1i1p1/pr/pr_Amon_NorESM1-M_amip_r1i1p1_197901-200512.nc']
	
	modfile = cmipfile
	intitle = cmiptitle
	
	modpath = '/Users/ailieg/Data/drought_model_eval_data/data/'
	
	
	obsnc = ncextractall(ofile)
	lon = obsnc['lon']
	lat = obsnc['lat']
	
	sigmod = np.zeros((len(modfile), lat.size, lon.size))
	
	for i in range(0,len(modfile)):
		
		mfile = modpath+modfile[i]
		it = intitle[i]

		d, sig, lat, lon = perc_compare_bsoblen(pc, wl, s, ofile, mfile, it)
		
		sig[sig < 0.0] = 0.0
		sigmod[i,:,:] = sig
		
	sumsig = np.sum(sigmod, axis=0)/len(modfile)
	
	#Set levels and colormap
	levels = np.arange(0,100,10)
	cmap = plt.get_cmap('Spectral')
	norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
	
	#Make the grid circular
	sumsig, lat, lon = fold_grid(sumsig, lat, lon)
	
	#Set axes and plot
	ax = plt.axes(projection=ccrs.PlateCarree())
	p=plt.pcolormesh(lon, lat, d, cmap=cmap, norm=norm)
	
	#Add a colorbar
	cbar = plt.colorbar(p, extend='both')
	cbar.ax.set_ylabel('%')
	ax.coastlines()
	
	#Create title for saved plot
	seasname = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec',\
	'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec',\
	'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
	title = 'PERC_MODELS_CMIP_'+str(wl)+'mth_'+seasname[s-1]+seasname[s+wl-2]+'_'+str(pc)+'th%ile'
	plt.title(title, fontsize=10)
	
	#Save the output
	savefile = title+'.png'
	outfile = '/Users/ailieg/Data/drought_model_eval_data/analysis/'+savefile
	plt.savefig(outfile, dpi=400, format='png',bbox_inches='tight')
	plt.close()	
		
		
#----------------------------------------------------------------------------------------
#SPI_QQ
'''
spi_qq

Purpose:
---------
To do a quantile/quantile (QQ) comparison of the lower tails of precipitation distributions
below particular SPI thresholds. A QQ comparison is made between the observed data and
that from a modelled Pareto and GEV distribution. The RMSE for each method is computed
to determine the model that best fits the data.

Input parameter:
-----------------
SPI number (3, 6, 12 or 24)
SEASON to be examined (1-12)

History:
---------
20160728 - Created by Ailie Gallant, Monash University

'''
def spi_qq(spin, season):
	import numpy as np
	from grid_tools import fill_missing
	from matplotlib import pyplot as plt
	

	#Read in observed data
	odata, olat, olon = rcs_gpcp(winlen)
	

	return 
	
