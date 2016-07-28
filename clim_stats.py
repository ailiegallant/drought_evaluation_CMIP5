#THIS MODULE CONTAINS TOOLS AND CODE FOR STATISTICAL ANALYSIS
#EACH SUBROUTINE IS SEPARATED BY A LINE --------
#-----------------------------------------------------------------------------

# PROGRAM TESTS THE GENERATION OF THE SYNTHETIC SERIES AND TESTS AGAINST THE 
# REAL DATA

'''
	Synthetic series are created using a model generated from a gamma distribution
	modelled on monthly time series. A separate gamma distribution is modelled for
	each month to ensure that any seasonality in the time series is captured.
	
	INPUT
	file - file name and path of the netcdf file containing the climate data from which
	       the synthetic time series are created
	vname - the variable name to read in the netcdf data file named 'file'
	
	OUTPUT
	A synthetic time series modelled on the input data
	
	HISTORY
	20160505 - created by Ailie Gallant
	
'''


def synth_climseries(file, vname, nsyn) :

	import numpy as np
	import scipy.stats as stats
	from netcdf_tools import ncextractall
	from netCDF4 import num2date

	#Import and extract the netcdf file, returns a dictionary
	datdict = ncextractall(file)
	
	#Read each variable
	data = datdict[vname]
	lon = datdict['lon']
	lat = datdict['lat']
	time = datdict['time']
	
	#convert time units to actual date
	time_u = datdict['time_units']
	if 'time_calendar' in datdict.keys(): 
		cal = datdict['time_calendar']
		time = num2date(time,units = time_u, calendar=cal)
	else: time = num2date(time,units = time_u)
	
	#convert missing numbers to NaN
	miss = datdict[vname+"_missing_value"] 
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
	
	#Create an empty array that will contain the synthetic series
	synthseries = np.zeros((lat.size, lon.size, time.size, nsyn))
	synthtest = np.zeros((lat.size, lon.size, time.size, nsyn))
	des = np.zeros((time.size, lat.size, lon.size))
	ac = np.zeros((lat.size, lon.size))
	
	for i in range(0,len(lat)):
		print ("Lat ", i, "of ", len(lat))
		for j in range(0,len(lon)):

			series= data[:,i,j]
			series = series.flatten()
			
			
			
			#The data must be stratified into months so that the SPI may be computed
			#for each month so data is deseasonalised. 
			#Reshape the array to stratify into months.
			lens = len(series)
			lens = lens/12
			lens = int(lens)
			series = series.reshape([lens,12])
			
			#compute auto-correlation of the deseasonalised time series (to remove auto-
			#correlation associated with seasonality, which will automatically be removed
			#when gamma distributions are computed seasonally)
			sm = series.mean(axis=0)
			seriessm = np.tile(sm,36)
			deseas = (series.flatten()) - seriessm
			des[:,i,j] = deseas
			lag = np.roll(deseas, -1)
			ac[i,j] = stats.pearsonr(deseas, lag)[0]
			
			#Compute NSYN number of synthetic series
			
			for m in range(0,nsyn-1):
				
				tmpsyn = np.zeros((lens,12))
				 
				#Compute the parameters for the gamma distribution, one month at a time
				for k in range(0,12):    
					
					tmp = series[:,k]
					tmpsave = tmp
					
					#remove any NaNs (i.e. missing numbers) from data so only real numbers exist
					tmp = tmp[~np.isnan(tmp)]
					
					
					if len(tmp) > 10:
						#compute the number of zeros
						numzeros = (tmp[np.where(tmp == 0.0)]).size
			
						#compute the probability of zeros based on the sample series
						q = numzeros/tmp.size
						#compute the probability of non-zeros based on the sample series
						p = 1.0 - q
						
						
						#compute the shape, scale and location parameters based on non-zero data only
						nonzerotmp = tmp[np.where(tmp > 0.0)]
						numnonzero = nonzerotmp.size
						A = np.log(np.mean(nonzerotmp)) - (np.sum(np.log(nonzerotmp))/len(nonzerotmp))
						shp = (1.0/(4*A)) * (1 + ((1 + ((4*A)/3) )**0.5))
						scl = np.mean(nonzerotmp)/shp 
		
						#test bit-------
						
						
						
						
						#--------------
		
		
		
						#Compute synthetic distribution of non-zero values 
						synthgam = stats.gamma.rvs(shp, scale=scl, size=numnonzero)
		
						if q > 0.0: 
							zeroarr = np.zeros(numzeros)
							synthgam = np.concatenate((zeroarr, synthgam))
							np.random.shuffle(synthgam)
						
						tmpsyn[:,k] = synthgam
						tmps = tmpsyn.flatten()
						
					else: tmps = np.zeros(len(tmpsave)) + miss
					
				tmps = tmps.flatten()
				synthseries[i,j,:,m] = tmps + ((ac[i,j])*tmps)
				synthtest[i,j,:,m] = tmps
		

	
	return synthseries, synthtest, data, des, lon, lat, time, time_u, miss, ac


#---------------------------------------------------------------
'''
	Synthetic series are created using a model generated from a gamma distribution
	modelled on monthly time series. A separate gamma distribution is modelled for
	each month to ensure that any seasonality in the time series is captured.
	
	INPUT
	file - file name and path of the netcdf file containing the climate data from which
	       the synthetic time series are created
	vname - the variable name to read in the netcdf data file named 'file'
	
	OUTPUT
	A synthetic time series modelled on the input data
	
	HISTORY
	20160505 - created by Ailie Gallant
	
'''


def synth_test(file, vname, nsyn) :

	import numpy as np
	import scipy.stats as stats
	from netcdf_tools import ncextractall
	from netCDF4 import num2date

	#Import and extract the netcdf file, returns a dictionary
	datdict = ncextractall(file)
	
	#Read each variable
	data = datdict[vname]
	lon = datdict['lon']
	lat = datdict['lat']
	time = datdict['time']
	
	#convert time units to actual date
	time_u = datdict['time_units']
	if 'time_calendar' in datdict.keys(): 
		cal = datdict['time_calendar']
		time = num2date(time,units = time_u, calendar=cal)
	else: time = num2date(time,units = time_u)
	
	#convert missing numbers to NaN
	miss = datdict[vname+"_missing_value"] 
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
	
	#Create an empty array that will contain the synthetic series
	synthseries = np.zeros((lat.size, lon.size, time.size, nsyn))
	des = np.zeros((time.size, lat.size, lon.size))
	
	for i in range(1,len(lat)):
		print ("Lat ", i, "of ", len(lat))
		for j in range(1,len(lon)):

			series= data[:,i,j]
			series = series.flatten()
			
			
			#The data must be stratified into months so that the SPI may be computed
			#for each month so data is deseasonalised. 
			#Reshape the array to stratify into months.
			lens = len(series)
			lens = lens/12
			lens = int(lens)
			series = series.reshape([lens,12])
			
			#compute auto-correlation of the deseasonalised time series (to remove auto-
			#correlation associated with seasonality, which will automatically be removed
			#when gamma distributions are computed seasonally)
			sm = series.mean(axis=0)
			seriessm = np.tile(sm,lens)
			deseas = (series.flatten()) - seriessm
			des[:,i,j] = deseas
			ac = np.zeros(6)
			for l in range(1,7):ac[l-1]=stats.pearsonr(deseas, np.roll(deseas,(-1*l)))[0]
			#reverse the autocorrelation function for use later
			ac = ac[::-1]
				
			tmpsyn = np.zeros((lens,12))
			shp = np.zeros(12)
			scl = np.zeros(12)
				
			#make an array of 1s for later - used to store zeros if necessary
			zeros = np.zeros((lens,12))+1.0
				 
			#Compute the parameters for the gamma distribution, one month at a time
			for k in range(0,12):    
					
				tmp = series[:,k]
				tmpsave = tmp
					
				#remove any NaNs (i.e. missing numbers) from data so only real numbers exist
				tmp = tmp[~np.isnan(tmp)]
					
				if len(tmp) > 10:
					#compute the number of zeros
					numzeros = (tmp[np.where(tmp == 0.0)]).size
			
					#compute the probability of zeros based on the sample series
					q = numzeros/tmp.size
					#compute the probability of non-zeros based on the sample series
					p = 1.0 - q
					
					#compute the shape, scale and location parameters based on non-zero data only
					nonzerotmp = tmp[np.where(tmp > 0.0)]
					numnonzero = nonzerotmp.size
					A = np.log(np.mean(nonzerotmp)) - (np.sum(np.log(nonzerotmp))/len(nonzerotmp))
					shp[k] = (1.0/(4*A)) * (1 + ((1 + ((4*A)/3) )**0.5))
					scl[k] = np.mean(nonzerotmp)/shp[k]
					
					tmpz = np.zeros(lens)+1.0
					tmpz[0:numzeros]=0.0
					np.random.shuffle(tmpz)
					
					zeros[:,k] = tmpz

						
				else: tmps = np.zeros(len(tmpsave)) + miss
			

			#Reshape shape and scale
			alpha = np.tile(shp,lens)
			beta = np.tile(scl, lens)
			zeros = zeros.flatten()
			
			#Compute NSYN number of synthetic (surrogate) series
			for m in range(0,nsyn):
				#noise = ((alpha[0]*beta[0]) + (alpha[0]*(beta[0]**2))*(stats.norm.rvs()))
				noise = stats.gamma.rvs(alpha[0], scale = beta[0], size=6) 
				surgam = np.zeros((time.size))
				surgam[0:6] = noise
				
				#Compute surrogate data using form x(t) = a0 + a1(xt-1) + sigma*e1
				#for w in range(7, time.size):
				#	if zeros[w] < 1.0: 
				#		surgam[w] = 0.0
				#		#noise = ((alpha[w]*beta[w]) + (alpha[w]*(beta[w]**2))*(stats.norm.rvs()))
				#		noise = stats.gamma.rvs(alpha[w], scale = beta[w]) 
				#	else: 
				#		noise = stats.gamma.rvs(alpha[w], scale = beta[w])  
				#		surgam[w] = noise + (np.sum(ac*surgam[w-7:w-1])) 
				#if surgam[100] > 1000: print(surgam[100])
				surgam = x
				
				synthseries[i,j,:,m] = surgam
			#synthtest[i,j,:,m] = tmps
		

	
	return synthseries, data, alpha,beta,zeros, ac, surgam


#----------------------------------------------------------------------------------	

'''
Creates a monthly-stratified synthetic series from a gamma distribution

INPUT:
Data - a monthly time series of length, n, beginning in January and ending in December
'''
def synth_gamma_seasonal(data) :

	import numpy as np
	import scipy.stats as stats

	series= data

	#The data must be stratified into months so that the synthetic distribution may be computed
	#for each month so data is deseasonalised. 
	#Reshape the array to stratify into months.
	lens = series.size
	lens = lens/12
	lens = int(lens)
	series = series.reshape([lens,12])
		
	syn = np.zeros((lens,12))
		 
	#Compute the parameters for the gamma distribution, one month at a time
	for k in range(0,12):    
		
		tmp = series[:,k]
		
		#remove any NaNs (i.e. missing numbers) from data so only real numbers exist
		tmp = tmp[~np.isnan(tmp)]
		
		if tmp.size > 10:
			#compute the number of zeros
			numzeros = (tmp[np.where(tmp == 0.0)]).size
	
			#compute the probability of zeros based on the sample series
			q = numzeros/tmp.size
			#compute the probability of non-zeros based on the sample series
			p = 1.0 - q
			
			#compute the shape and scale parameters based on non-zero data only
			nonzerotmp = tmp[np.where(tmp > 0.0)]
			numnonzero = nonzerotmp.size
			A = np.log(np.mean(nonzerotmp)) - (np.sum(np.log(nonzerotmp))/len(nonzerotmp))
			shp = (1.0/(4*A)) * (1 + ((1 + ((4*A)/3) )**0.5))
			scl = np.mean(nonzerotmp)/shp 
			
			#Compute synthetic distribution of non-zero values 
			synthgam = stats.gamma.rvs(shp, scale=scl, size=numnonzero)

			if q > 0.0: 
				zeroarr = np.zeros(numzeros)
				synthgam = np.concatenate((zeroarr, synthgam))
				np.random.shuffle(synthgam)
			
			syn[:,k] = synthgam
			
		else: syn = np.zeros(lens) + np.nan

	syn = syn.flatten()

	return syn

#----------------------------------------------------------------------------------	

'''
Creates a synthetic series from a gamma distribution

INPUT:
Data - a monthly time series of length, n, beginning in January and ending in December
'''
def synth_gamma(data) :

	import numpy as np
	import scipy.stats as stats

	series= data

	#The data must be stratified into months so that the synthetic distribution may be computed
	#for each month so data is deseasonalised. 
	#Reshape the array to stratify into months.
	lens = series.size
	lens = lens/12
	lens = int(lens)
	series = series.reshape([lens,12])
		
	syn = np.zeros((lens,12))
		 
	#Compute the parameters for the gamma distribution, one month at a time
	for k in range(0,12):    
		
		tmp = series[:,k]
		
		#remove any NaNs (i.e. missing numbers) from data so only real numbers exist
		tmp = tmp[~np.isnan(tmp)]
		
		if tmp.size > 10:
			#compute the number of zeros
			numzeros = (tmp[np.where(tmp == 0.0)]).size
	
			#compute the probability of zeros based on the sample series
			q = numzeros/tmp.size
			#compute the probability of non-zeros based on the sample series
			p = 1.0 - q
			
			#compute the shape and scale parameters based on non-zero data only
			nonzerotmp = tmp[np.where(tmp > 0.0)]
			numnonzero = nonzerotmp.size
			A = np.log(np.mean(nonzerotmp)) - (np.sum(np.log(nonzerotmp))/len(nonzerotmp))
			shp = (1.0/(4*A)) * (1 + ((1 + ((4*A)/3) )**0.5))
			scl = np.mean(nonzerotmp)/shp 
			
			#Compute synthetic distribution of non-zero values 
			synthgam = stats.gamma.rvs(shp, scale=scl, size=numnonzero)

			if q > 0.0: 
				zeroarr = np.zeros(numzeros)
				synthgam = np.concatenate((zeroarr, synthgam))
				np.random.shuffle(synthgam)
			
			syn[:,k] = synthgam
			
		else: syn = np.zeros(lens) + np.nan

	syn = syn.flatten()

	return syn
	
	
#----------------------------------------------------------------------------------	

'''
Compute correlations between seasonal rainfall and annual rainfall 

'''

def gpcp_corr(sone, stwo):
	from netcdf_tools import ncextractall
	import numpy as np
	from scipy import stats
	from netCDF4 import num2date
	
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
	
	corr = np.zeros((nlat,nlon)) +  obsmiss
	
	
	for i in range(0,nlon):
		for j in range(0,nlat):
		
			series = obsp[:,j,i]
			lens = series.size
			lens = lens/12
			lens = int(lens)
			series = series.reshape([lens,12])
			
			if sone[0] > sone[sone.size-1]:
				sidxone = np.arange(0,sone.size)
				series = np.roll(series,-1*(sone[0]-1))
				seriesone = np.mean(series[:,sidxone],axis=1)
			else:
				sidxone = sone - 1
				seriesone = np.mean(series[:,sidxone],axis=1)
			
			if stwo[0] > stwo[stwo.size-1]:
				sidxone = np.arange(0,sone.size)
				series = np.roll(series,-1*(sone[0]-1))
				seriestwo = np.mean(series[:,sidxtwo],axis=1)
			else:
				sidxtwo = stwo - 1
				seriestwo = np.mean(series[:,sidxtwo],axis=1)
			
			corr[j,i] = stats.pearsonr(seriesone, seriestwo)[0]
	
	return corr, lat, lon
	
	
	
	
	
	
	