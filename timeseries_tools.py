#THIS MODULE CONTAINS TOOLS AND CODE FOR TIME SERIES ANALYSIS
#EACH SUBROUTINE IS SEPARATED BY A LINE --------
#-----------------------------------------------------------------------------

# Function computes a moving average from a given input array
# and returns an array of moving averages
# 
#  a - the input array (1D)
#  n - the number of elements of the input array over which to average


def moving_average(a, n) :
    import numpy as np

#compute the sum over n points and divide by n. Result is "ma", the moving average.
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    ma = ret[n - 1:] / n
    
    #make empty array for nan, first n-1 values 
    tmp = np.zeros(n-1)+np.nan
    
    #input empty arrays at the start and end of the ma array
    ma = np.concatenate([tmp,ma])
    
    return ma
    
    
#-----------------------------------------------------------------------------
    '''
    statsfromgrid returns four scalars
    
    computes the first 3 moments (mean, stdev & skewness) and autocorrelation 
    (lag 1) of a gridded dataset where each grid contains a time series.

    Inputs
    ----------
    file : filename with full path of netcdf data set to read
    vname: a string containing the variable name to read

    Returns
    -------
    mean, stdev, skew, acorr: four 2D arrays containing moments from each grid point
        
    Code development
    -------
    20160429 - script developed by Ailie Gallant
    '''


def statsfromgrid(file, vname) :
	import numpy as np
	import scipy.stats as stats
	from netcdf_tools import ncextractall
	import collections

	#Import and extract the netcdf file, returns a dictionary
	datdict = ncextractall(file)
	
	#Read each variable
	data = datdict[vname]
	lon = datdict['lon']
	lat = datdict['lat']
	time = datdict['time']
	
	#Create empty arrays for output
	mean = np.zeros((len(lat), len(lon)))
	stdev = np.zeros((len(lat), len(lon)))
	skew = np.zeros((len(lat), len(lon)))
	acorr = np.zeros((len(lat), len(lon)))
	numzeros = np.zeros((len(lat), len(lon)))
	
	#Loop over lon and lat monthly data and record the statistics
	for i in range(0,len(lat)-1):
	    for j in range(0,len(lon)-1):
	         tmpdata = data[:,i,j]
	         tmpdata = tmpdata.flatten()
	         mean[i,j] = np.mean(tmpdata)
	         stdev[i,j] = np.std(tmpdata) 
	         skew[i,j] = stats.skew(tmpdata)
	         numzeros[i,j] = (tmpdata[np.where(tmpdata == 0.0)]).size
	         
	         tmp = collections.deque(tmpdata)
	         tmp.rotate(-1)
	         tmp = np.asarray(tmp)
	         acorr[i,j] = stats.pearsonr(tmpdata,tmp)[0]
	         if (stats.pearsonr(tmpdata,tmp)[0]) > 0.8: print(i,j,stats.pearsonr(tmpdata,tmp)[0])


	return numzeros.flatten(), acorr.flatten()
	
#--------------------------------------------------------------------------------------
'''
Does a block bootstrap of data, using blocks of size block_size

'''
def block_bootstrap(data, block_size):
	import numpy as np
	import random
	
	#Define sample size and block size
	nlen = data.size
	numblocks = int(nlen/block_size)
	leftover = nlen - (numblocks*block_size)
	
	#If the mean block size is not a factor of the length of the time series,
	#select some random choices from the data series as residuals
	resid = np.random.choice(data,size=leftover,replace=True)
	
	#Remaining random data will be block size *
	datalong = np.zeros((data.size, block_size))
	datalong[:,0] = data
	
	#Create an array with blocked data that contains every starting point in the time series
	i = 1
	while i < block_size: 
		datalong[:,i] = np.roll(data,(-1*i))
		i = i+1
	
	#Clip the data to be the length of the block size
	datalong = datalong[:(-1*(block_size-1)),:]
	
	#Extract random blocks of data
	blk = np.arange(0,datalong.shape[0])
	blk = np.random.choice(blk, size = numblocks, replace=True)
	blk = np.sort(blk)
	choice = datalong[blk,:]
	
	#Flatten and append the residual data to the end of the time series
	choice = choice.flatten()
	choice = np.concatenate((choice,resid))
	
	if choice.size != nlen: print("input and output sample sizes not the same!!")
	
	return choice
#--------------------------------------------------------------------------------------
'''
Gives you the average length of a run of a time series given a sub-index of an array
'''

def average_run(idx):
	import numpy as np
	
	idxsum = idx - np.roll(idx,1)
	idxsum = np.where(idxsum == 1,1,0)
	
	len = idxsum.size
	count = np.zeros(len)
	
	i = 1
	run = idxsum[0]
	c = 0
	
	while i < len:
		check = idxsum[i]
		if check > 0:
			run = run + check
			count[c] = run
		else: 
			run = 0
			c = c+1
		i = i+1
	
	count = count[np.where(count > 0)]
	count = count + 1
	
	return round(np.mean(count))
		
			 
		
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
'''
BLOCK_BOOTSTRAP_3D

Calculates a block bootstrap of data, using blocks of size block_size over a 3D time
series (input in form of [time, lat, lon]). The script will make a dummy index array 
to do the calculations and will then index data appropriately.

Created by Ailie Gallant 26/07/2016
'''

def block_bootstrap(data, block_size):
	import numpy as np
	import random
	
	#Define a 3D index array
	idx = np.arange(data.size).reshape(data.shape)
	
	#Compute randomised values along the time dimension (dim [0])
	
	#Define sample size and block size
	nlen = data.shape[0]
	numblocks = int(nlen/block_size)
	leftover = nlen - (numblocks*block_size)
	
	#If the mean block size is not a factor of the length of the time series,
	#select some random choices from the data series as residuals
	resid = np.random.choice(data,size=leftover,replace=True)
	
	#Remaining random data will be block size *
	datalong = np.zeros((data.size, block_size))
	datalong[:,0] = data
	
	#Create an array with blocked data that contains every starting point in the time series
	i = 1
	while i < block_size: 
		datalong[:,i] = np.roll(data,(-1*i))
		i = i+1
	
	#Clip the data to be the length of the block size
	datalong = datalong[:(-1*(block_size-1)),:]
	
	#Extract random blocks of data that are numblocks * nlat * nlon in length so that
	#indices for all grid points are extracted
	blk = np.arange(0,datalong.shape[0])
	ntote = numblocks * data.shape[1] * data.shape[2]
	blk = np.random.choice(blk, size = ntote, replace=True)
	choice = datalong[blk,:,:,:]
	
	#Flatten and append the residual data to the end of the time series
	choice = choice.flatten()
	choice = np.concatenate((choice,resid), axis=0)
	
	if choice.size != nlen: print("input and output sample sizes not the same!!")
	
	return choice
#--------------------------------------------------------------------------------------

