#THIS MODULE HAS SCRIPTS WHICH DO CONVERSIONS BETWEEN VARIOUS UNITS
#
#EACH SUBROUTINE IS SEPARATED BY A LINE --------
#-----------------------------------------------------------------------------
'''
MMD_MMM

This script converts monthly rainfall data from a 3D rainfall grid in measurements of 
mm/day into measurements of mm/month. It assumes traditional days in month and 28 days 
for February.

INPUT:
  data - a gridded data set of the observations to be compared*
  
*Inputs should be in the form of [time, lat, lon]

Created by Ailie Gallant 21/07/2016
'''

def mmd_mmm(data):
	import numpy as np
	
	done = data.shape[1]
	dtwo = data.shape[2]
	tdim = data.shape[0]
	
	mon = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
	mon = np.tile(mon, tdim/12)
	mon = np.repeat(mon[:, np.newaxis], done, axis=1)
	mon = np.repeat(mon[:,:,np.newaxis],dtwo,axis=2)
	
	return mon*data


