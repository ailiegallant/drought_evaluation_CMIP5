#PRECIP_MODEL_TESTS
#--------------------------------------------------------------------------------------
#Contains a series of code which with compare precipitation data to various statistical
#models, or, do synthetic modelling of precipitation.
#--------------------------------------------------------------------------------------


'''
qq_test_extremes(a, **kwargs)

Performs a quantile-quantile comparison of an input array, a and either a
Generalised Extreme Value or Generalised Pareto distribution. Returns the
RMSE of the residuals of each.

INPUT
-----
a - the data to test in terms of quantiles

OUTPUT
-----
a two-element array containing the RMSE from the GEV and GPD
'''

def qq_test_extremes(a):
	import numpy as np
	from scipy import stats
	
	#Fit parameters
	gevp = stats.genextreme.fit(a)
	parp = stats.genpareto.fit(a)
	
	#Sort the data
	a = np.sort(a)
	
	#Create data for the 
	
	
	
	
return


