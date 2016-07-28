#THIS MODULE CONTAINS TOOLS AND CODE FOR STATISTICAL METRICS
#EACH SUBROUTINE IS SEPARATED BY A LINE --------
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------

'''
perk_skillscore
Computes the Perkins et al. (2007) skill score via comparison of histograms
S = sum(0, n)min(Zx, Zy)
where n is the number of bins, Zx is the relative frequency in bin n from the X array and
Zy is the relative frequency in bin n from the Y array. 

INPUT
------
X - a data array for comparison
Y - a data array for comparison
BINS - an array containing the values of the bounds of the bins to compare. The number
of bins should be 10% of the size of the length of the X and Y array if X and Y are the same
(similar) lengths.

OUTPUT
--------
S, the Perkins Skill Score 

HISTORY
-------
20160524 - Created by Ailie Gallant
'''

def perk_skillscore(x, y, bins): 
	import numpy as np
	
	#Create histograms
	hx = np.histogram(x,bins=bins)
	hx = hx[0]/np.sum(hx[0])
	
	hy = np.histogram(y,bins=bins)
	hy = hy[0]/np.sum(hy[0])
	
	#Calculate skill scores
	S = np.sum(np.min(np.vstack((hx,hy)), axis=0))
	
	return S