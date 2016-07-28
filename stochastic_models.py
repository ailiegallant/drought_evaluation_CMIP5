#THIS MODULE CONTAINS TOOLS AND CODE FOR STOCHASTIC MODELLING
#EACH SUBROUTINE IS SEPARATED BY A LINE --------
#-----------------------------------------------------------------------------



'''
pgar returns stochastic simulations using a PGAR(1) process
as in Fernandez and Salas (1986) (page 1388-1389 for algorithm).
PGAR(1) process is a periodic autoregressive gamma process. Seasonality
is assumed in the input array (12-month). The script assumes monthly data
is provided.


Inputs
----------
series : a one-dimensional time series or array of monthly
nsurrogate: the number of surrogate series to produce

Returns
-------
Stochastic autoregressive gamma surrogate data set modelled from the input array.
    
Code development
-------
20160520 - script developed by Ailie Gallant
'''


def pgar_mth(series, nsurrogate):
	import numpy as np
	import scipy.stats as stats
	
	#Save the original array for later use
	orig = series
	
	
	#The data must be stratified into months
	#for each month so data is deseasonalised. 
	#Reshape the array to stratify into months.
	lens = len(series)
	lens = lens/12
	lens = int(lens)
	series = series.reshape([lens,12])
	
	'''
	Parameters required:
	  shp(m) - the shape parameter for the gamma distribution as a function of month (m)
	  scl(m) - the scale parameter for the gamma distribution as a function of month (m)
	  loc(m) - the location parameter for the gamma distribution as a function of month (m)
	  phi(m) - the periodic autoregressive coefficient as a function of month (m)
	  delta(m) - the periodic autoregressive exponent as a function of month (m)
	'''
	#Compute the parameters for each season, 1-month at a time.
	shp = np.zeros(12)
	scl = np.zeros(12)
	loc = np.zeros(12)
	ar = np.zeros(12)
	q = np.zeros(12)
	p = np.zeros(12)
	
	for k in range(0,12):    
			
		tmp = series[:,k]
			
		#remove any NaNs (i.e. missing numbers) from data so only real numbers exist
		tmp = tmp[~np.isnan(tmp)]
		
		#Only compute the distributions for the gamma distribution if the series is 
		#larger than 10, otherwise the sample size is too small to compute the parameters
		if tmp.size > 10:
			
			#compute the number of zeros
			numzeros = (tmp[np.where(tmp == 0.0)]).size
			q[k] = numzeros
			
			#compute the shape, scale and location parameters based on non-zero data only
			nonzerotmp = tmp[np.where(tmp > 0.0)]
			numnonzero = nonzerotmp.size
			p[k] = numnonzero
			
			A = np.log(np.mean(nonzerotmp)) - (np.sum(np.log(nonzerotmp))/len(nonzerotmp))
			
			#fit_shp, fit_loc, fit_scl=stats.gamma.fit(nonzerotmp)
			#shp[k] = fit_shp
			#scl[k] = fit_scl
			mn = np.mean(nonzerotmp)
			sd = np.std(nonzerotmp)
			skew = stats.skew(nonzerotmp)
			shp[k] = (1.0/(4*A)) * (1 + ((1 + ((4*A)/3) )**0.5))
			scl[k] = np.mean(nonzerotmp)/shp[k]
			loc[k] = mn - ((2*sd)/skew)
			ar[k] = abs(stats.pearsonr(tmp, np.roll(tmp,(1)))[0])
			
		else: print("Sample size is inadequate to compute robust parameters.")
	
	print(shp)
	print(scl)
	print(loc)	
	prevscl = np.roll(scl,(1))
	prevshp = np.roll(shp,(1))
	sclcomp = (np.roll(scl,(1)) > scl)
	
	delta = ar*((prevscl/scl)**0.5)
	delta = np.where(sclcomp == True,delta,0)
	
	phi = ar*(shp/prevshp)*(scl/prevscl)**0.5
	phi = np.where(sclcomp == False, phi, 0)
	
	#Reform the parameters to be the shape of the original series for computation
	shp = np.tile(shp, lens)
	scl = np.tile(scl, lens)
	loc = np.tile(loc, lens)
	phi = np.tile(phi, lens)
	delta = np.tile(delta, lens)
	ar = np.tile(ar, lens)
	q = np.tile(q,lens)
	p = np.tile(p,lens)
	nlen = shp.size

	
	#Empty array to save the data
	surr = np.zeros((nsurrogate, nlen))
	betaout = np.zeros((nsurrogate, nlen))
	Wout = np.zeros((nsurrogate, nlen))
	
	#Compute NSURROGATE number of surrogate (synthetic) series
	for m in range(0,nsurrogate):
		
			#Begin with a seed from a random distribution generated from the month before
		#the first month (i.e. if month 1, then use params from month 12). Here, we 
		#use month 12 (i.e. index 11).
		#Select a random value from a 1-in-lens number of years, which also includes
		#the probability of a zero. Note the X is two parameters only.
		seednonz = stats.gamma.rvs(shp[11],scale = scl[11], size=p[11])
		seedz = np.zeros(q[11])
		X = np.concatenate((seedz,seednonz))
		np.random.shuffle(X)
		X = X[0] - loc[0]
		
		#Make the first value in the stochastic series a random gamma variable
		sX = stats.gamma.rvs(shp[11], scale = scl[11], loc=loc[11])
		
		a = [shp[nlen-1], shp[0]]
		b = [scl[nlen-1], scl[0]]
		l = [loc[nlen-1], loc[0]]
		ph = [phi[nlen-1], phi[0]]
		d = [delta[nlen-1], delta[0]]
			
		W = pgar_noise(a, b, l, ph, d)
		
		surr[m,0] = (phi[0]*X)+((sX**delta[0])*W)
		
		for j in range(1, nlen):
			
			#sX is now the value at the previous timestep
			sX = surr[m,j-1]
			
			#X is a random variate from the two-param gamma distribution with params
			#taken from time step t-1.
			X = stats.gamma.rvs(shp[j-1], scale = scl[j-1]) - loc[j-1]
			
			#Remember that by going to j+1 this will give you j-1  and j only.
			a = shp[j-1:j+1]
			b = scl[j-1:j+1]
			l = loc[j-1:j+1]
			ph = phi[j-1:j+1]
			d = delta[j-1:j+1]

			W = pgar_noise(a, b, l, ph, d)
			

			surrreal =(phi[j]*X)+((sX**delta[j])*W)
			surrreal = surrreal.flatten()
			surrzero = 0.0
			surrchoice = [surrreal, surrzero]
			
			pnz = p[j]/lens
			pz = q[j]/lens
			select = np.random.choice(2, p=[pnz, pz])
			
			#print("Val:",surrreal,"Beta:",scl[j],"alpha:",shp[j],"W:",W, "b[1]<b[0]:",scl[j]<scl[j-1])
			
			surr[m,j] = surrchoice[select]
			betaout[m,j] = 1/scl[j]
			Wout[m,j] = W
	
	return surr, betaout, Wout
	
	
	
#----------------------------------------------------------------------------
'''
pgar_noise returns the noise component of the PGAR(1) process
as in Equation 20 (variable, W) in Fernandez and Salas (1986). 

Inputs parameters
----------
Two element arrays of the following parameters. Each element should be for time t
and time t-1.
alpha, beta, lam, phi, delta (shape, scale, location, periodic auto-regressive
coefficient and the periodic autoregressive exponent).

Returns
-------
A noise variable, W, for the PGAR(1) process.
    
Code development
-------
20160523 - script developed by Ailie Gallant
'''


def pgar_noise(alph, beta, lam, phi, delt) :
	import numpy as np
	import scipy.stats as stats

	if beta[1] < beta[0]:
		const = (alph[1]/(alph[0]**delt[1]))
		sone = beta[1]
		stwo = beta[0] - sone
		Y = (stats.beta.rvs(sone, stwo))**delt[1]
		
		gamma = -0.577216
		
		#Compute values for marginal probability, start with prob of delta and 1-delta
		prd = delt[1]
		prnd = 1 - delt[1]
		#Compute V' under marginal probability (Eq 13b)
		vdashd = gamma*(1-delt[1])
		vdashnd = gamma*(1-delt[1]) - (stats.expon.rvs(scale = 1/beta[1]))
		
		#Comput V under marginal probability (Eq 13c), iterate until tolerance at a low level
		vtol = 10
		vnum = 1
		vd = (1 - delt[1])/vnum
		vnd = ((1 - delt[1])/vnum) - stats.expon.rvs(scale = 1/(beta[1]+vnum))
		
		td = 0
		tnd = 0
		tol=0
		iter = 0
		while tol < 1:
			voldd = vd
			volnd = vnd
			
			vnum = vnum+1
			
			vd = voldd + ((1 - delt[1])/vnum)
			vnd = volnd + (((1 - delt[1])/vnum) - stats.expon.rvs(scale = 1/(beta[1]+vnum)))
			
			if abs(voldd-vd) < 0.01: td = td+1
			if abs(volnd-vnd) < 0.01: tnd = tnd+1
			if (td > 0) and (tnd > 0): tol=1
			iter = iter+1
		
		
		
		#Randomly select V' based on marginal probability
		sel = np.random.choice(2, p=[prd, prnd])
		vdash = [vdashd, vdashnd]
		vdash = vdash[sel]
		v = [vd, vnd]
		v = v[sel]
		
		Z = np.exp((vdash + v))
		
		#Compute vtau (note the delta power for Y has already been included)
		W = const*Z*Y
			
	else:
		
		#Compute e(0) as a three-parameter gamma variable with shape aez, scale bez
		#and location gez
		aez = beta[1] - beta[0]
		bez = alph[1]
		gez = lam[1] - (phi[1]*lam[0])
		ezero = stats.gamma.rvs(aez, scale = bez , loc = gez)
		
		#Compute innovation variable e(1), computed using Poisson random variable
		#theta.
		logval = phi[1]*(alph[0]/alph[1])
		if logval <= 0:
			theta = 0
		else:
			theta=-1*(beta[0]*np.log(logval))
		
		if theta < 0.0: theta = 0.0
		
		npoiss = stats.poisson.rvs(theta)
		
		if npoiss > 0:
			u = np.random.uniform(0,1)
			eone = (stats.expon.rvs(scale=alph[1]))*((phi[1]*(alph[0]/alph[1])**u))
			for o in range(0,npoiss):
				u = np.random.uniform(0,1)
				eone = eone + ((stats.expon.rvs(scale=alph[1]))*((phi[1]*(alph[0]/alph[1])**u)))
		else: eone = 0.0

		
		W = ezero + eone
	
	return W