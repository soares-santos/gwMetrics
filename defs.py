import math
import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt
import scipy
import scipy.interpolate
import scipy.optimize
import emcee
import scipy.stats as stats
from palettable.colorbrewer.qualitative import Paired_11
colors = Paired_11.hex_colors
import healpy as hp

def meanProb(mjd_obs, ra_obs, dec_obs, band_obs, event_prob_maps, event_mjds, delta_t=10.,band=None):

# check length of input gw event arrays
    nevents = len(event_prob_maps)
    if len(event_mjds) != nevents : 
        print "Error: event_mjds and event_prob_maps must have the same length." 
        exit

# set radius to camera's FOV diameter * 0.5
    radius=np.deg2rad(3.5*0.5)
    print radius

# add prob over all events
    total_prob = 0.
    for i in np.arange(nevents):
        nside=hp.npix2nside(len(event_prob_maps[i]))
        ix = ((mjd_obs < event_mjds[i] + delta_t) & (mjd_obs > event_mjds[i]))
        if len(ix) > 0:
            theta = 0.5 * np.pi - dec_obs[ix] 
            phi = ra_obs[ix]
            vec = hp.ang2vec(theta,phi)
            hpix = np.zeros(0,dtype=np.int)
            for j in np.arange(len(vec)):
                 hpix = np.append(hpix,hp.query_disc(nside,vec[j],radius))
            hpix = np.unique(hpix)
            prob = event_prob_maps[i][hpix]
            total_prob = total_prob + prob.sum()
            print i, prob.sum(), total_prob

            # make map of number of observations
            obs_nside = 16
            obs_map = np.zeros(hp.nside2npix(obs_nside),dtype=np.int) 
            hpix_obs = hp.ang2pix(obs_nside,theta,phi)
            for pix in hpix_obs: 
                obs_map[pix] = obs_map[pix] + 1
                ##print pix, obs_map[pix]
###            hp.mollview(obs_map)

    print total_prob, nevents, total_prob/nevents
    return total_prob/float(nevents)

def readFiles(spectraFile, transmissionFile):
    #open quasar spectrum and transmission curves
    spectra_data = ascii.read(spectraFile, delimiter = ',')
    trans_data = ascii.read(transmissionFile, delimiter = ',')
    
    #transmission wavelengths for all filters
    wl2 = trans_data['wavelength']
    
    #responses for all filters
    resp = trans_data['trans1p3p']
    
    #separate transmission wavelengths by filter
    wlresp = {'u': wl2[0:46], 'g' : wl2[47:135], 'r' : wl2[136:210], 'i' : wl2[211:299], 'z' : wl2[300:440]}
    #separate responses by filter
    response = {'u': resp[0:46], 'g' : resp[47:135], 'r' : resp[136:210], 'i' : resp[211:299], 'z' : resp[300:440]}
    
    return spectra_data['wl'], spectra_data['fl'], wlresp, response

def calcLambdaEff(band, z, wl, fl, wlresp, response):
    #takes the wavelength and flux for composite spectra,
    #wavelength and response for the filters
    #and the redshift to slide the composite spectra to
    #output the effective wavelength for the composite spectra in each filter band
    
    eff_band = []
    
    # Shift the spectra
    wl_z = wl*(1.0 + z)
    
    #Cut wavelength and flux for composite spectra to make the filter curves
    wl_cut = wl_z[(wl_z >= min(wlresp[band])) & (wl_z <= max(wlresp[band]))]
    fl_cut = fl[(wl_z >= min(wlresp[band])) & (wl_z <= max(wlresp[band]))]
        
    #Defines a function that will interpolate these x and y points
    # x is the wavelengths of the filter curve
    # y is the responses of the filter curve
    interp = scipy.interpolate.interp1d(wlresp[band], response[band])
        
    #Interpolate the above y points on this new x grid
    #new x grid is the wavelengths of the composite spectra
    newresp_band = interp(wl_cut)
        
    #Equation 4 in Kacz+09:
    
    #Numerator for effective wavelength
    # (filter response)*(flux)*ln(wavelength)
    num = newresp_band * fl_cut * np.log(wl_cut)
        
    #Denominator for effective wavelength
    # (filter response)*(flux)
    den = newresp_band * fl_cut
    
    #Reimann sum to get effective wavelength
    try:
        eff = np.exp(sum(num) / sum(den))
    except ZeroDivisionError:
         eff = np.nan
    
    #append result for each bandpass to output list 
    eff_band.append(eff)
    
    return np.array(eff_band)

def calcLambdaEffPL(band, pl, wlresp, response):
    #takes the wavelength and response for the filters
    #output the effective wavelength for a powerlaw in each filter band
    
    #Power-law spectrum (independent of redshift)
    wl = np.arange(2500,12000,1) #wavelengths
    fl = wl**pl #fluxes

    eff_band = []
    
    #Cut wavelength and flux vectors to make the filter curves
    wl_cut = wl[(wl >= min(wlresp[band])) & (wl <= max(wlresp[band]))]
    fl_cut = fl[(wl >= min(wlresp[band])) & (wl <= max(wlresp[band]))]
    
    #Defines a function that will interpolate these x and y points
    # x is the wavelengths of the filter curve
    # y is the responses of the filter curve
    interp = scipy.interpolate.interp1d(wlresp[band], response[band])
    
    #Interpolate the above y points on this new x grid
    #new x grid is the wavelengths of the composite spectra
    newresp_band = interp(wl_cut)
    
    #Equation 4 in Kacz+09:
    
    #Numerator for effective wavelength
    # (filter response)*(flux)*ln(wavelength)
    num = newresp_band * fl_cut * np.log(wl_cut)
    
    #Denominator for effective wavelength
    # (filter response)*(flux)
    den = newresp_band * fl_cut
    
    #Reimann sum to get effective wavelength
    eff = np.exp(sum(num) / sum(den))
    
    #append result for each bandpass to output list 
    eff_band.append(eff)
    
    return np.array(eff_band)
    
def calcR0(eff_band):
    #takes the effective wavelength in each filter
    #output R0, Equation 2 in Kacz+09
    
    #turn the list into numpy array
    leff = eff_band/1.0e4
    
    #Equation 3 in Kacz+09
    nm1e6 = 64.328 + 29498.1/(146.0-(1.0/leff)**2) + 255.4/(41.0-(1.0/leff)**2)
    
    #index of refraction that goes into Equation 2 of Kacz+09
    n = 1.0 + nm1e6/1.0e6
    
    #Equation 2 in Kacz+09
    #Needed to plug into Equation 1
    R0 = (3600.0*180.0/3.14159)*(n**2.0-1.0)/2.0*n**2.0
    
    return R0

def calcR(airmasses, filters, zshift=1.0):
    #takes array of airmasses and array of filters
    #returns the angular deflection
    
    #composite quasar spectrum
    spectraFile = 'datafiles/compositequasar.csv'
    #these are the SDSS filter curves ['u', 'g', 'r', 'i', 'z']
    transmissionFile = 'datafiles/alltransfinal.csv'
    wl, fl, wlresp, response = readFiles(spectraFile, transmissionFile)
    
    #powerlaw slope
    #same as the slope of the composite spectra's continuum
    pl=-1.5
    
    #redshifts of the composite quasar spectra
    #z_steps = np.arange(0.0, 4.6, 0.01)
    bands = ['u', 'g', 'r', 'i', 'z']
    
    Rdiffs = np.array([])
    #Rdiffs = np.append(Rdiffs, z_steps)
    #Determine offset
    for band in bands:
        R0pl = calcR0(calcLambdaEffPL(band, pl, wlresp, response))
        Rdiff = calcR0(calcLambdaEff(band, zshift, wl, fl, wlresp, response))-R0pl
        #R0spec = calcR0(calcLambdaEff(band, z_steps, wl, fl, wlresp, response))
        #Rdiff = R0spec[filter]-R0pl[filter]
        Rdiffs = np.append(Rdiffs, Rdiff)
    
    #print np.shape(Rdiffs)
    #Rdiffs = Rdiffs.reshape(len(bands)+1, -1)
    #print np.shape(Rdiffs)
    #print Rdiffs[0]
    #print np.asarray(Rdiffs[0])
    
    #print zip(['steps','u', 'g', 'r', 'i', 'z'],['float64']*(len(bands)+1))
    
    Rdiffs_rec = np.rec.fromarrays((Rdiffs[0], Rdiffs[1], Rdiffs[2], Rdiffs[3], Rdiffs[4]), dtype=zip(['u', 'g', 'r', 'i', 'z'],['float64']*len(bands))).copy()
    
    #print Rdiffs_rec['steps']
    
    #Determine zenith distance
    #Angle is (180.0/3.14159)*Z
    #Z = np.arccos(1.0/AM)
    
    #print np.shape(airmasses)
    #print np.shape(filters)
    #print np.shape(zip(airmasses, filters))
    
    #for (am ,f) in zip(airmasses, filters):
    #    print f, am
    #    print Rdiffs_rec[f], np.tan(np.arccos(1.0/am))
    #    print (Rdiffs_rec[f])*np.tan(np.arccos(1.0/am))
    
    tanZList = np.asarray([np.tan(np.arccos(1.0/am)) for am in airmasses])
    RList = np.asarray([(Rdiffs_rec[f])*np.tan(np.arccos(1.0/am)) for (am ,f) in zip(airmasses, filters)])
    #R = (Rdiffs_rec[filter])*np.tan(Z)
    return tanZList, RList

def fit_DCR_params(airmasses, filters, zshift=2.1):
	tanZList, RList = calcR(airmasses, filters, zshift=zshift)
	n = len(tanZList)
	print "number of visits: ", n
	def lnlike(theta, x, y, yerr):
		m, lnf = theta
		model = m * x + 0.0
		inv_sigma2 = 1.0/(yerr**2. + model**2.*np.exp(2.*lnf))
		return -0.5*(np.sum(((y-model)**2.*inv_sigma2 - np.log(inv_sigma2))))
		#return -0.5*(np.sum(((y-model)**2.)))
	def lnprior(theta):
		m, lnf = theta
		if (-1.0 < m < 1.0) and (-350.0 < lnf < 350.0):
			return 0.0
		return -np.inf
	def lnprob(theta, x, y, yerr):
		lp = lnprior(theta)
		if not np.isfinite(lp):
			return -np.inf
		return lp + lnlike(theta, x, y, yerr)
	nll = lambda *args: -lnprob(*args)
	x = np.copy(tanZList)
	y = np.copy(RList)
	yerr = np.sqrt((0.02**2.)+(0.02**2.))
	offset = yerr * np.random.randn(n) + 0.0
	pm = np.random.choice([-1.0,1.0], size=n, replace=True)
	y += offset
	result = scipy.optimize.minimize(nll, [-0.001, np.log(0.5)], args=(x, y, yerr), method="Nelder-Mead")
	m_ml, lnf_ml = result["x"]
	ndim, nwalkers = 2, 100
	pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
	sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr))
	sampler.run_mcmc(pos, 500)
	samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
	ms = samples[np.random.randint(len(samples), size=100)][:,0]
	m_mcmc, lnf_mcmc = map(lambda v: (v[2]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
	plt.plot(tanZList, RList, 'o')
	plt.errorbar(x, y, yerr=yerr, fmt='.')
	xs = np.arange(min(tanZList), max(tanZList), 0.01)
	plt.plot(xs, m_ml*xs + 0.0, color="r")
	for m, lnf in samples[np.random.randint(len(samples), size=100)]:
		plt.plot(xs, m*xs + 0.0, color="k", alpha=0.1)
	plt.xticks(np.arange(-0.5, 2.25, 0.25),size=14)
	plt.yticks(np.arange(-0.2, 0.25, 0.05),size=14)
	plt.xlim(-0.01,1.75)
	plt.ylim(-0.2, 0.2)
	plt.xlabel(r'$\tan (Z)$',size=14)
	plt.ylabel(r'$\Delta R_{||}$ (arcsec)',size=14)
	plt.savefig('auPar_tanZ_' + str(zshift) + '_' + str(n) + '.png')
	plt.clf()
	#b_ls, m_ls = fit_matrix(x, y, yerr)
	return x, y, yerr, ms, m_ml, m_mcmc, lnf_mcmc

def calculate_statistic(airmasses, filters, zshift1=0.8, zshift2=2.2):
	print zshift1, zshift2
	x1, y1, yerr1, ms1, m_ml1, m_mcmc1, quartiles1 = fit_DCR_params(airmasses, filters, zshift=zshift1)
	x2, y2, yerr2, ms2, m_ml2, m_mcmc2, quartiles1 = fit_DCR_params(airmasses, filters, zshift=zshift2)
	statistic, pvalue = stats.ks_2samp(ms1,ms2)
	return pvalue

def bayesian_odds_ratio(airmasses, filters, astrometric_error=0.020, zshift=2.1):
	plot_points = True
	plot_walkers = False
	intercept_fixed = True
	np.random.seed(0)
	tanZList, RList = calcR(airmasses, filters, zshift=zshift)
	n = len(tanZList)
	def lnlike(theta, x, y, yerr, type):
		if type=="flat":
			b = theta
			model = 0.0 * x + b
		if type=="slope":
			if intercept_fixed == True:
				m = theta
				model = m * x + 0.0
			else:
				m, b = theta
				model = m * x + b
		inv_sigma2 = 1.0/(yerr**2.)
		return -0.5*(np.sum(((y-model)**2.*inv_sigma2 - np.log(inv_sigma2))))
	def lnprior(theta, type):
		if type=="flat":
			b = theta
			if (-1.0 < b < 1.0):
				return 0.0
			return -np.inf
		if type=="slope":
			if intercept_fixed == True:
				m = theta
				if (-1.0 < m < 1.0):
					return 0.0
				return -np.inf
			else:
				m, b = theta
				if (-1.0 < m < 1.0) and (-1.0 < b < 1.0):
					return 0.0
				return -np.inf
	def lnprob(theta, x, y, yerr, type=None):
		if type=="flat" or type=="slope":
			lp = lnprior(theta, type)
			if not np.isfinite(lp):
				return -np.inf
			return lp + lnlike(theta, x, y, yerr, type)
		else:
			print "must specify flat or slope"
			return np.nan
	nll = lambda *args: -lnprob(*args)
	nsteps, nwalkers = 500, 100
	x = np.copy(tanZList)
	y = np.copy(RList)
	yerr = np.sqrt((astrometric_error**2.)+(astrometric_error**2.))
	offset = yerr * np.random.randn(n) + 0.0
	pm = np.random.choice([-1.0,1.0], size=n, replace=True)
	y += offset
	if plot_points == True:
		fig1 = plt.figure(1)
		#plt.plot(tanZList, RList, 'o', color=colors[0])
		plt.errorbar(x, y, yerr=yerr, fmt='.', color=colors[1])
	if intercept_fixed == True:
		ndim = 1
		result = scipy.optimize.minimize(nll, [-0.001], args=(x, y, yerr, "slope"), method="Nelder-Mead")
		m_ml = result["x"]
	else:
		ndim = 2
		result = scipy.optimize.minimize(nll, [-0.001, 0.0], args=(x, y, yerr, "slope"), method="Nelder-Mead")
		m_ml, b_ml = result["x"]
	pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
	sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr, "slope"))
	sampler.run_mcmc(pos, nsteps)
	samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
	ms = samples[np.random.randint(len(samples), size=100)][:,0]
	if intercept_fixed == True:
		m_mcmc_slope = map(lambda v: (v[1]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
	else:
		m_mcmc_slope, b_mcmc_slope = map(lambda v: (v[1]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
	if plot_walkers == True:
		fig2 = plt.figure(2)
		ax1 = plt.subplot(211)
		for i in range(nwalkers):
			ax1.plot(sampler.chain[i,:,0],color='k',alpha=0.05)
			ax1.axhline(y=m_mcmc_slope, xmin=0, xmax=nsteps, color='r')
		#ax2 = plt.subplot(412)
		#for i in range(nwalkers):
		#	ax2.plot(sampler.chain[i,:,1],color='k',alpha=0.05)
		#	ax2.axhline(y=lnf_mcmc_slope, xmin=0, xmax=nsteps, color='r')
	xs = np.arange(min(tanZList), max(tanZList), 0.01)
	if plot_points == True:
		plt.figure(1)
		if intercept_fixed == True:
			plt.plot(xs, m_mcmc_slope*xs + 0.0, color=colors[2], lw=3)
			for m in samples[np.random.randint(len(samples), size=100)]:
				plt.plot(xs, m*xs + 0.0, color=colors[2], lw=1, alpha=0.2)
		else:
			plt.plot(xs, m_mcmc_slope*xs + b_mcmc_slope, color=colors[2], lw=3)
			for m, b in samples[np.random.randint(len(samples), size=100)]:
				plt.plot(xs, m*xs + b, color=colors[2], lw=1, alpha=0.2)
	ndim = 1
	result = scipy.optimize.minimize(nll, [0.0], args=(x, y, yerr, "flat"), method="Nelder-Mead")
	b_ml = result["x"]
	pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
	sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr, "flat"))
	sampler.run_mcmc(pos, nsteps)
	samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
	ms = samples[np.random.randint(len(samples), size=100)][:,0]
	b_mcmc_flat = map(lambda v: (v[1]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
	xs = np.arange(min(tanZList), max(tanZList), 0.01)
	if plot_points == True:
		plt.figure(1)
		plt.plot(xs, 0.0*xs + b_mcmc_flat, color=colors[4], lw=3)
		for b in samples[np.random.randint(len(samples), size=100)]:
			plt.plot(xs, 0.0*xs + b, color=colors[4], lw=1, alpha=0.2)
	if plot_walkers == True:
		plt.figure(2)
		ax3 = plt.subplot(212)
		for i in range(nwalkers):
			ax3.plot(sampler.chain[i,:,0],color='k',alpha=0.05)
			ax3.axhline(y=b_mcmc_flat, xmin=0, xmax=nsteps, color='r')
		#ax4 = plt.subplot(414)
		#for i in range(nwalkers):
		#	ax4.plot(sampler.chain[i,:,1],color='k',alpha=0.05)
		#	ax4.axhline(y=lnf_mcmc_flat, xmin=0, xmax=nsteps, color='r')
		plt.savefig('walkers_modelcompariosn_test_' + filters[0] + '_' + str(max(airmasses)-min(airmasses)).replace(".", "") + '_' + str(astrometric_error).replace(".", "")  + '_' + str(zshift).replace(".", "")  + '_' + str(n) + '.png')
		plt.clf()
	if intercept_fixed == True:
		model_slope = m_mcmc_slope * x + 0.0
	else:
		model_slope = m_mcmc_slope * x + b_mcmc_slope
	model_flat = 0.0 * x + b_mcmc_flat
	inv_sigma2 = 1.0/(yerr**2.)
	slope_loglikelihood = (-0.5*(np.sum(((y-model_slope)**2.*inv_sigma2 - np.log(inv_sigma2)))))
	flat_loglikelihood = (-0.5*(np.sum(((y-model_flat)**2.*inv_sigma2 - np.log(inv_sigma2)))))
	bayes_ratio = np.e**(slope_loglikelihood - flat_loglikelihood)
	if plot_points == True:
		plt.figure(1)
		#plt.text(0.1, 0.05, str(bayes_ratio), ha='left', va='center')
		#plt.text(0.1, 0.00, str(slope_loglikelihood), ha='left', va='center')
		#plt.text(0.1, -0.05, str(flat_loglikelihood), ha='left', va='center')
		plt.xticks(np.arange(-0.5, 2.5, 0.25),size=14)
		plt.yticks(np.arange(-0.4, 0.4, 0.10),size=14)
		plt.xlim(-0.01,2.3)
		plt.ylim(-0.35, 0.35)
		plt.xlabel(r'$\tan (Z)$',size=14)
		plt.ylabel(r'$\Delta R_{||}$ (arcsec)',size=14)
		plt.savefig('offset_tanZ_modelcompariosn_' + filters[0] + '_' + str(max(airmasses)-min(airmasses)).replace(".", "") + '_' + str(astrometric_error).replace(".", "")  + '_' + str(zshift).replace(".", "")  + '_' + str(n) + '.png')
		plt.clf()
	if not np.isfinite(bayes_ratio):
		print "Redshift: ", zshift, " Airmass Range: ", min(airmasses), " - ", max(airmasses), " Number of Observations: ", n, " Astrometric Error: ", astrometric_error
		print min(x), max(x), min(y), max(y)
		print slope_loglikelihood, flat_loglikelihood
		print bayes_ratio
	return bayes_ratio
