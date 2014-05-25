import myhazmaps as mhp
import hazmap as hmp
import BASScast as bcp
import eqcatalog as eqp
import ANSStools as atp
import matplotlib.pyplot as plt
import matplotlib.dates as mpd
import mpl_toolkits.basemap as bmp
import scipy
import operator
import math
import glob
import os
import datetime as dtm
import pytz

from geographiclib.geodesic import Geodesic as ggp
import rbTectoFigs as rfp
import linefit

nContours1=15


def get_rbratios_tohoku(mc=5.0, todt=mhp.dtm.datetime.now(mhp.pytz.timezone('UTC')), winlen=None, avlen=None, targmag=9.0, ndithers=10, nContours=nContours1, bigmag=7.0, lons=[135., 148.5], lats=[30., 45.25], refreshcat=True, dt0=mhp.dtm.datetime(1990,1,1, tzinfo=mhp.pytz.timezone('UTC')), catname='cats/tohoku_rfits.cat', mt=7.6):
	#
	if refreshcat==True:
		#cl1=yp.catfromANSS(lon=[92.0, 106.0],lat=[-9.0, 10.0], minMag=3.5, dates0=[yp.dtm.datetime(1990,1,1, tzinfo=pytz.timezone('UTC')), yp.dtm.datetime(2012,6,19, tzinfo=pytz.timezone('UTC'))], fout=catname)
		cl1=atp.catfromANSS(lon=lons,lat=lats, minMag=mc, dates0=[dtm.datetime(1990,1,1, tzinfo=pytz.timezone('UTC')), dtm.datetime.now(pytz.timezone('UTC'))], fout=catname)
	#
	c1=eqp.eqcatalog([])
	c1.mt=mt
	c1.loadCatFromFile(catname)
	c1.cat.sort(key = lambda x:x[0])
	dlambda=1.76
	N_sequence = rfp.getNofm(targmag, mc)
	if avlen==None: avlen = max(1,int(N_sequence/10))
	tohokuEvent=c1.getMainEvent()
	print "n_seq: %d" % N_sequence
	tohokuLatLon = tohokuEvent[1], tohokuEvent[2]
	Lr = 10.0**(.5*targmag - 1.76)
	lfactor=.5
	#
	c1.subcats+=[['r1', rfp.circularcat(c1.getcat(0), latlon=[tohokuEvent[1], tohokuEvent[2]], Rkm=Lr*lfactor)]]
	#
	#rbratios = c1.plotIntervalRatiosAx(winlen=N_sequence, cat=c1.getcat(1), hitThreshold=1.0, bigmag=9.0, thisAx=myaxes[i], ratios=None, delta_t=1, avlen=rbavelen, mainEv=mainshock, logZ=None, rbLegLoc=0, reverse=False)
	rbratios = c1.plotIntervalRatiosAx(winlen=N_sequence, cat=c1.getcat(1), hitThreshold=1.0, bigmag=9.0, thisAx=None, ratios=None, delta_t=1, avlen=avlen, mainEv=tohokuEvent, logZ=None, rbLegLoc=0, reverse=False)
	plt.figure(1)
	log_fact = math.log10(N_sequence)
	#
	#plt.plot(map(operator.itemgetter(0), rbratios), map(operator.itemgetter(5), rbratios), '-')
	plt.semilogy(map(operator.itemgetter(1), rbratios), [x[4]*1.0 for x in rbratios], 'k-', alpha=.8, zorder=11)
	
	return rbratios

def lin_funct(X=None, a=0., b=1.):
	# a linear function for fitting.
	#
	#X=scipy.array(X)
	return a + b+X

def get_ratio_fits(ratios=None, Nfits=[1,5], ratiocol=4, meancol=5, t_col=1):
	# ratios: return from rbratios (or any time series)
	# Nfits=[] list of Ns over which to do fitting
	# ratiocol: data column
	#
	# fits will include: mean value, linear fits {chi-square, a, b, sig_a, sig_b? (maybe later}
	# append these vals in a dict. to the end of each row: [{}, {}, {}..]
	#
	r_vals = map(operator.itemgetter(ratiocol), ratios)
	r_vals_mean = map(operator.itemgetter(meancol), ratios)
	X_vals = map(operator.itemgetter(t_col), ratios)
	fitsets = {}
	for n in Nfits:
		fitsets[n]={'means':[], 'fit_prams':[], 'mean_fit_prams':[], 'fitses':[]}
		for i in xrange(n, len(ratios)):
			these_r = r_vals[(i-n):i]
			these_rm = r_vals_mean[(i-n):i]
			these_x = map(mpd.date2num, X_vals[(i-n):i])
			these_p = (0., 0.)
			#
			meanval = scipy.mean(these_r)	# or, it's faster to do a moving mean; mean+=(x_i-x_{i-n})
			# now, get a line-fit for this segment...
			#
			'''
			fit_return = scipy.optimize.curve_fit(lin_funct, these_x, these_r, p0=these_p)
			#print "fit_return: ", fit_return
			#print fit_return[1]
			fit_vals = fit_return[0]
			cov = fit_return[1]/float(n) #scipy.trace(fit_return[1])
			#
			these_p = (0., 0.)
			fit_return_mean = scipy.optimize.curve_fit(lin_funct, these_x, these_rm, p0=these_p)
			mean_fit_vals = fit_return_mean[0]
			cov_mean = fit_return_mean[1]/float(n) #scipy.trace(fit_return_mean[1])
			y_theory = [lin_funct(x, mean_fit_vals[0], mean_fit_vals[1]) for x in these_x]
			chisqr = (scipy.sum(scipy.array(these_rm)-scipy.array(y_theory))/(n-2.0))**2.0
			#
			'''
			lf=linefit.linefit([these_x, these_rm])
			lf.doFit()
			
			fitsets[n]['means']          += [meanval]
			#fitsets[n]['fit_prams']      += [list(fit_vals) + [cov]]
			#fitsets[n]['mean_fit_prams'] += [list(mean_fit_vals) + [chisqr]]
			fitsets[n]['fitses'] += [[lf.a, lf.b, lf.meanVar()]]
			#
			#print 'fitses: ', fitsets[n]['fit_prams'][-1]
		#
	#
	return fitsets
	#
	

