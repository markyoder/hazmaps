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

def doTohoku():
	A = get_rbratios_tohoku()
	B = get_ratio_fits(ratios=A, Nfits=[22])
	C = plot_ratio_fits(B, new_figs=False, fignum0=2)
	#
	return (A,B,C)

def get_rbratios_tohoku(mc=5.0, todt=mhp.dtm.datetime.now(mhp.pytz.timezone('UTC')), winlen=None, avlen=None, targmag=9.0, ndithers=10, nContours=nContours1, bigmag=7.0, lons=[135., 148.5], lats=[30., 45.25], refreshcat=True, dt0=mhp.dtm.datetime(1990,1,1, tzinfo=mhp.pytz.timezone('UTC')), catname='cats/tohoku_rfits.cat', mt=7.6):
	#
	if refreshcat==True:
		#cl1=yp.catfromANSS(lon=[92.0, 106.0],lat=[-9.0, 10.0], minMag=3.5, dates0=[yp.dtm.datetime(1990,1,1, tzinfo=pytz.timezone('UTC')), yp.dtm.datetime(2012,6,19, tzinfo=pytz.timezone('UTC'))], fout=catname)
		cl1=atp.catfromANSS(lon=lons,lat=lats, minMag=mc, dates0=[dtm.datetime(1990,1,1, tzinfo=pytz.timezone('UTC')), dtm.datetime.now(pytz.timezone('UTC'))], fout=catname)
	#
	plt.figure(1)
	plt.clf()
	
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
		fitsets[n]={'means':[], 'fit_prams':[], 'mean_fit_prams':[], 'fitses':[], 'X':X_vals}
		for i in xrange(n, len(ratios)):
			these_r = r_vals[(i-n):i]
			these_rm = r_vals_mean[(i-n):i]
			these_x = map(mpd.date2num, X_vals[(i-n):i])
			these_p = (0., 0.)
			#
			#meanval = scipy.mean(these_r)	# or, it's faster to do a moving mean; mean+=(x_i-x_{i-n})
			meanval = these_rm[-1]
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
			#
			
			lf  = linefit.linefit([these_x, these_r])
			lf2 = linefit.linefit([these_x, these_rm])
			lf.doFit()
			lf2.doFit()
			#
			fitsets[n]['means']          += [meanval]
			#fitsets[n]['fit_prams']      += [list(fit_vals) + [cov]]
			#fitsets[n]['mean_fit_prams'] += [list(mean_fit_vals) + [chisqr]]
			fitsets[n]['fit_prams'] += [[lf.a, lf.b, lf.meanVar()]]
			fitsets[n]['mean_fit_prams'] += [[lf2.a, lf2.b, lf2.meanVar()]]
			#
			#print 'fitses: ', fitsets[n]['fit_prams'][-1]
		#
	#
	return fitsets
	#
#
def get_averaging_set(ratios=None, ave_lens=[5, 500], x_col=1, ratio_col=4):
	# get plot of fit quality as a function of ave-len
	if ratios==None: return None
	#
	if ave_lens[0]>ave_lens[1]: ave_lens.reverse()
	#
	r_vals = map(operator.itemgetter(ratio_col), ratios)
	x_vals = map(operator.itemgetter(x_col), ratios)
	if hasattr(x_vals[0], 'year'):
		# it's a date type.
		x_vals = map(mpd.date2num, x_vals)
	#
	chi_sqrs = []		# will be a list of lists [[chi-sqr_r, aic, totalVar, ndof]]
	total_chi_sqrs=[]	# [ [ave_len, mean_chi_sqr_r, mean_aic, N_points] ]
	#
	ave_len = ave_lens[0]
	for ave_len in xrange(*ave_lens):
		mean_r_vals = [scipy.mean(r_vals[i-ave_len:i]) for i in xrange(ave_len, len(r_vals))]
		chi_sqrs += [ [] ]
		#aic_vals += [ [] ]
			
		for i in xrange(ave_len, len(mean_r_vals)):
			# note: 2*ave_len so we have a full set of means.
			#these_r = r_vals[(i-ave_len):i]
			these_r = mean_r_vals[(i-ave_len):i]
			these_x = x_vals[i:i+ave_len]
			these_p = (0., 0.)
			#
			#meanval = scipy.mean(these_r)	# or, it's faster to do a moving mean; mean+=(x_i-x_{i-n})
			
			# now, get a line-fit for this segment...
			#
			#
			lf  = linefit.linefit([these_x, these_r])
			lf.doFit()
			#
			chi_sqrs[-1]+=[ [lf.meanVar(), lf.AIC, lf.totalVar, lf.Ndof] ]
		print "fit(%d)" % (ave_len)
		#
		my_arrays = zip(*chi_sqrs[-1])
		total_chi_sqrs += [ [ave_len, scipy.mean(my_arrays[0]), scipy.mean(my_arrays[1]), len(chi_sqrs[-1])] ]
		#
	#
	return [total_chi_sqrs, chi_sqrs]
#
def get_averaging_set2(ratios=None, ave_lens=[5, 100], x_col=1, ratio_col=4):
	# get plot of fit quality as a function of ave-len
	if ratios==None: return None
	#
	if ave_lens[0]>ave_lens[1]: ave_lens.reverse()
	#
	r_vals = map(operator.itemgetter(ratio_col), ratios)
	x_vals = map(operator.itemgetter(x_col), ratios)
	if hasattr(x_vals[0], 'year'):
		# it's a date type.
		x_vals = map(mpd.date2num, x_vals)
	#
	chi_sqrs = []		# will be a list of lists [[chi-sqr_r, aic, totalVar, ndof]]
	total_chi_sqrs=[]	# [ [ave_len, mean_chi_sqr_r, mean_aic, N_points] ]
	fit_len=22
	#
	ave_len = ave_lens[0]
	for ave_len in xrange(*ave_lens):
		mean_r_vals = [scipy.mean(r_vals[i-ave_len:i]) for i in xrange(ave_len, len(r_vals))]
		chi_sqrs += [ [] ]
		#aic_vals += [ [] ]
		#
		#for i in xrange(ave_len, len(mean_r_vals)):
		for i in xrange(fit_len, len(mean_r_vals)):
			# note: 2*ave_len so we have a full set of means.
			these_r = mean_r_vals[(i-fit_len):i]
			these_x = x_vals[(i + ave_len - fit_len):i+ave_len]
			these_p = (0., 0.)
			#
			#meanval = scipy.mean(these_r)	# or, it's faster to do a moving mean; mean+=(x_i-x_{i-n})
			
			# now, get a line-fit for this segment...
			#
			#
			lf  = linefit.linefit([these_x, these_r])
			lf.doFit()
			#
			chi_sqrs[-1]+=[ [lf.meanVar(), lf.AIC, lf.totalVar, lf.Ndof] ]
		print "fit(%d)" % (ave_len)
		#
		my_arrays = zip(*chi_sqrs[-1])
		total_chi_sqrs += [ [ave_len, scipy.mean(my_arrays[0]), scipy.mean(my_arrays[1]), len(chi_sqrs[-1])] ]
		#
	#
	return [total_chi_sqrs, chi_sqrs]

#
def plot_total_chi_sqrs(totals_list=None, ave_lens=[5,100], x_col=1, ratio_col=4, fignum=0):
	#
	# totals_list should come from the get_averaging_set[0]
	#
	if totals_list == None: totals_list = get_averaging_set(ratios=None, ave_lens=ave_lens, x_col=x_col, ratio_col=ratio_col)
	#
	my_arrays = zip(*totals_list)
	#
	plt.figure(fignum)
	plt.clf()
	plt.plot(my_arrays[0], my_arrays[1], '-', label='chi-squares')
	plt.plot(my_arrays[0], my_arrays[2], '-', label='aic vals')
	#
	plt.legend(loc=0, numpoints=1)
	

	
#
def plot_ratio_fits(fit_dict, new_figs=False, fignum0=2):
	#
	fnum=fignum0
	plt.figure(fnum)
	plt.clf()
	#
	for fit_set in fit_dict.keys():
		this_set=fit_dict[fit_set]
		X   = this_set['X']
		y_r = this_set['means']
		y_a = map(operator.itemgetter(0), this_set['fit_prams'])
		y_b = map(operator.itemgetter(1), this_set['fit_prams'])
		y_chi = map(operator.itemgetter(2), this_set['fit_prams'])
		#
		plt.plot(X[-len(y_r):], y_r, '-', label='ratios')
		plt.plot(X[-len(y_chi):], y_chi, '--', label='chisqr: %.4f/%.4f/%.4f' % (sum(y_chi), scipy.mean(y_chi), len(y_chi)))
		plt.legend(loc=0, numpoints=1)
		plt.title('raw fits')
		#
		'''
		fnum+=1
		plt.figure(fnum)
		plt.clf()
		plt.plot(X[-len(y_a):], y_a, '-')
		plt.title("a values")
		#
		fnum+=1
		plt.figure(fnum)
		plt.clf()
		plt.plot(X[-len(y_b):], y_b, '-')
		plt.title("b values")
		'''
		
		#
		if new_figs==True or 1==1:
			fnum+=1
			plt.figure(fnum)
			plt.clf()
		#
		this_set=fit_dict[fit_set]
		y_r = this_set['means']
		y_a = map(operator.itemgetter(0), this_set['mean_fit_prams'])
		y_b = map(operator.itemgetter(1), this_set['mean_fit_prams'])
		y_chi = map(operator.itemgetter(2), this_set['mean_fit_prams'])
		#
		#
		
		plt.plot(X[-len(y_r):], y_r, '-', label='ratios')
		plt.plot(X[-len(y_chi):], y_chi, '--', label='chisqr: %.4f/%.4f/%.4f' % (sum(y_chi), scipy.mean(y_chi), len(y_chi)))
		plt.legend(loc=0, numpoints=1)
		plt.title('mean-fits')
		#
		'''
		fnum+=1
		plt.figure(fnum)
		plt.clf()
		plt.plot(X[-len(y_a):], y_a, '-')
		plt.title("a values (mean fits)")
		#
		fnum+=1
		plt.figure(fnum)
		plt.clf()
		plt.plot(X[-len(y_b):], y_b, '-')
		plt.title("b values (mean fits)")
		'''
		#
	return None
	
	
	
	
	
	
	
	
	
	
	
	
	

