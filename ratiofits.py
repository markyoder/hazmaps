import myhazmaps as mhp
import hazmap as hmp
import BASScast as bcp
import eqcatalog as eqp
import ANSStools as atp
import matplotlib.pyplot as plt
import matplotlib.dates as mpd
import mpl_toolkits.basemap as bmp
import scipy
import numpy
import operator
import math
import glob
import os
import datetime as dtm
import pytz
import random

from geographiclib.geodesic import Geodesic as ggp
import rbTectoFigs as rfp
import linefit

nContours1=15

tohoku_prams = {'mc':5.0, 'todt':mhp.dtm.datetime.now(mhp.pytz.timezone('UTC')), 'winlen':None, 'avlen':None, 'targmag':9.0, 'ndithers':10, 'nContours':nContours1, 'bigmag':7.0, 'lons':[135., 148.5], 'lats':[30., 45.25], 'refreshcat':True, 'dt0':mhp.dtm.datetime(1990,1,1, tzinfo=mhp.pytz.timezone('UTC')), 'catname':'cats/tohoku_rfits.cat', 'mt':7.6}

sumatra_prams = {'mc':5.0, 'targmag':9.1, 'rbavelen':None, 'bigmag':9.5, 'intlist':None, 'catname':'cats/sumatra.cat', 'refreshcat':False, 'plotevents':False, 'mt':7.55, 'lons':[92.0, 106.0], 'lats':[-9.0, 10.0], 'lfactor':.5}

def doTohoku():
	A = get_rbratios_tohoku()		# returns a rbratios() set like [[index, dtm, r>, r<, ratio, <ratio>]]
	B = get_ratio_fits(ratios=A, Nfits=[22])
	C = plot_ratio_fits(B, new_figs=False, fignum0=2)
	#
	return (A,B,C)

def doSumatra(mc=5.0, targmag=9.1, rbavelen=None, bigmag=9.5, intlist=None, catname='cats/sumatra.cat', refreshcat=False, plotevents=False, mt=7.55, lons=[92.0, 106.0],lats=[-9.0, 10.0], lfactor=.5):
	# get a bunch of sumatra fits from rbTectoFigs. we'll want analysis on 0, 2, 3
	#sumatra_catalog = rfp.sumatraQuad()
	
	cl1=atp.catfromANSS(lon=lons,lat=lats, minMag=mc, dates0=[dtm.datetime(1990,1,1, tzinfo=pytz.timezone('UTC')), dtm.datetime.now(pytz.timezone('UTC'))], fout=None)
	catalogs=[]
	#
	c1=eqp.eqcatalog(cl1)
	c1.mc=mc
	c1.mt=mt
	c1.targmag = targmag
	#
	c1.cat.sort(key = lambda x:x[0])
	dlambda=1.76
	Lr=10.0**(targmag/2.0 - dlambda)
	latlonMS=[3.316, 95.854]	#mainshock
	#mainshock = [dtm.datetime(2004, 12, 26, 0, 58, 53, 449995, tzinfo=pytz.timezone('UTC')),  3.295,  95.982,  9.0,  30.0]
	mainshock=c1.getMainEvent()
	mainshock[3]=9.1		#ANSS seems to be listing sumatra as 9.0 lately.
	mevIndex = mainshock[-1]
	#
	foreshocks = []
	foreshocks+=[mainshock[:]]
	foreshocks+=[mainshock[:]]
	#foreshocks+=[[dtm.datetime(2000, 6, 4, 16, 28, 26, 170001, tzinfo=pytz.timezone('UTC')), -4.721, 102.087, 7.9, 33.0]]
	#foreshocks+=[[dtm.datetime(2001, 2, 13, 19, 28, 30, 260001, tzinfo=pytz.timezone('UTC')), -4.68, 102.562, 7.4, 36.0]]
	#foreshocks+=[[dtm.datetime(2002, 11, 2, 1, 26, 10, 699996, tzinfo=pytz.timezone('UTC')), 2.824, 96.085, 7.4, 30.0]]
	#
	# adjust mainshock:
	foreshocks[0][1]+=0.
	foreshocks[0][2]-=0.
	#
	for i in xrange(mevIndex+1,len(c1.getcat(0))):
		if c1.getcat(0)[i][3]>=8.0:
			foreshocks+=[c1.getcat(0)[i]]
			print "large aftershock: ", c1.getcat(0)[i]
		
	#
	interesting_quads = [0,2,3]		# we know these are interesting...
	interesting_quads = range(len(foreshocks))
	for i in interesting_quads:
		#thiscatnum=i
		#thisfignum=i
		thismag=foreshocks[i][3]
		thisLr = 10.0**(thismag/2.0 - dlambda)
		#
		#if i>0:
		fsLatLon=[foreshocks[i][1], foreshocks[i][2]]
		#if i>0: c1.subcats+=[['r%d' %i, circularcat(c1.getcat(0), latlon=fsLatLon, Rkm=Lr*lfactor)]]
		catalogs+=[eqp.eqcatalog(rfp.circularcat(c1.getcat(0), latlon=fsLatLon, Rkm=Lr*lfactor))]
		catalogs[-1].mainshock = catalogs[-1].getMainEvent()
	#
	# some stuff we know:
	catalogs[0].mainshock[3] = targmag		# or 9.1... but anss sometimes reports this differently in think.
	#
	#
	fignum0=11
	current_fig = plt.figure(fignum0)
	#
	for i_catalog, catalog in enumerate(catalogs):
		#
		#mev = catalog.getMainEvent()
		mev = catalog.mainshock		# which we've added in this script. otherwise, use catalog.getMainEvent()
		print "mev: ", mev
		interval_len = rfp.winlen(m=mev[3], mc=mc, mt=mt, doInt=True)
		avlen = max(1,int(interval_len/10))
		#
		current_fig = plt.figure()
		plt.clf()
		catalog.rbomoriQuadPlot(targmag=mev[3], mc=mc, weighted=False, fignum=current_fig.number)
		#
		current_fig = plt.figure()
		plt.clf()
		catalog.rbomoriQuadPlot(targmag=mev[3], mc=mc, weighted=True, fignum=current_fig.number)
		#
		# poisson ratio limits:
		these_axes = current_fig.axes
		if len(these_axes)>=3:
			ax_rb = these_axes[2]		
			# poisson ratio bit:
			# this doesn't vary much, so for now just do a regular gr guess:
			n_rb_poisson = 10**(mev[3]-mc-2.0)
			poisson_ratio_0 = list(ratio_sigma(n_rb_poisson, 1.0))
			poisson_ratio_1 = list(ratio_sigma(n_rb_poisson, .5))
			print 'sigma_ratios: ', poisson_ratio_0, poisson_ratio_1
			#plt.figure(fignum0-1)
			for y in poisson_ratio_0 + poisson_ratio_1:
				#print catalog.getcat(0)[0]
				ax_rb.plot([catalog.getcat(0)[0][0], catalog.getcat(0)[-1][0]], [y, y], 'm-', lw=2)
		#
		#rbratios = catalog.getNRBratios(intervals=None, winlen=interval_len, delta_t=1, reverse=False)
		#rb_values = map(operator.itemgetter(4), rbratios)
		#mean_rbs = [numpy.mean(rb_values[max(i-interval_len, 0):i+1]) for i, x in enumerate(rb_values)]
		#for i,x in enumerate(mean_rbs): rbratios[i]+=[x]
		#
		#return rbratios
		#
		#fits = get_ratio_fits(ratios=rbratios, Nfits=[interval_len])
		#plotses = plot_ratio_fits(fits, new_figs=False, fignum0=fignum0)
		#fignum0+=4
	#
	#return (rbratios, rbratio_fits, rbf_plot)
	return None
	#
#	
def get_rbratios_tohoku(mc=5.0, todt=mhp.dtm.datetime.now(mhp.pytz.timezone('UTC')), winlen=None, avlen=None, targmag=9.0, ndithers=10, nContours=nContours1, bigmag=7.0, lons=[135., 148.5], lats=[30., 45.25], refreshcat=True, dt0=mhp.dtm.datetime(1990,1,1, tzinfo=mhp.pytz.timezone('UTC')), catname='cats/tohoku_rfits.cat', mt=7.6):
	#
	# fetch and return a standard rb-ratio set.
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
			if type(X_vals[0])==type(dtm.datetime.now()):
				these_x = map(mpd.date2num, X_vals[(i-n):i])
			else:
				these_x = X_vals[(i-n):i]
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
def get_averaging_set(ratios=None, ave_lens=[5, 500], x_col=1, ratio_col=4, fout='ratio_fit_data/aveset.csv'):
	# get plot of fit quality as a function of ave-len
	if ratios==None:
		ratios = get_rbratios_tohoku()
		# return None
	#
	if ave_lens[0]>ave_lens[1]: ave_lens.reverse()
	#
	if fout!=None:
		f=open(fout, 'w')
		f.write('#rb-ration fit sets for... some data set. ave_lens: %s, x_col=%d, ratio_col=%d\n' % (str(ave_lens), x_col, ratio_col))
		f.write('#ave_len\tmean_var\tAIC\ttotal_var\ts/n\tNdoF\n')
		f.close()
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
			chi_sqrs[-1]+=[ [lf.meanVar(), lf.AIC, lf.totalVar, lf.totalVar/these_r[-1], lf.Ndof] ]
		print "fit(%d)" % (ave_len)
		#
		my_arrays = zip(*chi_sqrs[-1])
		total_chi_sqrs += [ [ave_len, scipy.mean(my_arrays[0]), scipy.mean(my_arrays[1]), scipy.mean(my_arrays[2]), len(chi_sqrs[-1])] ]
		if fout!=None:
			f=open(fout, 'a')
			out_string = ''
			for col in total_chi_sqrs[-1]:
				out_string += '%s\t' % col
			f.write('%s\n' % out_string[:-1])
			f.close()
		#
	#
	return [total_chi_sqrs, chi_sqrs]
#
def get_averaging_set2(ratios=None, ave_lens=[5, 100], x_col=1, ratio_col=4):
	# get plot of fit quality as a function of ave-len
	if ratios==None:
		ratios = get_rbratios_tohoku()
		# return None
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
			chi_sqrs[-1]+=[ [lf.meanVar(), lf.AIC, lf.totalVar, lf.totalVar/these_r[-1], lf.Ndof] ]
		print "fit(%d)" % (ave_len)
		#
		my_arrays = zip(*chi_sqrs[-1])
		total_chi_sqrs += [ [ave_len, scipy.mean(my_arrays[0]), scipy.mean(my_arrays[1]), scipy.mean(my_arrays[2]), len(chi_sqrs[-1])] ]
		#
	#
	return [total_chi_sqrs, chi_sqrs]
#
def plot_randomized_total_chi_sqrs(rbratios=None, ave_lens=[5, 100], x_col=0, ratio_col=4, fignum0=1, fout='ratio_fit_data/rand_ave_set.csv'):
	#
	if rbratios==None:
		rbratios = get_rbratios_tohoku()
	#
	# now, randomize the ratio_col (we can either shuffle the existing values or just replace them with random numbers)
	# let's shuffle...
	x_temp = map(operator.itemgetter(ratio_col), rbratios)
	random.shuffle(x_temp)
	for i, x in enumerate(x_temp):
		rbratios[i][ratio_col]=x
	#
	avset_rand = get_averaging_set(ratios=rbratios, ave_lens=[5, 500], x_col=x_col, ratio_col=ratio_col, fout=fout)
	plot_return = plot_total_chi_sqrs(totals_list=avset_rand[0], ave_lens=ave_lens, x_col=x_col, ratio_col=ratio_col, fignum=fignum0, fout=None)
	#
	return avset_rand
#
def plot_total_chi_sqrs(totals_list=None, ave_lens=[5,100], x_col=1, ratio_col=4, fignum=0, fout=None):
	#
	# totals_list should come from the get_averaging_set[0]
	# run this with defaults.
	#
	if totals_list == None: totals_list = get_averaging_set(ratios=None, ave_lens=ave_lens, x_col=x_col, ratio_col=ratio_col, fout=fout)[0]
	#
	my_arrays = zip(*totals_list)
	sig_to_noise = [my_arrays[3][i]*math.sqrt(my_arrays[1][i]) for i in xrange(len(my_arrays[0]))]
	#
	plt.figure(fignum)
	plt.clf()
	ax1=plt.gca()
	
	ax1.plot(my_arrays[0], my_arrays[2], '-', label='aic vals')
	ax1.plot(my_arrays[0], my_arrays[3], '-', label='lacunarity')
	ax1.plot(my_arrays[0], sig_to_noise, '-', label = 's/n')
	plt.legend(loc=0, numpoints=1)
	#
	ax2=ax1.twinx()
	#
	ax2.plot(my_arrays[0], my_arrays[1], '.-', label='chi-squares')
	plt.legend(loc=0, numpoints=1)
	#
	return totals_list

def plot_total_chi_sqrs_files(f1='ratio_fit_data/tohoku_total_chi_sqrs.csv', f2='ratio_fit_data/rand_ave_tohoku.csv', fignum=0):
	#
	plt.figure(fignum)
	plt.clf()
	plot_func = plt.semilogy
	#plot_func = plt.plot
	#
	cols = [0, 1, 2]
	datas1 = []
	if f1!=None:
		fin = open(f1, 'r')
		for rw in fin:
			if rw[0]=='#': continue
			rws = rw.split()
			datas1 += [map(float, rws[0:3])]
		fin.close()
		#
		zdatas = zip(*datas1)
		plot_func(zdatas[0], zdatas[1], '-', label='f1')

	if f2!=None:
		datas2=[]
		fin=open(f2)
		for rw in fin:
			if rw[0]=='#': continue
			rws=rw.split()
			datas2 += [map(float, rws[0:3])]
		fin.close()
		zdatas = zip(*datas2)
		plot_func(zdatas[0], zdatas[1], '-', label='f2')
	plt.legend(numpoints=1, loc=0)
	#
	if f2!=None and f1!=None:
		delta_y = [abs(datas1[i][1]-datas2[i][1]) for i in range(len(datas1))]
		fact_y  = [(datas2[i][1]/datas1[i][1]) for i in range(len(datas1))]
		#
		#plot_func(zdatas[0], delta_y, '--')
		plot_func(zdatas[0], fact_y, '--')
	#
#
def plot_ratio_fits(fit_dict, new_figs=False, fignum0=2):
	#
	# plot data from get_ratio_fits()
	#
	#return fit_dict
	#
	fnum=fignum0
	plt.figure(fnum)
	plt.clf()
	#
	for fit_set in fit_dict.iterkeys():
		this_set=fit_dict[fit_set]
		#
		X   = this_set['X']
		y_r = this_set['means']
		y_a = map(operator.itemgetter(0), this_set['fit_prams'])
		y_b = map(operator.itemgetter(1), this_set['fit_prams'])
		y_chi = map(operator.itemgetter(2), this_set['fit_prams'])
		#
		plt.semilogy(X[-len(y_r):], y_r, '-', label='ratios')
		plt.semilogy(X[-len(y_chi):], y_chi, '--', label='chisqr: %.4f/%.4f/%.4f' % (sum(y_chi), scipy.mean(y_chi), len(y_chi)))
		plt.semilogy([X[-len(y_chi)], X[-1]], [1., 1.], 'k-', lw=2, zorder=1)
		plt.legend(loc=0, numpoints=1)
		plt.title('raw fits')
		#
		fnum+=1
		plt.figure(fnum)
		plt.clf()
		plt.title('raw lacunarity')
		#y_lac = numpy.array(y_r[-len(y_chi):])/numpy.array(y_chi)
		y_lac = [math.log10(y_r[-len(y_chi):][i])/ychi for i, ychi in enumerate(y_chi)]
		plt.plot(X[-len(y_chi):], y_lac, '-', lw=2, label='r/chi')
		plt.plot([X[-len(y_chi)], X[-1]], [0., 0.], 'k-', lw=2, zorder=1)
		plt.legend(loc=0)
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
		plt.semilogy(X[-len(y_r):], y_r, '-', label='ratios')
		plt.semilogy(X[-len(y_chi):], y_chi, '--', label='chisqr: %.4f/%.4f/%.4f' % (sum(y_chi), scipy.mean(y_chi), len(y_chi)))
		plt.semilogy([X[-len(y_chi)], X[-1]], [1., 1.], 'k-', lw=2, zorder=1)
		plt.legend(loc=0, numpoints=1)
		plt.title('mean-fits')
		#
		#
		fnum+=1
		plt.figure(fnum)
		plt.clf()
		plt.title('mean lacunarity')
		#y_lac = numpy.array(y_r[-len(y_chi):])/numpy.array(y_chi)
		y_lac = [math.log10(y_r[-len(y_chi):][i])/ychi for i, ychi in enumerate(y_chi)]
		plt.semilogy(X[-len(y_chi):], y_lac, '-', lw=2, label='r/chi')
		plt.semilogy([X[-len(y_chi)], X[-1]], [1., 1.], 'k-', lw=2, zorder=1)
		plt.legend(loc=0)
		#
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
		#
		# poisson ratio bit:
		rb_len=fit_set*10
		poisson_ratio_0 = list(ratio_sigma(rb_len, 1.0))
		poisson_ratio_1 = list(ratio_sigma(rb_len, .5))
		print 'sigma_ratios: ', poisson_ratio_0, poisson_ratio_1
		plt.figure(fignum0-1)
		for y in poisson_ratio_0: # + poisson_ratio_1:
			#y_prime = 10.**(math.log10(y)/math.sqrt(float(fit_set)))
			y_prime = y**(1./math.sqrt(float(fit_set)))		# note fit_set ~ .1*rb_len
			plt.plot([X[0], X[-1]], [y_prime, y_prime], 'm-', lw=2)
		#
	return None
	
def ratio_sigma(N, sigma_factor=1.0):
	r = (numpy.log(N) + sigma_factor*numpy.sqrt(numpy.log(N)))/(numpy.log(N) - sigma_factor*numpy.sqrt(numpy.log(N)))
	return (r, 1./r)	
	
	
	
	
	
	
	
	
	
	
	

