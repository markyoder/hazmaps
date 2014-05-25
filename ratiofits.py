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

nContours1=15


def get_cat_tohoku(mc=5.0, todt=mhp.dtm.datetime.now(mhp.pytz.timezone('UTC')), winlen=None, targmag=9.0, ndithers=10, nContours=nContours1, bigmag=7.0, lons=[135., 148.5], lats=[30., 45.25], refreshcat=True, dt0=mhp.dtm.datetime(1990,1,1, tzinfo=mhp.pytz.timezone('UTC')), catname='cats/tohoku_rfits.cat', mt=7.6):
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
	tohokuEvent=c1.getMainEvent()
	print "n_seq: %d" % N_sequence
	tohokuLatLon = tohokuEvent[1], tohokuEvent[2]
	Lr = 10.0**(.5*targmag - 1.76)
	lfactor=.5
	#
	c1.subcats+=[['r1', rfp.circularcat(c1.getcat(0), latlon=[tohokuEvent[1], tohokuEvent[2]], Rkm=Lr*lfactor)]]
	#
	#rbratios = c1.plotIntervalRatiosAx(winlen=N_sequence, cat=c1.getcat(1), hitThreshold=1.0, bigmag=9.0, thisAx=myaxes[i], ratios=None, delta_t=1, avlen=rbavelen, mainEv=mainshock, logZ=None, rbLegLoc=0, reverse=False)
	rbratios = c1.plotIntervalRatiosAx(winlen=N_sequence, cat=c1.getcat(1), hitThreshold=1.0, bigmag=9.0, thisAx=None, ratios=None, delta_t=1, avlen=1, mainEv=tohokuEvent, logZ=None, rbLegLoc=0, reverse=False)
	
	return rbratios

def get_ratio_fits(ratios=None, Nfits=[1,5], ratiocol=4):
	# ratios: return from rbratios (or any time series)
	# Nfits=[] list of Ns over which to do fitting
	# ratiocol: data column
	#
	# fits will include: mean value, linear fits {chi-square, a, b, sig_a, sig_b? (maybe later}
	# append these vals in a dict. to the end of each row: [{}, {}, {}..]
	#
	r_vals = map(operator.itemgetter(ratiocol), ratios)
	for n in Nfits:
		for i in xrange(n, len(ratios)):
			meanval = numpy.mean(rvals[i-n:n])	# or, it's faster to do a moving mean; mean+=(x_i-x_{i-n})
			# now, get a line-fit for this segment...
	#
	

