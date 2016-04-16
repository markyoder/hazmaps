import math
import scipy
#
import pylab as plt
import matplotlib
import matplotlib.mpl as mpl
import numpy
from PIL import Image as ipp
#import numpy.fft as nft
import scipy.optimize as spo
#from matplotlib import pyplot as plt
#from matplotlib import rc
from matplotlib.patches import Ellipse
#
import string
import sys
#from matplotlib import *
#from pylab import *
import os
import random
import time
#
# gamma function lives here:
#import scipy.special
from scipy.special import gamma
#from scipy.optimize import leastsq
from matplotlib import axis as aa
#
from threading import Thread
#
import datetime as dtm
import pytz
import calendar
import operator
import urllib

import eqcatalog as yp 
import rbIntervals as rbi

#
# python3 vs python2 issues:
# a bit of version 3-2 compatibility:
if sys.version_info.major>=3:
	xrange=range

class hazCrack(object):
	# record-breaking hazard map on crack.
	
	hazmaps=[]
	plotstyle=0		#{0: little squares, 1: little squares mean-val, 2: mean-val contours, 3: overlapping big squares, 4: mean-val big squares}
	gridsize=1.0	# sometimes called L; it's the length of the grid size.
	ndithers=3	# so the total number of elements will be ndithers^2 (aka, 3-> dither into thirds on each side).
	dL=None	# will be L/ndithers
	hmdate=None
	#
	#
	def __init__(self, catalog=None, catnum=0, mc=2.5, gridsize=.5, ndithers=3, winlen=128, avlen=1, bigmag=5.0, fignum=0, logZ=1.0, mapres='i'):
		self.dL=gridsize/float(ndithers)
		self.hazmaps=[]
		self.mc=mc
		self.gridsize=gridsize
		self.ndithers=ndithers
		self.winlen=winlen
		self.avlen=avlen
		self.bigmag=bigmag
		self.fignum=fignum
		self.mapres=mapres
		if logZ==None: logZ=math.log10(winlen)
		self.logZ=logZ		# for log-normalizing rb-ratios 
		#self.logZ=1.0
		#self.forecastDate=dtm.datetime.now(pytz.timezone('UTC'))	# or we might change this.
		self.forecastDate=dtm.datetime(2011, 3,5, 0, 0, 0, 0, pytz.timezone('UTC'))
		self.hmdate=dtm.datetime.now(pytz.timezone('UTC'))
		self.rb=rbi.intervalRecordBreaker(None)
		#
		for i in xrange(ndithers*ndithers):
			xphase=(i%ndithers)*self.dL				# lon	phase
			yphase=(int(i/ndithers))*self.dL		# lat pahse
			# in the end, it is going to be easier to rewrite the hazmap object. make a lean class excluding the graphical parts - just
			# an array of subcats. we'll load a single own rb object to this class or rewrite the RB counting script.
			# catalog=None, catnum=0, mc=2.5, gridsize=.5, phi=[0,0], winlen=128, avlen=1, bigmag=5.0
			thishm=rbHazMapLean(catalog=catalog, catnum=catnum, mc=mc, gridsize=gridsize, phi=[yphase, xphase], winlen=winlen, avlen=avlen, bigmag=bigmag, logZ=logZ) # but we'll need to suppress fig-plotting and i don't know if the phases are switched.
			a=thishm.hazMapTo(self.forecastDate)
			self.hazmaps+=[thishm]
	
	def getcat(self, catindex=0):
		if len(self.hazmaps)<1:
			return None
		if catindex>(len(self.hazmaps)-1):
			catindex=len(self.hazmaps)-1
		#
		return self.hazmaps[catindex].catalog
	
	def getNlon(self):
		#nlon=0
		#for hm in self.hazmaps:
		#	nlon+=hm.Nlon
		return self.hazmaps[0].Nlon*self.ndithers
		#return nlon
		
	def getNlat(self):
		#nlat=0
		#for hm in self.hazmaps:
		#	nlat+=hm.Nlat
		return self.hazmaps[0].Nlat*self.ndithers
		#return nlat
	
	def boxesTo(self, hmdate=dtm.datetime.now(pytz.timezone('UTC')), fignum=0):
		self.hazMapTo(hmdate)
		return self.simpleBoxes(fignum)
	
	def hazMapTo(self, hmdate=dtm.datetime.now(pytz.timezone('UTC'))):
		self.hmdate=hmdate
		for i in xrange(len(self.hazmaps)):
			self.hazmaps[i].hazMapTo(hmdate)
		
	def contourData(self, dobinary=0):
		# for our first trick, just combine all the hazard maps into one big map and contour it.
		# later, we'll do some averaging schemes.
		#
		# first, get one big array of X,Y,Z:
		#cX=[]
		#cY=[]
		#cZ=[]
		cdata=[]
		#for i in xrange(len(self.hazmaps[0].rbRatios)):
		#	for hm in self.hazmaps:
		for hm in self.hazmaps:
			for i in xrange(len(hm.rbRatios)):
				#cX+=[hm.getLon(i) + (self.gridsize-self.dL)/2.0]
				#cY+=[hm.getLat(i) + (self.gridsize-self.dL)/2.0]
				#cZ+=[hm.rbRatios[i]]
				#
				thisval=hm.rbRatios[i]
				#if thisval==None: thisval=0
				if thisval!=None:
					thisval=-thisval		# because we've defined r<1 to be high risk. the easiest way to color these "hot" colors is to reverse them.
					if dobinary:
						# return binary (or arguably tri-mary, {-1, None, 1} (noting at this point that the parity has been reversed)
						if thisval>=0: thisval=1
						if thisval<0: thisval=-1
				# so these should be bin-centered:
				cdata+=[[hm.getLon(i) + (self.gridsize-self.dL)/2.0, hm.getLat(i) + (self.gridsize-self.dL)/2.0, thisval]]
		#
		# so now we have a [[x, y, rb]] array from all the hazard maps.
		# sort them by x,y (row vectors if you will). [x0, y0], [x1, y0], [x2, y0]..., [x0, y1], [x1, y1], [x2, y1], ...]
		# sort first by the secondary key (x), then on primary (y):
		cdata.sort(key=operator.itemgetter(0))
		cdata.sort(key=operator.itemgetter(1))
		#
		# give all the clusters a 0-value border (for prettiness).
		#width=self.hazmaps[0].Nlon*self.ndithers
		width=self.getNlon()
		maxi=len(cdata)	# actually, max i + 1, so tread carefully.
		#
		# pad the clusters to make better contours:
		# spin through cdata. if for every element, if an adjacent cell is None and not over an edge, make it zero. note this is made easy because
		# the lat/lon are provided...
		for i in xrange(maxi):
			if cdata[i][2]==0. or cdata[i][2]==None: continue	#(0, None) because, if we start placing 0 values, we fill the whole space with 0's.
			#thislat=cdata[1]
			#thislon=cdata[0]
			# on the left border?
			if i%width>0:
				if cdata[i-1][2]==None: cdata[i-1][2]=0.
			# right:
			if i%width<(width-1):
				if cdata[i+1][2]==None: cdata[i+1][2]=0.
			# top:
			#if i/width!=0:
			if (i-width)>=0:
				if cdata[i-width][2]==None: cdata[i-width][2]=0.
			if (i+width)<maxi:
				if cdata[i+width][2]==None: cdata[i+width][2]=0.
			#
			
		#
		# though we might eventually return a later form of this data (as per see contourPlot()...
		return cdata

	def contourArrays(self):
		#
		cdata=self.contourData()
		# now, borrowing from ffgridplot.py contouring methods, we "shape" the array:
		# can we just contour() from here? no.
		#
		Xs=map(operator.itemgetter(0), cdata)
		Ys=map(operator.itemgetter(1), cdata)
		Xvals=list(set(Xs))
		Yvals=list(set(Ys))
		Xvals.sort()
		Yvals.sort()
		#
		#Xax=list(set(map(operator.itemgetter(0), cdata))).sort()
		#Yax=list(set(map(operator.itemgetter(1), cdata))).sort()
		#
		X,Y=numpy.meshgrid(Xvals, Yvals)
		#z=map(operator.itemgetter(2), cdata)
		Z=map(operator.itemgetter(2), cdata)
		Z2d=scipy.array(Z)
		Z2d.shape=(len(Yvals), len(Xvals))
		#
		# for now, this is a fully class-contained function, so we can set some class-scope
		# variables.
		self.X_i = X
		self.Y_i = Y
		self.Z2d = Z2d
		#
		return [X,Y,Z2d]
		
	def contourPlot(self, fignum=1):
		
		cary=self.contourArrays()
		X = cary[0]
		Y = cary[1]
		Z2d = cary[2]
		
		plt.ion()
		plt.figure(fignum+1)
		plt.clf()
		plt.figure(fignum)
		plt.clf()
		#plt.contourf(map(operator.itemgetter(0), cdata), map(operator.itemgetter(1), cdata), map(operator.itemgetter(2), cdata))
		#return X, Y, Z2d
		self.conts = plt.contourf(X,Y,Z2d, alpha=.3)
		plt.colorbar()
		plt.spectral()
		#return [cX, cY, cZ]
		#return cdata
		return None
	
	def plotEvent(self, coords=None, symbolstr='r*', msize=18, fignum=0, zorder=5):
		cm=self.hazmaps[0]
		if coords==None: coords=cm.getMainEvent()
		x,y=cm.catalog.catmap(coords[0], coords[1])
		plt.figure(fignum)
		plt.plot([x], [y], symbolstr, ms=msize, zorder=zorder)
		
		return None
	
	def plotEvents(self, events=None, symbolstr=',', msize=1, fignum=0, alpha=1.0, zorder=5):
		cm=self.hazmaps[0]
		print("doing events. %d, %d" % (len(events), fignum))
		#
		#fignum=1
		# assume events come in 'catalog' form [dtm, lat, lon, mag, ...]
		lats=map(operator.itemgetter(1), events)
		lons=map(operator.itemgetter(2), events)
		X,Y=cm.catalog.catmap(lons, lats)
		plt.figure(fignum)
		plt.plot(X,Y, symbolstr, ms=msize, alpha=alpha, zorder=zorder)
		
		return None
	
	def contourMap(self, fignum=1, nconts=12, mapres=None, maptransform=True, alpha=.85):
		# maptransform: transform to matplotlib.basemap coords. for KML export, set this to false (coords in LL).
		# (though we might not use this feature; we might do it another way...)
		if mapres==None: mapres=self.mapres
		cdata=self.contourData()
		if nconts==None: nconts=12	# but we can also pass a list of contour values (to standardize output)[c1, c2, c3...]
		#
		# set up the map. we'll need it to process the contour data.
		#
		catdate=self.hazmaps[0].catalog.cat[-1][0]	# most recent date in catalog.
		hmdatestr='%d-%d-%d, %d:%d:%d (UTC)' % (self.hmdate.year, self.hmdate.month, self.hmdate.day, self.hmdate.hour, self.hmdate.minute, self.hmdate.second)
		#
		plt.ion()
		fig1=plt.figure(fignum+1)
		plt.clf()
		#plt.title("Record-breaking contour hazard map\ndtm=%s, $d\\lambda=%.4f$, nDith=%d, mc=%.2f, rblen=%d, avlen=%d\n\n" % (str(self.hmdate), self.gridsize, self.ndithers, self.mc, self.winlen, self.avlen))
		plt.title("Record-breaking contour hazard map\ndtm=%s, $d\\lambda=%.4f$, nDith=%d, mc=%.2f, rblen=%d, avlen=%d\n\n" % (hmdatestr, self.gridsize, self.ndithers, self.mc, self.winlen, self.avlen))
		ax1  = fig1.add_subplot(111)
		#
		fig=plt.figure(fignum)
		fig.clf()
		ax  = fig.add_subplot(111)
		#plt.title('(diagnostic) Record-breaking intervals Hazard Map')
		plt.title("Record-breaking contour haz-map\ndtm=%s, $d\\lambda=%.4f$, nDith=%d, mc=%.2f, rblen=%d\n\n" % (hmdatestr, self.gridsize, self.ndithers, self.mc, self.winlen))
		#
		llrange=self.hazmaps[0].catalog.getLatLonRange(self.hazmaps[0].thiscat)
		self.llrange = llrange
		self.deltaLat=abs(float(self.llrange[1][0])-self.llrange[0][0])
		self.deltaLon=abs(float(self.llrange[1][1])-self.llrange[0][1])
		Nlat=1+int(self.deltaLat/self.gridsize)
		Nlon=1+int(self.deltaLon/self.gridsize)	#number of cells in row/col
		#
		lat0=self.llrange[0][0]
		lon0=self.llrange[0][1]
		#
		#thismev=self.hazmaps[0].catalog.getMainEvent(self.hazmaps[0].catalog.getcat(0))
		#epicen=[thismev[2], thismev[1]]
		try:
			epicen=self.epicen
		except:
			c=self.hazmaps[0].catalog.getcat(0)
			i=0
			m0=c[i][3]
			i0=0
			while i<len(c):
				if c[i][3]>m0 and c[i][0]>self.hmdate: 
					m0=c[i][3]
					i0=i
				i+=1
			self.epicen=[c[i0][2], c[i0][1]]
		
		plt.figure(fignum)
		plottingCat=yp.eqcatalog(self.hazmaps[0].catalog.getcat(0))
		# now, make two sub-cats, one before, one after the haz-date:
		# addTimeRangeCat(self, subcatname='dtSubcat', fullcat=None, dtFrom=None, dtTo=None)
		plottingCat.addTimeRangeCat(subcatname='before', fullcat=None, dtFrom=self.hmdate-dtm.timedelta(days=1000), dtTo=self.hmdate)
		plottingCat.addTimeRangeCat(subcatname='after', fullcat=None, dtFrom=self.hmdate, dtTo=self.hmdate+dtm.timedelta(days=1000))
		plotcat=[[plottingCat.subcats[0][0], plottingCat.subcats[0][1]], [plottingCat.subcats[1][0], plottingCat.subcats[1][1][-1000:]]]
		#cm=self.hazmaps[0].catalog.plotCatMap(self.hazmaps[0].catalog.getcat(0), True, False, None, epicen, 'best', True, 'g,', ax, fignum)
		self.hazmaps[0].catalog.mapres=mapres
		cm=self.hazmaps[0].catalog.plotCatMap(self.hazmaps[0].catalog.getcat(0), True, False, None, None, 'best', True, 'g,', ax, fignum, plotevents=1)
		
		#cm1=plottingCat.plotCatsMap(catalogses=[plottingCat.subcats[0]], legendLoc='best', doCLF=True, fignum=fignum+1, bigmag=6.0)
		#cm1=self.hazmaps[0].catalog.plotCatsMap(catalogses=plottingCat.subcats, legendLoc='best', doCLF=False, fignum=fignum+1, bigmag=self.bigmag)
		#cm1=self.hazmaps[0].catalog.plotCatsMap(catalogses=[plottingCat.subcats[0], plottingCat.subcats[1]], legendLoc='best', doCLF=False, fignum=fignum+1, bigmag=self.bigmag)
		#cm1=self.hazmaps[0].catalog.plotCatsMap(catalogses=plotcat, legendLoc='best', doCLF=False, fignum=fignum+1, bigmag=self.bigmag)
		#### (this (above) produces a skewed catalog-map. instead, plot a map with no events and add the events separately.
		cm1=self.hazmaps[0].catalog.plotCatMap(self.hazmaps[0].catalog.getcat(0), True, False, None, None, 'best', True, 'g,', ax1, fignum+1, plotevents=0)
		#cm=self.hazmaps[0].catalog.plotCatMap(self.hazmaps[0].catalog.getcat(0), True, False, None, epicen, 'best', True, 'g,', ax, fignum)
		####
		plt.figure(fignum)
		#
		Xs=map(operator.itemgetter(0), cdata)
		Ys=map(operator.itemgetter(1), cdata)
		Xvals=list(set(Xs))
		Yvals=list(set(Ys))
		Xvals.sort()
		Yvals.sort()
		#
		X,Y=numpy.meshgrid(Xvals, Yvals)
		X1, Y1=numpy.meshgrid(Xvals, Yvals)
		# and now, we have to mash these into map coordinates:
		for i in xrange(len(Y)):
			for j in xrange(len(Y[0])):
				myx, myy = X[i][j], Y[i][j]
				myxm, myym = cm(myx, myy)		# map coords...
				#
				if maptransform==True:
					myxm1, myym1 = cm1(myx, myy)
				else:
					myxm1, myym1 = myx, myy
				#
				X[i][j]=myxm
				Y[i][j]=myym
				#
				X1[i][j]=myxm1
				Y1[i][j]=myym1
				#
				#X[i][j],Y[i][j] =cm(myx, myy)
		#
		Z=map(operator.itemgetter(2), cdata)
		Z2d=scipy.array(Z)
		Z2d.shape=(len(Yvals), len(Xvals))
		#
		# and pyplot has a weird, or stupid, feathre such that if we truncate the upper/lower contour range, it just leaves
		# a hole (instead of inferring all z>zmax -> zmax). so, we "clip".
		myZ = numpy.array(Z2d, copy=True)
		#if nconts.insert:
		try:
			nconts.insert
			# nconts has an "insert" function or property anyway (more properly see if
			# type(a.insert).__name__== 'builtin_function_or_method'
			#zmin=nconts[0]
			#zmax=nconts[-1]
			#myZ = Z2d.clip(nconts[0], nconts[-1])
			# the problem with .clip() is that it replaces all the None values with zmin, so jus do it
			# old-school:
			for iy in xrange(len(myZ[0])):
				if myZ[ix][iy]==None: continue
				if myZ[ix][iy]>nconts[-1]: myZ[ix][iy] = nconts[-1]
				if myZ[ix][iy]<nconts[0]: myZ[ix][iy] = nconts[0]
		except:
			# jut move on...
			dummyvar=None
			
			
		#
		plt.figure(fignum)
		#plt.clf()
		#self.conts = plt.contourf(X,Y,Z2d, nconts, alpha=alpha, zorder=10)
		self.conts = plt.contourf(X,Y,myZ, nconts, alpha=alpha, zorder=10)
		plt.colorbar()
		plt.spectral()
		
		plt.figure(fignum+1)
		self.conts2 = plt.contourf(X1,Y1,myZ, nconts, alpha=alpha, zorder=10)
		plt.colorbar()
		plt.spectral()
		
		return self.conts2
	
	def makeColorbar(self, cset=None, colorbarname=None, reffMag=5.0, contfact=5.0):
		#
		# reffMag:
		if reffMag==None:
			#reffMag=self.reffMag
			reffMag=5.0
		#
		if colorbarname==None: colorbarname='scale.png'
		if cset==None: cset=self.conts
		cs=cset.collections
		#resolution = 5
		#LevelsNumber = 5 * resolution
		warnings = ['No','Low','Guarded','Elevated','High','Severe']
		#
		#resolution = int(len(cs)/self.contfact)
		resolution = int(len(cs)/contfact)
		#
		startindex=0
		#
		fig = plt.figure(figsize=(1.5,3))
		axe = fig.add_axes([0.05, 0.05, 0.2, 0.9])
		#
		fig.figurePatch.set_alpha(0)
		#
		matplotlib.rc('font', family='monospace', weight='black')
		#
		cmap = mpl.cm.spectral
		#norm = mpl.colors.Normalize(vmin=Z_i.min(), vmax=Z_i.max())
		#ratefactor=10**(self.mc-self.reffMag)
		ratefactorExp=self.mc-reffMag	# so the raw rate probably should be with mc_catalog = mc_etas.
													# but, we'll want to report some rate of m>reffMag
		
		#
		#tics = [0,10,20,40,80]
		#tics = [round(self.Z2d.min(),3), round(self.Z2d.max(),3)]
		#timefactExp = math.log10(year2secs)
		# min, max functions puking on None types, so do it manually:
		zmin=None
		zmax=None
		for i in xrange(len(self.Z2d)):
			for ii in xrange(len(self.Z2d[0])):
				z0=self.Z2d[i][ii]
				if z0==None: continue
				if zmin==None: zmin=z0
				if zmax==None: zmax=z0
				#
				if z0>zmax: zmax=z0
				if z0<zmin: zmin=z0
				#
		norm = mpl.colors.Normalize(vmin=zmin, vmax=zmax)
		#
		t1='%.2f' % (zmin)
		t2='%.2f' % (zmax)
		#
		# poisson noise thresholds:
		poissonSigma = (.5*math.log10(self.winlen))/self.logZ
		tsn = '-%.2f' % poissonSigma
		tsp = '%.2f' % poissonSigma
		#print("max set: ", self.Z2d.max(), timefactExp, ratefactorExp)
		tics = [zmin, -poissonSigma, 0.0, poissonSigma, zmax]
		tlbls = [t1, tsn, 0.0, tsp, t2]
		# clean up labels/tics in case max/min values are near poisson threshold.
		if abs(poissonSigma - zmax)<.1:
			tics.pop(-1)
			tlbls.pop(-1)
		if abs(-poissonSigma - zmin)<.1:
			tics.pop(0)
			tlbls.pop(0)
		#
		#cb1 = mpl.colorbar.ColorbarBase(axe, norm=norm, ticks=tics, format="%g%%", orientation='vertical')
		cb1 = mpl.colorbar.ColorbarBase(axe, norm=norm, ticks=tics, format="%.2f", orientation='vertical')
		cb1.set_ticklabels(tlbls, update_ticks=True)
		cb1.set_label('Local normalized RB ratio\n$r=\\frac{\\log (n_{rb-short}/n_{rb-long})}{\\log(N)}$')
		#cb1.set_label('ETAS rate')
		plt.savefig(colorbarname)	# vertical
		#
		#suffix=''
		i=len(colorbarname)-1
		#while colorbarname[i]!='.':
		#	#suffix.insert(0, colorbarcopy[-i])
		#	# this will get us the last '.'
		#	i-=1
		i=colorbarname.rfind('.')
		hname=colorbarname[:i] + '-h' + colorbarname[i:]
		#
		im = ipp.open(colorbarname, 'r')
		im.load()							# just to be sure. open() only creates a pointer; data are not loaded until analyzed.
		im=im.rotate(-90)
		#im.show()
		im.save(hname)
		#
		# depricated method:
		#os.system('convert -rotate 90 %s %s' % (colorbarname, hname))
		#os.system('cp %s %s' % (colorbarname, hname))
		#
		return [colorbarname, hname]	
			
	def simplePlot(self, fignum=0, thresh=0.):
		# do a simple matplotlib plot of hazmap. this is probably diagnostic, so it may just be points of some sort.
		# note: this plots the ll-corner of each box. for a science-interesting plot, we should center these points (aka, [x+L/2, y+L/2]
		#
		if thresh==None or type(thresh).__name__ not in ('float', 'int'): thresh=0.0
		plt.figure(fignum)
		plt.clf()
		plt.ion()
		#
		plt.title("diagnostic plot of dithered rb-hazmaps")
		hmindex=0
		for hm in self.hazmaps:
			llrange=hm.llrange
			lonsRed=[]
			latsRed=[]
			#rbvals=[]
			lonsBlue=[]
			latsBlue=[]
			#
			lonsCyan=[]
			latsCyan=[]	# for None , at least temporarily
			#
			for i in xrange(len(hm.rbRatios)):
				#
				#x=llrange[0][1]+self.gridsize*(i%hm.Nlon)			#lon
				#y=llrange[0][0]+self.gridsize*int(i/hm.Nlon)		# lat
				x=hm.getLon(i)
				y=hm.getLat(i)
				#if x<0: x+=180.	# this seems to be more trouble than it's worth...
				if hm.rbRatios[i]==None:
					continue
					#lonsCyan+=[x]
					#latsCyan+=[y]
				#
				if hm.rbRatios[i]<=thresh:
					lonsRed+=[x]
					latsRed+=[y]
				elif hm.rbRatios[i]>thresh:
					lonsBlue+=[x]
					latsBlue+=[y]
			#plt.plot(lonsRed, latsRed, 'r%s' % yp.pyicon(hmindex))
			#plt.plot(lonsBlue, latsBlue, 'b%s' % yp.pyicon(hmindex))
			plt.plot(lonsRed, latsRed, 'r.')
			plt.plot(lonsBlue, latsBlue, 'b.')
			
			#plt.plot(lonsCyan, latsCyan, '%s' % yp.pyicon(hmindex), color='y')
			hmindex+=1
	
	def simpleBoxes(self, fignum=0, thresh=0.0, plotevents=1, mapres=None):
		# do a simple matplotlib plot of hazmap. this is probably diagnostic, so it may just be points of some sort.
		# note: this plots the ll-corner of each box. for a science-interesting plot, we should center these points (aka, [x+L/2, y+L/2]
		# eventually, turn this into a map plot.
		#
		if mapres==None: mapres=self.mapres
		#if thresh==None or type(thresh).__name__ not in ('float', 'int'): thresh=0.0
		noisethresh=thresh
		if thresh==None:
			thresh=((self.winlen/2.)-.577 - math.log(self.winlen))/((self.winlen/2.)+.577 + math.log(self.winlen))
			thresh=abs(math.log10(thresh))
			# print("noisethresh: %f" % (noisethresh))
		#
		catdate=self.hazmaps[0].catalog.cat[-1][0]
		hmdatestr='%d-%d-%d, %d:%d:%d (UTC)' % (self.hmdate.year, self.hmdate.month, self.hmdate.day, self.hmdate.hour, self.hmdate.minute, self.hmdate.second)
		#
		plt.ion()
		self.fig=plt.figure(fignum)
		self.fig.clf()
		self.ax  = self.fig.add_subplot(111)
		#plt.title('(diagnostic) Record-breaking intervals Hazard Map')
		#plt.title("Record-breaking hazardmap\ndtm=%s, $d\\lambda=%.4f$, nDith=%d, mc=%.2f, rblen=%d, avlen=%d\n\n" % (str(self.hmdate), self.gridsize, self.ndithers, self.mc, self.winlen, self.avlen))
		plt.title("Record-breaking hazardmap\ndtm=%s, $d\\lambda=%.4f$, nDith=%d, mc=%.2f, rblen=%d, avlen=%d\n\n" % (hmdatestr, self.gridsize, self.ndithers, self.mc, self.winlen, self.avlen))
		#
		self.llrange=self.hazmaps[0].catalog.getLatLonRange(self.hazmaps[0].thiscat)
		self.deltaLat=abs(float(self.llrange[1][0])-self.llrange[0][0])
		self.deltaLon=abs(float(self.llrange[1][1])-self.llrange[0][1])
		self.Nlat=1+int(self.deltaLat/self.gridsize)
		self.Nlon=1+int(self.deltaLon/self.gridsize)	#number of cells in row/col
		#
		self.lat0=self.llrange[0][0]
		self.lon0=self.llrange[0][1]
		#
		#ilat=0
		#ilon=0	# index counters
		#self.cm=self.catalog.plotCatMap(self.thiscat, True, False, None, None, 'best', True, 'g,', None, self.fignum)
		#self.cm=self.catalog.plotCatMap(self.thiscat, True, False, None, None, 'best', True, 'g,', self.ax)
		#minicat=self.catalog.getcat(0)[0:10]
		#
		#thismev=self.hazmaps[0].catalog.getMainEvent(self.hazmaps[0].catalog.getcat(0))
		#epicen=[thismev[2], thismev[1]]
		#
		#self.cm=self.hazmaps[0].catalog.plotCatMap(self.hazmaps[0].catalog.getcat(0), True, False, None, self.epicen, 'best', True, 'g,', self.ax, fignum, plotevents)
		self.hazmaps[0].catalog.mapres=mapres
		self.cm=self.hazmaps[0].catalog.plotCatMap(self.hazmaps[0].catalog.getcat(0), True, False, None, self.epicen, 'best', True, 'g,', self.ax, fignum, plotevents)
		#
		plt.figure(fignum)
		hmindex=0
		for hm in self.hazmaps:
			llrange=hm.llrange
			lonsRed=[]
			latsRed=[]
			#rbvals=[]
			lonsBlue=[]
			latsBlue=[]
			#
			lonsCyan=[]
			latsCyan=[]	# for None , at least temporarily
			#
			for i in xrange(len(hm.rbRatios)):
				#
				#x=llrange[0][1]+self.gridsize*(i%hm.Nlon)			#lon
				#y=llrange[0][0]+self.gridsize*int(i/hm.Nlon)		# lat
				x=hm.getLon(i) + self.gridsize/2.0 - self.dL/2.0
				y=hm.getLat(i) + self.gridsize/2.0 - self.dL/2.0	 # to center the little boxes on their area of influence. note, this visually
																					 # excludes the left side of the cluster.
				# let's try always casting all lon values as positive:
				# if x<0: x+=180.	# basemap appears to not understand the int. dt. line.
																
				if hm.rbRatios[i]==None:
					continue
					#lonsCyan+=[x]
					#latsCyan+=[y]
				#
				#mapX, mapY=self.cm(map(operator.itemgetter(0), thissquare), map(operator.itemgetter(1), thissquare))
				#self.gridSquares+=self.ax.fill(mapX, mapY, alpha=.1, color='b', lw=2, zorder=2)
				thissq=hm.getSquare([x,y], float(self.gridsize)/self.ndithers)
				#thissq=hm.getSquare([x,y], float(self.gridsize))
				X=map(operator.itemgetter(0), thissq)
				Y=map(operator.itemgetter(1), thissq)
				X,Y=self.cm(X,Y)
				#thisalpha=abs(hm.rbRatios[i])
				if hm.rbRatios[i]<(-thresh):
					# (note: ratios are "loggoe": r[i]=log(nrb1/nrb2)
					# getSquare(self, xy=[], dx=.1):
					#plt.fill(X,Y, alpha=thisalpha, color='r', lw=2, zorder=9)
					plt.fill(X,Y, alpha=.55, color='r', lw=2, zorder=9)
					#lonsRed+=[x]
					#latsRed+=[y]
				elif hm.rbRatios[i]>=thresh:
				#elif hm.rbRatios[i]>thresh:
					#lonsBlue+=[x]
					#latsBlue+=[y]
					#plt.fill(X,Y, alpha=thisalpha/2., color='b', lw=2, zorder=8)
					plt.fill(X,Y, alpha=.2, color='b', lw=2, zorder=8)
				elif hm.rbRatios[i]>=-thresh and hm.rbRatios[i]<0.:
					plt.fill(X,Y, alpha=.4, color='y', lw=2, zorder=8)
				elif hm.rbRatios[i]>=0. and hm.rbRatios[i]<=thresh:
					#lonsBlue+=[x]
					#latsBlue+=[y]
					plt.fill(X,Y, alpha=.2, color='c', lw=2, zorder=8)
				
			#plt.plot(lonsRed, latsRed, 'r%s' % yp.pyicon(hmindex))
			#plt.plot(lonsBlue, latsBlue, 'b%s' % yp.pyicon(hmindex))
			
			#plt.plot(lonsRed, latsRed, 'r.')
			#plt.plot(lonsBlue, latsBlue, 'b.')
			
			#plt.plot(lonsCyan, latsCyan, '%s' % yp.pyicon(hmindex), color='y')
			hmindex+=1	
	
	def getClusters(self, beta=2):
		#cdata=self.contourData()
		clusts=self.findClusters(map(operator.itemgetter(2), self.contourData()), self.getNlon())
		'''
		clusts=[]
		for rw in clusts0:
			if len(rw[1])==0: continue
			clusts+=rw[1]
		'''
		return clusts
		#
	def getClustStats(self, clusts=None, beta=2):
		if clusts==None: clusts=self.getClusters()
		clustStats=[]
		for i in xrange(len(clusts)):
		#for rw in clusts:
			rw=clusts[i]
			#clustStats+=[[len(rw), self.guessMag(len(rw), beta=beta)]]
			clustStats+=[[i, len(rw), self.guessMag(len(rw), beta=beta)]]
		#clustStats.sort()	# i think this will, by default, sort by first element. if not, use a key, itemgetter(), etc.
		clustStats.sort(key=lambda rw: rw[1])
		return clustStats
	
	def plotClustRBratios(self, clustlist=None, rblen=None, minmag=1.5, fignum=1, avlen=1, justTop=False):
		# notes: the idea here is to see if a cluster is coherent, in other words, if two adjacen sites show r<1, is r<1 for both sites together?
		# a problem arises for small clusters in dithered maps. for example, if ndithers=3, each square represents only 1/9 of a full sequence.
		# we estimate the area and sequence length basically from n(rblen/9). in the simple example of 1 square, we end up trying to compute rb-statistics
		# as n=rblen/9, when in fact the statistics for that cell were computed over rblen intervals, so the result it toally meaningless. not until
		# N_clusts is large does this approximation work. we could try to re/de-construct the clusters to figure out the actual area covered and
		# actual number of events in the calculation, but as per shape-statistics, this will be nasty.
		#
		# given this problem of estimating cluster areas, maybe the thing to do is to look at the source haz-maps themselves. basically filter non-coherent
		# clusters in hazmaps[i], then plot only the coherent clusters. of course, the next problem is that our coherent cluster might be connected to other 
		# clusters by noisy pixels. maybe this problem will be mitigated in the dithering, but we might end up filtering out our target. finding an efficient
		# way to look for coherent clusters within larger clusters will be a big job, but it might be tractable. of course, we'll be back to the question of,
		# how do we define "coherent". all red? mostly red?
		#
		if clustlist==None: clustlist=self.findClusters()	# so do all of them?
		if type(clustlist[0]).__name__ == 'int':
			# we've got a list of indeces.
			tmpclist=self.findClusters()
			newlist=[]
			for i in clustlist:
				newlist+=[tmpclist[i]]
			clustlist=newlist
			tmpclist=None
			newlist=None
		#
		c1=yp.eqcatalog(self.catFromClusters(clustlist, justTop))
		if rblen==None:
			# optimal (??) rblen:
			nsquares=0
			for clst in clustlist:
				nsquares+=len(clst)
			rblen=2*int(self.winlen*nsquares/(self.ndithers**2))
			print("rblen = %d" % rblen)
		self.rb.fullCat=c1.getcat(0)
		self.rb.shockCat=c1.getcat(0)
		#
		self.rb.plotIntervalRatios(windowLen=rblen, minmag=minmag, fignum=fignum, avlen=avlen)
	
	def catFromClusters(self, clustlist=[0], justTop=False):
		# form a sub-catalog from one or more clusters.
		# justTop -> take only top self.winlen events from each cluster. this weights each cell equally, rather than focusing on the most active sites.
		#
		# clustlist is a collection (list) of clusters from findClusters(). each entry is a sub-cat in a hazmap.catalog object.
		#
		# since we've dithered, we have to choose a way to choose which events. aka, do we choose all events inside the dithered boundary as we've drawn it,
		# or all elements in the catalogs used to draw the cluster... or, a set of sub-cats for each dithered set.
		# for now, combine all the catalogs (they will overlap a lot) and then remove duplicates.
		# 
		# de-index the cluster-map (nLon, nLat are for for dithered set; ndithers is dither-length (ndithers=3 -> 3x3), / implies integer division).
		# see getMapIndex() for index to mapID mapping. as per dithering, we're arranging our maps like (where the mapID's below represent the indeces
		# of the maps as they were created. it would probably be a good idea to actually name the maps to make this process easier).
		# 0  1  2| 9 10 11 | 18 19 20
		# 3  4  5|12 13 14 | 21 22 23
		# 6  7  8|15 16 17 | 23 25 26
		# 27 o  o  *  o  o  * o o
		# o  o  o  o  o  o  o o o
		# o  o  o  o  o  o  o o o
		# *  o  o  *  o  o  * o o
		# 
		# ndithers=3
		# nLon=9
		# then, the mapID is mapIndex%(ndithers**2),
		#       the gridIndex (for that map) is mapID/(ndithers**2)
		#
		# so combine all the indeces from the entries in clustlist, make one big catalog from all the events; strip duplicates.
		newcat=[]
		i0=0
		if justTop==True: i0=-self.winlen
		for clust in clustlist:
			for elem in clust:
				#print("clustLen: %d" % len(clust[i0:]))
				(hazmapIndex, subcatIndex) = self.getMapIndex(elem)
				#newcat+=self.hazmaps[hazmapIndex].catalog.subcats[subcatIndex][1]	# because eacy point on the haz-map represents a stastic from a set of events in a subcat.
				newcat+=self.hazmaps[hazmapIndex].catalog.subcats[subcatIndex][1][i0:]	# because eacy point on the haz-map represents a stastic from a set of events in a subcat.
				#print(hazmapIndex, subcatIndex)
		#
		# now, strip duplicates:
		#rcat=list(set(newcat))
		#rcat=newcat[[0]]
		rcat=[]
		i=0
		while i<len(newcat):
			#if newcat[i]==rcat[-1]: continue
			if newcat[i] not in rcat: rcat+=[newcat[i]] 		# this is the slow way, but it's straight forward. otherwise, sort and look at rcat[-1]
			i+=1
		#
		rcat.sort()
		return rcat
		
	def getMapIndexMod(self, clustIndex=0, ndithers=None, nLon=None):
		# given an index from a dithered map, returns [mapID (aka, dithered hazmap object index), gridIndex (grid index from that hazmap)]
		# (a funky way of doing it) (and both of these methods seem to work... at lest they are identically flawed).
		if ndithers==None: ndithers=self.ndithers
		if nLon==None: nLon=self.getNlon()
		#
		i=clustIndex
		nd=ndithers
		mapIndex = (i%nd) + nd*int(int(i/nLon)%nd) + nd*nd*int((i%nLon)/nd)
		# but i think this is giving me only the x-position, not the proper y-pos.
		fullypos=i/self.getNlon()			# y position for dithered map.
		# add this separately to return value or lump into mapIndex? dunno just yet...
		ydithered=fullypos/self.ndithers	# y position for each map.
		mapGrid=mapIndex/(ndithers*ndithers) + (self.getNlon()/self.ndithers)*ydithered
		return [mapIndex%(ndithers*ndithers), mapGrid]

	def getMapIndex(self, clustIndex=0, ndithers=None, nLon=None):
		# given an index from a dithered map, returns [mapID (aka, dithered hazmap object index), gridIndex (grid index from that hazmap)]
		# (less elegant but more straight-forward) (and both of these methods seem to work... at lest they are identically flawed).
		if ndithers==None: ndithers=self.ndithers
		if nLon==None: nLon=self.getNlon()
		#
		x=clustIndex%nLon
		y=clustIndex/nLon		#x,y on dithered haz-map.
		#
		mapindex=x%ndithers + ndithers*(y%ndithers)	# hazardmap index (haxmaps[mapindex])
		#
		# and the index on the haz-map...
		mapx=x/ndithers
		mapy=y/ndithers
		gridindex=mapx + nLon*mapy/ndithers
		#
		return [mapindex, gridindex]	# so this is [the hazard map being used, grid index on that hazard map]
		
		
	def findClusters(self, gridlist=[], L=None):
		if gridlist==[] or gridlist==None:gridlist=map(operator.itemgetter(2), self.contourData())
		if L==None: L=self.getNlon()
		
		# find clusters and their size.
		# two passes total (i think) through the data:
		# - walk gridk "left" to "right" (aka, i -> i+1). L/R adjacent elements are elements of a cluster; record clusters in clustList=[[j0, [indeces]], [j1, [ind.]] ]
		# - if there is an element "below" (aka z[i-L]):
		#		- z[i] -> that cluster
		#		- members of z[i]'s cluster are tagged for renaming (add an entry to clustMap -> [j, j']. use new clust number going forward.
		# - after walking the grid, update clustList indeces from top to bottom according to clustMap, also top to bottom.
		#
		mygridlist=scipy.array(gridlist)
		clusterIDs=[]	# will be same dim. as gridlist, but will be populated with clusterID indices.
		clustList=[[0, []]]	#[[j0, [indices for j0]], [j1, [indeces for j1]...]	# and we might not need the index name j0, j1, ... we can probably just use order as index.
		clustMap=[]		# when we join two clusters: [[j0, j0'], [j1, j1'], ...]
		#
		# 1D or  2D array?
		#if type(gridlist[0]).__name__ in ('float', 'int')
		# more robustly?
		'''
		if type(gridlist[0])==type(5.0) or type(gridlist[0])==(5):
			# in case type names change in later pythons? aka, float64, intlong, or whatever. of course, this dos not guarantee that the flavor
			# of float/int we provide is the same.
			#
			# anyway, this is a 1D array. infer 2D from L.
			# let's just cast this into a 2D array (it's got to be one or the other).
			mygridlist.shape=(len(gridlist)/L, L)	# assuming config. [x0,y0; x1, y0; x2,y0... we should also trap for not integer division, etc.
		#
		# otherwise, it's a 2D array (or something that will make a big mess and break or produce nonsense).
		'''
		# but i seem to be happier with 1D array management, so let's do it this way:
		if type(gridlist[0]).__name__ in ('tuple', 'list', 'ndarray'):
			# elemtens are list-like, aka a 2+D array. assume 2D.
			L=len(gridlist[0])
			mygridlist.shape=(len(mygridlsit))
		#
		gridlen=len(gridlist)
		'''
		if L==None:
			# before giving up, guess a square:
			squareL=((float(gridlen))**.5)
			if squareL%1==0:
				# it could be a square...
				L=squareL
			else:
				# we don't know how to parse it.
				print("array width not defined.")
				return None
		'''
		#
		clustIndex=0
		#
		i=0
		while i<gridlen:
			if mygridlist[i]==None or mygridlist[i]<=0:
				clusterIDs+=[None]
				i+=1
				continue	# though we might just use NONE and to allow for 0-value perimeter padding.
			#
			# at this point, we've found an occupied site.
			#
			# is it a (potentially) new cluster on the left edge?
			#if i%L==0:
			# for starters, assume it is:
			myclustID=clustList[-1][0]+1
			if i%L>0:
				# not on left edge
				if clusterIDs[i-1]!=None:
					# join cluster to left
					myclustID=clusterIDs[i-1]
				
			if i>=L:		# a more efficient (aka, not always checking i>=L) approach might be to process the first row separately.
				# we've passed the first row. is it part of a prev-row cluster?
				if clusterIDs[i-L]!=None:
					myclustID=clusterIDs[i-L]	# join the cluster below.
					# was it part of a cluster to the left?
					if clusterIDs[i-1]!=None and clusterIDs[i-1]!=myclustID:
						# i suppose we could change this guy while we're here, but it is not necessary. write the label change:
						#clustMap+=[[myclustID, clusterIDs[i-1]]] #[[i_from, i_to]] as executed from the top of the list.
						clustMap+=[[clusterIDs[i-1], myclustID]] #[[i_from, i_to]] as executed from the top of the list.
				
					#
				#
			#
			# write to cluterIDS:
			clusterIDs+=[myclustID]
			#
			# add this element to its cluter list.
			if len(clustList)<(myclustID+1):
				# we've added a new cluster.
				clustList+=[[myclustID, []]]
			clustList[myclustID][1]+=[i]
			#
			i+=1
		# now, stitch together the clusters:	
		newclusts={}
		for i in xrange(len(clustList)):
			newclusts[i]=i
		clustMap.reverse()
		for rw in clustMap:
			srcindex=newclusts[rw[0]]
			targetindex=newclusts[rw[1]]
			#
			#clustList[rw[1]][1]+=clustList[rw[0]][1]
			#clustList[rw[0]][1]=[]
			if targetindex==srcindex:
				#print("target/source problem")
				continue
			clustList[targetindex][1]+=clustList[srcindex][1]
			clustList[srcindex][1]=[]
			newclusts[srcindex] = targetindex
		#
		rlist=[]
		for rw in clustList:
			if len(rw[1])==0: continue
			rlist+=[rw[1]]
		return rlist	
	
	def guessMagL(self, nx=None, ny=None, dithers=1.0, lat=None, gridsize=None, llkm=111.0, gamma=1.0, degs=1):
		# gridsize can be calculated from haz-map info.. eventually.
		rlat=lat
		if degs==1: rlat=2*math.pi*lat/360.
		#
		xsq=(nx*math.cos(rlat)*gridsize*llkm/dithers)**2. + (ny*gridsize*llkm/dithers)**2
		L=xsq**.5
		print("length: %f" % L)
		#
		m=2.*math.log10(L*gamma)+3.25
		return m
	
	def guessMag(self, nsquares=None, mc=None, deltas=2.0, beta=1.0, rblen=None, dithers=None):
		# we must have mc, rblen, dithers, nsquares. other prams we can guess.
		#if mc==None or nsquares==None or rblen==None or dithers==None: return None
		if mc==None: mc=self.mc
		if rblen==None: rblen=self.winlen
		if dithers==None: dithers=self.ndithers*self.ndithers
		if nsquares==None: nsquares=1	# minimun estimate, though maybe we should break it.
		#
		return mc + deltas + math.log10(beta*rblen*nsquares/float(dithers))

class rbHazMapLean(object):
	# 
	# a very lean hazard map - aka, just the data. primarily, this class is to be used in a dithered haz-map. all the graphical/mapping
	# bits will be in the parent class.
	
	#plt.ion()
	
	def __init__(self, catalog=None, catnum=0, mc=2.5, gridsize=.5, phi=[0,0], winlen=128, avlen=1, bigmag=5.0, logZ=1.0):
		# catalog: yp.eqcatalog object, phi: spatial phase shift(s) [lat,lon]
		self.catnum=catnum
		self.mc=mc
		self.gridsize=gridsize
		self.phi=phi
		self.winlen=winlen
		self.avlen=avlen
		self.bigmag=bigmag
		self.logZ=logZ
		#self.fignum=fignum
		#self.catalog=catalog
		#self.catalog=yp.eqcatalog(catalog.getcat(catnum))	# move to after "if self.catalog==None"
		#self.avlen=10
		self.precatindex=None
		#
		self.rbRatios=[]	#rb ratios for grid. this will either be a single array [r_0, r1, r2...] or an array of lists; we'll take the last val or average over the
								# last rblen vals [[r00, r01, r02...], [r10, r11, r12...], ... ]
		#
		if catalog==None:
			#catalog=getMexicaliCat(mc)
			self.catalog=getSocalCat(mc)
			# clean up subcats:
		else:
			self.catalog=yp.eqcatalog(catalog.getcat(catnum))
		#
		while len(self.catalog.subcats)>0: self.catalog.subcats.pop()
		#
		try:
			self.catalog.rb
		#except NameError:
		except:
			self.catalog.rb=rbi.intervalRecordBreaker(None)
		#
		#self.thiscat0=self.catalog.getcat(catnum)
		self.thiscat0=self.catalog.getcat(catnum)
		self.thiscat=self.thiscat0	# simplify this? do we need these two subcats?
		#self.thiscat=self.thiscat0
		#print("mainEvent: %s" % mev)
		#
		self.llrange=self.catalog.getLatLonRange(self.thiscat)
		# maps will stack better if we add phi to both sides of llrange:
		self.llrange[0]=[self.llrange[0][0]+phi[0], self.llrange[0][1]+phi[1]]
		self.llrange[1]=[self.llrange[1][0]+phi[0], self.llrange[1][1]+phi[1]]
		#
		self.deltaLat=abs(float(self.llrange[1][0])-self.llrange[0][0])
		self.deltaLon=abs(float(self.llrange[1][1])-self.llrange[0][1])
		self.Nlat=1+int(self.deltaLat/self.gridsize)
		self.Nlon=1+int(self.deltaLon/self.gridsize)	#number of cells in row/col
		#
		self.lat0=self.llrange[0][0]
		self.lon0=self.llrange[0][1]
		#
		ilat=0
		ilon=0	# index counters
		#
		# add empty sub-cats:
		#Z0=100	# z-order for squares
		for i in xrange(self.Nlat*self.Nlon):
			self.catalog.subcats+=[[i, []]]
			self.rbRatios+=[[None]]
		#	thissquare=getSquare([self.getLon(i), self.getLat(i)], self.gridsize)	# where [lon, lat] are the ll corner.
		
		#
		return None
		#
	def deleteSubcats(self):
		while len(self.catalog.subcats)>0: self.catalog.subcats.pop()
	
	def emptySubcats(self):
		for i in xrange(len(self.catalog.subcats)):
			self.catalog.subcats[i][1]=[]
			
	def getLat(self, indx):
		# clat=int(rw[0]/Nlon)*gridsize+lat0
		# clon=int(rw[0]%Nlon)*gridsize+lon0
		return float(int(indx/self.Nlon))*self.gridsize+self.lat0
	def getLon(self, indx):
		return float(int(indx%self.Nlon))*self.gridsize+self.lon0
	#	
	def hazMapTo(self, forecastDate=None, mev=None):
		# mev: provide a main-event if the default is not working for us.
		if forecastDate==None: return hazMap()
		#
		# otherwise, make a precat and a post cat.
		precat=[]
		postcat=[]
		#
		for rw in self.catalog.getcat(0):
			if rw[3]<self.mc: continue
			#print(rw[0])
			#print(forecastDate)
			if rw[0]<forecastDate: precat+=[rw]
			if rw[0]>=forecastDate: postcat+=[rw]
		#
		#self.precatindex=len(precat)-1
		self.hazMap(precat, postcat, mev)
		
		return None
	#
	def hazMap(self, preCat=None, postCat=None, mev=None):
		# set up the hazard map. divide the catalog into lat/lon bins; assign values to subcats (each bin is a subcat).
		if preCat==None: preCat=self.thiscat
		if postCat==None:
			if mev==None: mev=self.catalog.getMainEvent(self.thiscat0)
			postCat=self.thiscat0[mev[4]:]
		#if mev==None: mev=self.mev
		#
		self.precatindex=len(preCat)-1
		# any lingering subcats?
		self.emptySubcats()
		#
		# xy bin the catalog (note: nominally we can improve the speed by starting at the recent end of the catalog. of course, then we have to
		# determine how far back we go in the catalog. exactly how do we exclude areas with sparse seismicity?)
		#thiscat=self.thiscat
		rwcount=0
		#for rw in thiscat:
		# bin up the catalog... the problem here, of course, is that (i think), this
		# is not properly accounting for the latitude correction in bin size. this basically
		# needs to be rewritten and generalized to either determine N from the actual bin area
		# or to properly select all earthquakes r<R0 at a set of locations. note that this permits earthquakes to occupy
		# multiple bins for some spacings. it is probably fastest to use fairly sparse spacing so that after an earthquake
		# is placed into a bin it is removed from the catalog... then dither and sort like before.
		#
		for rw in preCat:
			#if rwcount>len(preCat): 
			#	print("stepping off pre-cat")
			#	continue
			if rw[3]<self.mc: continue
			ilat=(rw[1]-self.lat0)/self.gridsize
			ilat=int(ilat)
			ilon=(rw[2]-self.lon0)/self.gridsize
			ilon=int(ilon)
			scindex=ilat*self.Nlon + ilon	# sub-catalog index
			#
			# when we put phase into this, we can end up with events off the grid. for now, just dump these events. in the future,
			# cast a broader set of subcats.
			if scindex>=len(self.catalog.subcats): continue
	
			self.catalog.subcats[scindex][1]+=[rw]
			rwcount+=1
		#
		# catalog is xy binned. now, for each bin, get RB ratios and plot.
		self.updateHazMap()
	#	
	def updateHazMap(self, catalog=None):
		#
		# now, calc. statistics from each square.
		# [index, totalInterval, r_nrb]
		#gridStats=[]
		if catalog==None: catalog=self.catalog
		reds=[]
		blues=[]
		#print('as of %s' % str(self.thiscat[-1][0]))
		#cm=catalog.plotCatMap(thiscat, True, False, None, None, 'best', True, 'g,', None, self.fignum)
		#
		#return catalog
		for ct in self.catalog.subcats:
			catindex=ct[0]
			if len(ct[1])<=(self.winlen+self.avlen):
				self.rbRatios[catindex]=None
				# do we need to draw/undraw the square?
				continue
			#
			thisratios=catalog.rb.getIntervalRatios(minmag=self.mc, windowLen=self.winlen, cat0=ct[1][-(self.winlen+self.avlen):], logZ=self.logZ)
			rs=map(operator.itemgetter(-1), thisratios)
			#meanr=sum(rs[-self.avlen:])/float(self.avlen)
			# we need to take the log before we average:
			for i in xrange(len(rs)):
				rs[i]=math.log10(rs[i])
			meanr=numpy.average(rs[-self.avlen:])
			#self.rbRatios[catindex]=math.log10(meanr)
			self.rbRatios[catindex]=meanr
			#
			#rw=[ct[0], plt.date2num(ct[1][-1][0])-plt.date2num(ct[1][0][0]), meanr]
			#
			# boxes are already drawn we only have to update the color/alpha
			#clat=self.getLat(rw[0]) #int(rw[0]/self.Nlon)*self.gridsize+self.lat0
			#clon=self.getLon(rw[0]) #int(rw[0]%self.Nlon)*self.gridsize+self.lon0
			#
	#
	#
	def simplePlot(self, fignum=0):
		# do a simple matplotlib plot of hazmap. this is probably diagnostic, so it may just be points of some sort.
		#
		plt.figure(fignum)
		plt.clf()
		plt.ion()
		#
		plt.title("diagnostic plot of LEAN rb-hazmaps")
		
		llrange=self.llrange
		lonsRed=[]
		latsRed=[]
		#rbvals=[]
		lonsBlue=[]
		latsBlue=[]
		#
		lonsCyan=[]
		latsCyan=[]	# for None , at least temporarily
		#
		for i in xrange(len(self.rbRatios)):
			#
			#x=llrange[0][1]+self.gridsize*(i%hm.Nlon)			#lon
			#y=llrange[0][0]+self.gridsize*int(i/hm.Nlon)		# lat
			x=self.getLon(i)
			y=self.getLat(i)
			if self.rbRatios[i]==None:
				continue
				#lonsCyan+=[x]
				#latsCyan+=[y]
			#
			if self.rbRatios[i]<=0:
				lonsRed+=[x]
				latsRed+=[y]
			elif self.rbRatios[i]>0:
				lonsBlue+=[x]
				latsBlue+=[y]
		plt.plot(lonsRed, latsRed, 'r%s' % yp.pyicon(0))
		plt.plot(lonsBlue, latsBlue, 'b%s' % yp.pyicon(0))
		#plt.plot(lonsCyan, latsCyan, '%s' % yp.pyicon(0), color='y')
		#print(hmindex)
			
	def getSquare(self, xy=[], dx=.1):
		# xy = ll corner.
		return [xy, [xy[0], xy[1]+dx], [xy[0]+dx, xy[1]+dx], [xy[0]+dx, xy[1]], xy]	
	


class rbHazMap(object):
	# 
	# i think i'm going to move this to it's own module common/hazmap.py
	#plt.ion()
	
	def __init__(self, catalog=None, catnum=0, mc=2.5, gridsize=.5, phi=[0,0], winlen=128, avlen=1, bigmag=5.0, fignum=0):
		# catalog: yp.eqcatalog object, phi: spatial phase shift(s) [lat,lon]
		self.catnum=catnum
		self.mc=mc
		self.gridsize=gridsize
		self.phi=phi
		self.winlen=winlen
		self.avlen=avlen
		self.bigmag=bigmag
		self.fignum=fignum
		#self.catalog=catalog
		#self.catalog=yp.eqcatalog(catalog.getcat(catnum))	# move to after "if self.catalog==None"
		self.avlen=10
		self.precatindex=None
		#
		self.gridSquares=[]
		#
		self.fig=plt.figure(self.fignum)
		self.fig.clf()
		self.ax  = self.fig.add_subplot(111)
		plt.title('Record-breaking intervals Hazard Map')
		if catalog==None:
			#catalog=getMexicaliCat(mc)
			self.catalog=getSocalCat(mc)
			# clean up subcats:
		else:
			self.catalog=yp.eqcatalog(catalog.getcat(catnum))
		#
		while len(self.catalog.subcats)>0: self.catalog.subcats.pop()
		#
		try:
			self.catalog.rb
		#except NameError:
		except:
			self.catalog.rb=rbi.intervalRecordBreaker(None)
		#
		#self.thiscat0=self.catalog.getcat(catnum)
		self.thiscat0=self.catalog.getcat(0)
		mev=self.catalog.getMainEvent(self.thiscat0)
		print(mev)
		self.mev=mev
		self.thiscat=self.thiscat0[0:mev[4]-1]	# eliminate all events after mainshock...
		self.precatindex=len(self.thiscat)-1
		#self.thiscat=self.thiscat0
		#print("mainEvent: %s" % mev)
		#
		self.llrange=self.catalog.getLatLonRange(self.thiscat)
		# most things should work better if we add phi to both sides of llrange (not just the source side).
		self.llrange[0]=[self.llrange[0][0]+phi[0], self.llrange[0][1]+phi[1]]
		self.llrange[1]=[self.llrange[1][0]+phi[0], self.llrange[1][1]+phi[1]]
		#
		self.deltaLat=abs(float(self.llrange[1][0])-self.llrange[0][0])
		self.deltaLon=abs(float(self.llrange[1][1])-self.llrange[0][1])
		self.Nlat=1+int(self.deltaLat/self.gridsize)
		self.Nlon=1+int(self.deltaLon/self.gridsize)	#number of cells in row/col
		#
		self.lat0=self.llrange[0][0]
		self.lon0=self.llrange[0][1]
		#
		ilat=0
		ilon=0	# index counters
		#self.cm=self.catalog.plotCatMap(self.thiscat, True, False, None, None, 'best', True, 'g,', None, self.fignum)
		#self.cm=self.catalog.plotCatMap(self.thiscat, True, False, None, None, 'best', True, 'g,', self.ax)
		#minicat=self.catalog.getcat(0)[0:10]
		thismev=self.catalog.getMainEvent(self.catalog.getcat(0))
		epicen=[thismev[2], thismev[1]]
		self.cm=self.catalog.plotCatMap(self.catalog.getcat(catnum), True, False, None, epicen, 'best', True, 'g,', self.ax)
		#self.cm=self.catalog.plotCatMap([], True, False, None, epicen, 'best', True, 'g,', self.ax)
		self.mapEvents=self.ax.plot([], [], 'g,', zorder=1)[0]	# line object of currently plotted data.
		self.plottedCat=[]	#[x,y]	# data currently plotted
		
		# add empty sub-cats:
		#Z0=100	# z-order for squares
		for i in xrange(self.Nlat*self.Nlon):
			self.catalog.subcats+=[[i, []]]
			thissquare=getSquare([self.getLon(i), self.getLat(i)], self.gridsize)
			mapX, mapY=self.cm(map(operator.itemgetter(0), thissquare), map(operator.itemgetter(1), thissquare))
			self.gridSquares+=self.ax.fill(mapX, mapY, alpha=.1, color='b', lw=2, zorder=2)
		self.legend=self.ax.legend(loc='best', numpoints=1)
		self.fig.canvas.draw()
		# and manipulate the map like this:
		#for i in xrange(len(self.gridSquares)):
		#	self.gridSquares[i].set_color('b')
		#	self.gridSquares[i].set_alpha(.1)
		#self.fig.canvas.draw()
		#
		return None
		#
	def deleteSubcats(self):
		while len(self.catalog.subcats)>0: self.catalog.subcats.pop()
	
	def emptySubcats(self):
		for i in xrange(len(self.catalog.subcats)):
			self.catalog.subcats[i][1]=[]
			
	def getLat(self, indx):
		# clat=int(rw[0]/Nlon)*gridsize+lat0
		# clon=int(rw[0]%Nlon)*gridsize+lon0
		return float(int(indx/self.Nlon))*self.gridsize+self.lat0
	def getLon(self, indx):
		return float(int(indx%self.Nlon))*self.gridsize+self.lon0
	#	
	def hazMapTo(self, forecastDate=None, mev=None):
		# mev: provide a main-event if the default is not working for us.
		if forecastDate==None: return hazMap()
		#
		# otherwise, make a precat and a post cat.
		precat=[]
		postcat=[]
		#
		for rw in self.catalog.getcat(0):
			if rw[3]<self.mc: continue
			if rw[0]<forecastDate: precat+=[rw]
			if rw[0]>=forecastDate: postcat+=[rw]
		#
		#self.precatindex=len(precat)-1
		self.hazMap(precat, postcat, mev)
		
		return None
	#
	def hazMap(self, preCat=None, postCat=None, mev=None):
		if preCat==None: preCat=self.thiscat
		if postCat==None:
			if mev==None: mev=self.catalog.getMainEvent(self.thiscat0)
			postCat=self.thiscat0[mev[4]:]
		if mev==None: mev=self.mev
		
		self.precatindex=len(preCat)-1
		# any lingering subcats?
		self.emptySubcats()
		#
		# xy bin the catalog (note: nominally we can improve the speed by starting at the recent end of the catalog. of course, then we have to
		# determine how far back we go in the catalog. exactly how do we exclude areas with sparse seismicity?)
		#thiscat=self.thiscat
		rwcount=0
		#for rw in thiscat:
		for rw in preCat:
			#if rwcount>len(preCat): 
			#	print("stepping off pre-cat")
			#	continue
			if rw[3]<self.mc: continue
			ilat=(rw[1]-self.lat0)/self.gridsize
			ilat=int(ilat)
			ilon=(rw[2]-self.lon0)/self.gridsize
			ilon=int(ilon)
			scindex=ilat*self.Nlon + ilon	# sub-catalog index
			#
			# when we put phase into this, we can end up with events off the grid. for now, just dump these events. in the future,
			# cast a broader set of subcats.
			if scindex>=len(self.catalog.subcats): continue
	
			self.catalog.subcats[scindex][1]+=[rw]
			rwcount+=1
		#
		# catalog is xy binned. now, for each bin, get RB ratios and plot.
		self.updateHazMap()
		# and where did the earthquakes happen?
		bigeqs=[]
		#for rw in self.thiscat0[self.mev[4]-1:]:
		for rw in postCat:
			if rw[3]>=self.bigmag:
				# bigeqs+=[rw]
				x,y=self.cm(rw[2], rw[1])
				self.cm.plot(x,y, '*', ms=15, label='m%f, %s' % (rw[3], str(rw[0])))
				#print('m%f, %s' % (rw[3], str(rw[0])))
		#x,y=self.cm(self.mev[2], self.mev[1])
		#print(x,y)
		x,y=self.cm(mev[2], mev[1])
		#print(x,y)
		self.cm.plot(x,y, 'r*', ms=25,  label='ms:%f, %s' % (rw[3], str(rw[0])))		
		#
		#plt.gca().set_title('dtm: %s' % str(preCat[-1][0]))
		#plt.figure(self.fignum)
		#plt.title('dtm: %s' % str(preCat[-1][0]))
		self.ax.set_title('dtm: %s' % str(preCat[-1][0]))
		#
		plt.legend(loc='best')
		
	
	def hazMapMovie(self, nRes=100, cat=None, imagesDir=None):
		#thiscat=self.thiscat
		if cat==None: cat=self.thiscat0
		if nRes==None: nRes=100	# number of events to skip.
		#
		# any lingering subcats?
		self.emptySubcats()
		#
		# and where did the earthquakes happen?
#		bigeqs=[]
#		#for rw in self.thiscat0[self.mev[4]-1:]:
#		for rw in postCat:
#			if rw[3]>=self.bigmag:
#				# bigeqs+=[rw]
#				x,y=self.cm(rw[2], rw[1])
#				self.cm.plot(x,y, '*', ms=15, label='m%f, %s' % (rw[3], str(rw[0])))
#		x,y=self.cm(self.mev[2], self.mev[1])
#		self.cm.plot(x,y, 'r*', ms=25,  label='ms:%f, %s' % (mev[3], str(mev[0])))		
#		#
#		self.ax.legend(loc='best', numpoints=1)
		# xy bin the catalog (note: nominally we can improve the speed by starting at the recent end of the catalog. of course, then we have to
		# determine how far back we go in the catalog. exactly how do we exclude areas with sparse seismicity?)
		#
		# catalog is xy binned. now, for each bin, get RB ratios and plot.
#		self.updateHazMap()
		#
		# set empty marker-place-holders for "big" earthquakes.
		NbigEqs=5
		self.nextEqs=[]
		for i in xrange(NbigEqs):
			self.nextEqs+=self.ax.plot([], [], '*', ms=15, label='', zorder=3)
		self.ax.legend(loc='best', numpoints=1)
		
				
		rwcount=0
		#for rw in thiscat:
		maxFrames=1+int(len(cat))
		numLen=1+int(math.log10(maxFrames))
		for rw in cat:
			if rw[3]<self.mc: 
				continue
				rwcount+=1	# but we have to advance rwcount because we use it as an index relative to the full catalog.
			self.plottedCat+=[[rw[2], rw[1]]]	# earthquake position
			ilat=(rw[1]-self.lat0)/self.gridsize	#lat-index
			ilat=int(ilat)
			ilon=(rw[2]-self.lon0)/self.gridsize	#lon-index
			ilon=int(ilon)
			scindex=ilat*self.Nlon + ilon	# sub-catalog index
	
			self.catalog.subcats[scindex][1]+=[rw]		# add earthquake to appropriate grid-subcat
			#
			if rwcount%nRes==0:
				#plt.title('date=%s' % str(rw[0]))
				#self.ax.text(.4, -.2, 'dtm: %s' % str(rw[0]))
				x,y=self.cm(map(operator.itemgetter(0), self.plottedCat), map(operator.itemgetter(1), self.plottedCat))
				self.mapEvents.set_data(x,y)
				self.ax.set_title('dtm: %s' % str(rw[0]))
				self.updateHazMap()
				#self.updateBigEqs(cat[rwcount:rwcount+500], self.bigmag)
				#self.updateBigEqs2(cat, rwcount, self.bigmag, NbigEqs)
				self.updateBigEqList(cat, rwcount, self.bigmag)
				if imagesDir!=None:
					framenum=str(rwcount)
					while len(framenum)<numLen: framenum='0'+framenum
					fname='%s/hazMapMovie-%s.png' %(imagesDir, framenum)
					fname.replace('//', '/')
					self.fig.savefig(fname)
					os.system('convert %s %s' % (fname, fname.replace('.png', '.jpg')))
					
			rwcount+=1
			#
		#
		# os.system('mencoder mf:///home/myoder/Documents/Research/Rundle/econoPhysics/econoDensityWave/images1/*.jpg -mf w=800:h=600:fps=10:type=jpg -ovc copy -o econoDwaveMovie.avi')
		if imagesDir!=None: os.system('mencoder mf://%s/%s/*.jpg -mf w=800:h=600:fps=10:type=jpg -ovc copy -o hazmapmovie2.avi' % (os.getcwd(), imagesDir))

		
	def updateBigEqs(self, cat=None, bigmag=None):
		if cat==None: cat=self.thiscat
		if bigmag==None: bigmag=self.bigmag
		#
		for rw in cat:
			if rw[3]>=bigmag:
				x,y=self.cm(rw[2], rw[1])
				strMags=str(round(rw[3],2)).split('.')
				strLegend='m%s.%s, %d/%d/%d' % (strMags[0], strMags[1][0:2], rw[0].year, rw[0].month, rw[0].day)
				self.cm.plot(x,y, '*', ms=15, label='m%f, %s' % (rw[3], str(rw[0])))
		#self.ax.legend(loc='best', numpoints=1)
	
	def updateBigEqs2(self, cat=None, rwindex=0, bigmag=None, Nbigeq=5):
		if cat==None: cat=self.thiscat
		if bigmag==None: bigmag=self.bigmag
		#
		neq=0
		for rw in cat:
			if rw[3]>=bigmag:
				x,y=self.cm(rw[2], rw[1])
				strMags=str(round(rw[3],2)).split('.')
				strLegend='m%s.%s, %d/%d/%d' % (strMags[0], strMags[1][0:2], rw[0].year, rw[0].month, rw[0].day)
				self.cm.plot(x,y, '*', ms=15, label='m%f, %s' % (rw[3], str(rw[0])))
				neq+=1
			if neq>=Nbigeq: break
	
	def updateBigEqList(self, cat=None, rwindex=0, bigmag=None):
		if cat==None: cat=self.thiscat0
		if bigmag==None: bigmag=self.bigmag
		Nbigeq=len(self.nextEqs)
		#
		neq=0
		for i in xrange(rwindex, len(cat)):
		#for rw in cat:
			rw=cat[i]
			if rw[3]>=bigmag:
				x,y=self.cm(rw[2], rw[1])
				strMags=str(round(rw[3],2)).split('.')
				strLegend='m%s.%s, %d/%d/%d' % (strMags[0], strMags[1][0:2], rw[0].year, rw[0].month, rw[0].day)
				#self.cm.plot(x,y, '*', ms=15, label='m%f, %s' % (rw[3], str(rw[0])))
				self.nextEqs[neq].set_data(x,y)
				self.nextEqs[neq].set_label(strLegend)
				neq+=1
			if neq>=Nbigeq: break
		self.ax.legend(loc='lower left')
	
	def updateHazMap(self, catalog=None):
		#
		# now, calc. statistics from each square.
		# [index, totalInterval, r_nrb]
		self.rbRatios=[]
		gridStats=[]
		if catalog==None: catalog=self.catalog
		reds=[]
		blues=[]
		#print('as of %s' % str(self.thiscat[-1][0]))
		#cm=catalog.plotCatMap(thiscat, True, False, None, None, 'best', True, 'g,', None, self.fignum)
		#
		#return catalog
		for ct in self.catalog.subcats:
			catindex=ct[0]
			if len(ct[1])<=(self.winlen+self.avlen):
				self.gridSquares[catindex].set_alpha(0)
				# do we need to draw/undraw the square?
				continue
			#
			thisratios=catalog.rb.getIntervalRatios(self.mc, self.winlen, ct[1][-(self.winlen+self.avlen-1):])
			rs=map(operator.itemgetter(-1), thisratios)
			for i in xrange(len(rs)):
				rs[i]=math.log10(rs[i])
			#meanr=sum(rs[-self.avlen:])/float(self.avlen)
			meanr=numpy.average(rs[-self.avlen:])
			self.rbRatios=meanr
			#####
			#gridStats+=[[ct[0], plt.date2num(ct[1][-1][0])-plt.date2num(ct[1][0][0]), meanr]]
			rw=[ct[0], plt.date2num(ct[1][-1][0])-plt.date2num(ct[1][0][0]), meanr]
			#
			# boxes are already drawn we only have to update the color/alpha
			#clat=self.getLat(rw[0]) #int(rw[0]/self.Nlon)*self.gridsize+self.lat0
			#clon=self.getLon(rw[0]) #int(rw[0]%self.Nlon)*self.gridsize+self.lon0
			#
			if meanr<=0: boxcolor='r'
			if meanr>0: boxcolor='b'
			
			self.gridSquares[catindex].set_alpha(.4)
			self.gridSquares[catindex].set_color(boxcolor)
		#
		self.fig.canvas.draw()

######
# utilities:
#
def getWinlen(mc=4.5, gridsize=1.0, lat=45., sigmaExp=1.68):
	#
	# this is in lat/lon sizes
	# note: R_earth = 6378 km
	# dy= R sin dLat = R*2*pi*dLat_deg/360
	# dx = (R*2*pi/360)*dLon_deg*cos(Lat)
	#
	Rearth=6378.1	#km
	myfactor=Rearth*2.*math.pi/360.	# this is basically (mostly) the degrees to radians conversion.
	lat=2*math.pi*lat/360.
	
	thisarea=((gridsize*myfactor)**2)*math.cos(lat)	# in km^2 for equal lat/lon size bins.
	winlen=round((10**(sigmaExp-mc))*thisarea)
	#
	return winlen

def getgridsize(mc=4.5, winlen=16, lat=45., sigmaExp=1.68):
	# same as getWinlen but inverted...
	Rearth=6378.1	#km
	myfactor=360./(Rearth*2.*math.pi)	# this is basically (mostly) the degrees to radians conversion.
	lat=2*math.pi*lat/360.	# in radians...
	
	thisGS=myfactor*((winlen*10**(mc-sigmaExp))*(1./math.cos(lat)))**.5
	
	return thisGS

def getMidLat(objcat, mc=None):
	#
	lats=map(operator.itemgetter(1), objcat.cat)
	midlat=numpy.mean(lats)
	#
	return midlat

def crackMap(objCat=None, catnum=0, mc=2.5, gridsize=None, ndithers=4, winlen=20, sigma=1.68, avlen=None, todt=dtm.datetime(2010, 4, 1, 0, 0, 0, 0, pytz.timezone('UTC')), epicen=None, thresh=0, logZ=None):
	# one mc, gridsize combination:mc=2.5, gridsize=.125, winlen=20 (approximately)
	# choose an optimal combination:
	# get gridsize/winlen:
	# getgridsize(mc=4.5, winlen=16, lat=45., sigmaExp=1.68)
	# def getWinlen(mc=4.5, gridsize=1.0, lat=45., sigmaExp=1.68):
	#
	if avlen==None:
		avlen=int(winlen/10)
		if avlen==0: avlen=1
	# get gridsize/winlen:
	# getgridsize(mc=4.5, winlen=16, lat=45., sigmaExp=1.68)
	# def getWinlen(mc=4.5, gridsize=1.0, lat=45., sigmaExp=1.68):
	if gridsize==None and winlen==None:
		# can't do anything with that, can we?
		return None
	midlat=getMidLat(objCat, mc)
	if gridsize==None:
		gridsize=getgridsize(mc, winlen, midlat, sigma)
	if winlen==None:
		winlen=getWinlen(mc, gridsize, midlat, sigma)
	#
	bigmag=None
	fignum=0	# but we don't use it; just keep it for a place holder
	#
	# def __init__(self, catalog=None, catnum=0, mc=2.5, gridsize=.5, ndithers=3, winlen=128, avlen=1, bigmag=5.0, fignum=0, logZ=1.0):
	objHM=hazCrack(catalog=objCat, catnum=catnum, mc=mc, gridsize=gridsize, ndithers=ndithers, winlen=winlen, avlen=avlen, bigmag=bigmag, fignum=fignum, logZ=logZ)
	#
	#objHM.epicen=epicen
	objHM.hazMapTo(todt)
	#
	#
	magres=objHM.guessMag(nsquares=1, mc=mc, deltas=2.0, beta=1.0, rblen=winlen, dithers=ndithers**2)
	magresDith=objHM.guessMag(nsquares=ndithers**2, mc=mc, deltas=2.0, beta=1.0, rblen=winlen, dithers=ndithers**2)
	#print("mag-res: %f/%f (dithering=%d)" % (magres, magresDith, ndithers))
	#
	#plt.figure(fignum)
	#objHM.simpleBoxes(fignum=fignum, thresh=thresh)
	#
	# el mayor (baja, mexicali) earthquake:
	#hmx, hmy=objHM.hazmaps[0].catalog.catmap(-115.3, 32.13)
	#plt.plot([hmx], [hmy], 'b*', ms=18, zorder=10, alpha=.8)
	
	#z2=objHM.contourMap(fignum+2)
	#plt.figure(fignum+2)
	#plt.plot([hmx], [hmy], 'b*', ms=18)
	#
	return objHM
