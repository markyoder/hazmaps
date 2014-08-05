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

plt.ion()

nContours=[-.9, -.8, -.7, -.6, -.5, -.4, -.3, -.2, -.1, 0.0, .1, .2, .3, .4, .5, .6, .7]
nContours1=[-.8, -.7, -.6, -.5, -.4, -.3, -.2, -.1, 0.0, .1, .2, .3, .4, .5]

def socalRegionalQuad(targmag=7.2, mc=2.5, rfactor=.5, fnum=0, refreshcat=True):
	N0 = getNofm(targmag, mc)
	emcEvent = [dtm.datetime(2010, 4, 4, 22, 40, 42, 360004, tzinfo=mhp.pytz.timezone('UTC')), 32.2862, -115.2953, 7.2, 10.0]
	emcFS    = [dtm.datetime(2009, 12, 30, 18, 48, 57, 330004, tzinfo=mhp.pytz.timezone('UTC')), 32.464, -115.1892, 5.8, 6.0]	#"foreshock"?
	a=mhp.elmayorhm(mc=mc, todt=mhp.dtm.datetime.now(mhp.pytz.timezone('UTC')), ndithers=5, nContours=nContours1, bigmag=9.0, refreshcat=refreshcat, dt0=mhp.dtm.datetime(1990,1,1, tzinfo=mhp.pytz.timezone('UTC')))
	c1=a.getcat()
	while len(c1.subcats)>0:c1.subcats.pop()
	#
	# now, plot all the big earthquakes, little earthquakes, etc.
	c1=a.getcat()
	c1.rbomoriQuadPlot(mc=mc, targmag=targmag, bigmag=8.0, catnum=len(c1.subcats), fignum=fnum, plotevents=False, intlist=[int(N0*.5), int(N0*.666), N0, int(N0*1.5)], thislw=2.25)
	
	f=plt.figure(fnum)
	thisAxes = f.get_axes()
	#biggies = []
	for ev in c1.getcat(0):
		if ev[0]<dtm.datetime(1990,1,1, tzinfo=mhp.pytz.timezone('UTC')) or ev[3]<7.0: continue
		#biggies+=[ev]
		X,Y = c1.catmap(ev[2], ev[1])
		ax = thisAxes[3]
		ax.plot([X], [Y], '*', ms=15, alpha=.8, zorder=11, label='m=%.2f, %s' % (ev[3], shortdate(ev[0])))
		#
		ax = thisAxes[0]
		ax.plot([ev[0]], [ev[3]], '*', ms=15, label='m=%.2f, %s' % (ev[3], shortdate(ev[0])))
		#
		ax = thisAxes[1]
		ax.plot([ev[0], ev[0]], [0., 1.8], 'c-', lw=2.5, ms=15)
		ax = thisAxes[2]
		ax.plot([ev[0], ev[0]], [1./9., 9.], 'c-', lw=2.5, ms=15)
		
	#
	f=plt.figure(fnum)
	thisAxes = f.get_axes()
	cat_transposed = zip(*c1.getcat())
	Xev, Yev = c1.catmap(cat_transposed[2], cat_transposed[1])
	ax = thisAxes[3]
	ax.plot(Xev, Yev, 'b.', ms=3, alpha=.65, zorder=10)
	#
	for myaxes in [0,2,3]:
		ax=thisAxes[myaxes]
		ax.legend(loc=0, numpoints=1)
	
def shortdate(dtm):
	rstring = '%d-%d-%d'	% (dtm.month, dtm.day, dtm.year)
	return rstring

def EMCTriple(targmag=7.2, mc=2.5, rfactor=.5, fignum=0, lats=[31.5, 33.0], lons=[-116.15, -114.5]):
	# a triple-plot + catalog
	
	#lats=[31.5, 33.5]
	#lons=[-116.65, -114.5]
	
	plt.figure(fignum)
	plt.clf()
	#ax0=plt.axes([.1,.1,.85, .35])
	# define axis (subplots) boundaries:
	ydelim=.03
	xdelim=.05
	xTS=[0.05, .5]
	yTS0=0.05
	dyTS=.3
	#
	# map variables:
	lat_sep=.25
	lon_sep=.25
	#
	rb_len = winlen(targmag, mc, mt=7.6, doInt=False)
	rb_len=int(round(rb_len,-1))
	if rb_len<1: rb_len=1
	#
	rbavelen=int(rb_len/10.)
	if rbavelen==0: rbavelen=1
	#
	# some earthquake trivia (mainshock and large aftershocks):
	# get a catalog:
	catlist=atp.catfromANSS(dates0=[dtm.datetime(1990, 1, 1, 0, 0, 0, 0, tzinfo=mhp.pytz.timezone('UTC')), dtm.datetime.now(mhp.pytz.timezone('UTC'))], minMag=mc, lat=lats, lon=lons)
	quakes=[]
	quakes+=[[dtm.datetime(2010, 4, 4, 22, 40, 42, 360004, tzinfo=pytz.timezone('UTC')), 32.2862, -115.2953, 7.2, 10.0]]
	quakes+=[[dtm.datetime(2010, 4, 4, 22, 43, 0, 429995, tzinfo=pytz.timezone('UTC')), 32.4533, -115.632, 5.26, 10.0]]
	quakes+=[[dtm.datetime(2010, 4, 4, 22, 50, 17, 120002, tzinfo=pytz.timezone('UTC')), 32.0987, -115.0482, 5.7, 10.0]]
	quakes+=[[dtm.datetime(2010, 4, 4, 23, 15, 14, 239995, tzinfo=pytz.timezone('UTC')), 32.3, -115.2595, 5.43, 10.0]]
	quakes+=[[dtm.datetime(2010, 4, 4, 23, 25, 7, 189997, tzinfo=pytz.timezone('UTC')), 32.2662, -115.2925, 5.38, 10.0]]
	quakes+=[[dtm.datetime(2010, 4, 8, 16, 44, 25, 100004, tzinfo=pytz.timezone('UTC')), 32.1647, -115.2683, 5.29, 10.0]]
	quakes+=[[dtm.datetime(2010, 6, 15, 4, 26, 58, 480000, tzinfo=pytz.timezone('UTC')), 32.7002, -115.9213, 5.72, 5.41]]
	quakes+=[[dtm.datetime(2011, 2, 18, 17, 47, 35, 769999, tzinfo=pytz.timezone('UTC')), 32.047, -115.0622, 5.07, 15.0]]
	mainshock=quakes[0]
	objCat=eqp.eqcatalog()
	#
	myaxes=[]
	nax=0
	# mags:
	x0=xTS[0]
	y0=yTS0+dyTS*nax
	Dy=.275
	dy=.05
	#myaxes+=[plt.axes([xTS[0], y0, xTS[1], dyTS])]
	myaxes+=[plt.axes([.05, .05, .4, .925])]
	nax+=1
	myaxes+=[plt.axes([.5, 3.0*dy+2.0*Dy, .47, Dy])]
	myaxes[-1].text(.25, .8, 'North-west catalog', horizontalalignment='center', verticalalignment='center', transform=myaxes[-1].transAxes, fontsize=15)
	nax+=1
	myaxes+=[plt.axes([.5, 2.0*dy+Dy, .47, Dy], sharex=myaxes[1])]
	myaxes[-1].text(.35, .8, 'Center catalog', horizontalalignment='center', verticalalignment='center', transform=myaxes[-1].transAxes, fontsize=15)
	nax+=1
	myaxes+=[plt.axes([.5, dy, .47, Dy], sharex=myaxes[1])]
	myaxes[-1].text(.75, .8, 'South-east catalog', horizontalalignment='center', verticalalignment='center', transform=myaxes[-1].transAxes, fontsize=15)
	nax+=1
	#
	axescolor  = None #'#f6f6f6'
	myaxes+=[None]
	#myaxes+=[plt.axes([.51, 3.0*dy+2.0*Dy + Dy/3, .25, .2], frameon=False, axisbg=axescolor)]
	nax+=1
	myaxes+=[None]
	#myaxes+=[plt.axes([.75, 2.0*dy+1.0*Dy + dy, .25, Dy/3.0], frameon=False, axisbg=axescolor)]
	nax+=1
	myaxes+=[plt.axes([.45, 1.0*dy+0.0*Dy + Dy/2., .25, .2], frameon=True, axisbg=axescolor)]
	myaxes[-1].text(.25, .8, 'South-east (zoom)', horizontalalignment='center', verticalalignment='center', transform=myaxes[-1].transAxes, fontsize=15)
	nax+=1
	#
	lat0=min(lats)
	lat1=max(lats)
	lon0=min(lons)
	lon1=max(lons)
	bm=bmp.Basemap(llcrnrlon=lon0, llcrnrlat=lat0, urcrnrlon=lon1, urcrnrlat=lat1, projection='tmerc', resolution='i', ax=myaxes[0], lat_0=(lat0+.5*(lat1-lat0)), lon_0=(lon0+.5*(lon1-lon0)))
	bm.drawcoastlines()
	#bm.drawlsmask(land_color='beige', ocean_color='cyan', lakes=True, resolution='i')
	bm.fillcontinents(color='beige', lake_color='aqua')
	#bm.bluemarble()
	#bm.shadedrelief()
	#bm.etopo(scale=1.0)
	bm.drawstates()
	bm.drawcountries()
	bm.drawrivers(color='aqua')
	first_parallel = lat0-lat0%lat_sep
	last_parallel =  lat1-lat1%lat_sep
	first_merid = lon0-lon0%lon_sep
	last_merid  = lon1-lon1%lon_sep
	bm.drawparallels(scipy.arange(first_parallel, last_parallel, lat_sep), labels=[1,0,0,0])
	bm.drawmeridians(scipy.arange(first_merid, last_merid, lon_sep), labels=[0,0,0,1])
	#
	befores,afters=[],[]
	for rw in catlist:
		if rw[0]>mainshock[0]: afters+= [rw]
		if rw[0]<mainshock[0]: befores+=[rw]
	#
	x,y = bm(zip(*befores)[2], zip(*befores)[1])
	bm.plot(x,y, 'b.', ms=3, label='befores', alpha=.5, zorder=8)
	x,y = bm(zip(*afters)[2], zip(*afters)[1])
	bm.plot(x,y, 'g.', ms=3, label='afters', alpha=.5, zorder=9)
	#
	x,y = bm(zip(*quakes)[2], zip(*quakes)[1])
	bm.plot(x,y, 'm*', ms=25, alpha=.7, zorder=10)
	#
	center_SE = [-115.065487, 32.093364]
	center_0  = [-115.295300, 32.286200]
	center_NW = [-115.601718, 32.543315]
	#center_NW = [-115.632, 32.4533]
	center_NW2 = [-115.9213, 32.7002]
	centers=zip(*[center_NW, center_0, center_SE])
	#
	x,y = bm(centers[0], centers[1])
	bm.plot(x,y, 'ro', ms=12, alpha=1.0, zorder=11)
	#
	N0=int(getNofm(m=targmag, mc=mc, mt=7.6, dmstar=1.0, dms=1.0))
	Lr = 10.0**(targmag/2.0 - 1.76)
	#
	i=1
	centers=zip(*centers)
	for center in centers:
		circ_cat = circularcat(incat=catlist, latlon=[center[1], center[0]], Rkm=.5*Lr)
		#c1=eqp.eqcatalog(circ_cat)
		ratios = objCat.plotIntervalRatiosAx(winlen=rb_len, cat=circ_cat, hitThreshold=1.0, bigmag=9.0, thisAx=myaxes[i], ratios=None, delta_t=1, avlen=rbavelen, mainEv=mainshock, logZ=None, rbLegLoc=0, reverse=False)
		#print "ratios: ", len(ratios), len(ratios[0])
		#return ratios
		rs=map(operator.itemgetter(5), ratios)
		myaxes[i].set_ylim(bottom=min(rs), top=max(rs))
		#
		#myaxes[i].plot([mainshock[0]+dtm.timedelta(seconds=100), mainshock[0]+dtm.timedelta(seconds=100)], [.1, 10.], 'k-', lw=3, alpha=.8, zorder=11)
		myaxes[i].plot([mainshock[0]], [1.15], 'rd', ms=8, alpha=.8, zorder=11)
		myaxes[i].plot([mainshock[0]], [1.15], 'kd', ms=10, alpha=.7, zorder=10)
		#
		#
		circXY = getLLcircle(center[0], center[1], r=.5*Lr*1000., rtype=0, N=500)
		Xcirc, Ycirc=bm(scipy.array(circXY[0]), scipy.array(circXY[1]))
		bm.plot(Xcirc, Ycirc, 'r--', ax=myaxes[0], zorder=8, lw=3, alpha=.8)
		#ax.plot(Xcirc, Ycirc, 'r--', zorder=8, lw=2, alpha=.8)
		#llcenter = c1.catmap(thisCenter[0], thisCenter[1])
		#llepicen = c1.catmap(emcEvent[2], emcEvent[1])
		#llfs     = c1.catmap(emcFS[2], emcFS[1])
		#ax.plot([llcenter[0]], [llcenter[1]], 'ro', ms=10, zorder=9)
		#ax.plot([llepicen[0]], [llepicen[1]], 'm*', ms=18, zorder=10, label='Mainshock')
		#ax.plot([llfs[0]], [llfs[1]], 'c*', ms=13, zorder=10, label='Foreshock')
		#ax.legend(loc=0, numpoints=1
		#
		newratios=[]
		tickDates=[]
		tickLabels=[]
		if i==1:
			myax=myaxes[4]
			dt1=dtm.datetime(2009,11,1, tzinfo=pytz.timezone('UTC'))
			dt2=dtm.datetime(2010,5,1, tzinfo=pytz.timezone('UTC'))			
			for rw in ratios:				
				if rw[1]>dt1 and rw[1]<dt2: newratios+=[rw]
		#		
		if i==2:
			myax=None
			# do nothing
			#for rw in ratios:
			#	if rw[1]>dtm.datetime(2009,11,1, tzinfo=pytz.timezone('UTC')) and rw[1]<dtm.datetime(2010,5,1, tzinfo=pytz.timezone('UTC')):
			#		newratios+=[rw]
			#	myax=myaxes[5]
			#	myax.set_xticks([])
			#	myax.set_yticks([])
			#	myax.set_xticklabels([])
			#	myax.set_yticklabels([])
				#
		#
		if i==3:
			myax=myaxes[6]
			dt1=dtm.datetime(2010,2,1, tzinfo=pytz.timezone('UTC'))
			dt2=dtm.datetime(2010,5,1, tzinfo=pytz.timezone('UTC'))
			for rw in ratios:
				if rw[1]>dt1 and rw[1]<dt2: newratios+=[rw]
		#	
		if myax!=None:
			myax.set_ylim(auto=True)
			r2=objCat.plotIntervalRatiosAx(winlen=None, cat=circ_cat, hitThreshold=1.0, bigmag=9.0, thisAx=myax, ratios=newratios, delta_t=1, avlen=1, mainEv=mainshock, logZ=None, rbLegLoc=0, reverse=False)
			myax.plot([mainshock[0]], [1.15], 'rd', ms=7, alpha=.8, zorder=11)
			myax.plot([mainshock[0]], [1.15], 'kd', ms=9, alpha=.7, zorder=10)
			#
			miny=min(map(operator.itemgetter(5), newratios))
			maxy=max(map(operator.itemgetter(5), newratios))
			myax.set_ybound(lower=.9*miny, upper=maxy*1.1)
			myax.set_xbound(lower=dt1, upper=dtm.datetime(2010,5,1, tzinfo=pytz.timezone('UTC')))
			#
			tickDates  += [dt1]
			tickLabels += ['%d/%d' % (dt1.month, dt1.year)]
			#
			while tickDates[-1]<=dt2:
				newmonth=1+(tickDates[-1].month)%12
				#if newmonth==0: newmonth=12
				dYear = tickDates[-1].year-tickDates[0].year + (tickDates[-1].month)/12
				#
				#print "new date: ", dt1.year+dYear,newmonth,1
				tickDates  += [dtm.datetime(dt1.year+dYear,newmonth,1, tzinfo=pytz.timezone('UTC'))]
				tickLabels += ['%d/%d' % (tickDates[-1].month, tickDates[-1].year)]
			#myax.set_xticks(map(mpd.date2num, tickDates))
			myax.set_xticks(tickDates)
			myax.set_xticklabels(tickLabels)

		#
		i+=1
	 
		
	return None

def EMCQuads(targmag=7.2, mc=2.5, rfactor=.5, lats=[31.5, 33.0], lons=[-116.15, -114.5], weighted_ratios=False):
	# rfactor=.5, mc=3.0 looks great.
	Nsequence = getNofm(targmag, mc)
	emcEvent = [dtm.datetime(2010, 4, 4, 22, 40, 42, 360004, tzinfo=mhp.pytz.timezone('UTC')), 32.2862, -115.2953, 7.2, 10.0]
	emcFS    = [dtm.datetime(2009, 12, 30, 18, 48, 57, 330004, tzinfo=mhp.pytz.timezone('UTC')), 32.464, -115.1892, 5.8, 6.0]	#"foreshock"?
	a=mhp.elmayorhm(mc=mc, todt=mhp.dtm.datetime.now(mhp.pytz.timezone('UTC')), ndithers=10, nContours=nContours1, bigmag=6.0, refreshcat=True, dt0=mhp.dtm.datetime(1990,1,1, tzinfo=mhp.pytz.timezone('UTC')), lats=lats, lons=lons)
	c1=a.getcat()
	while len(c1.subcats)>0:c1.subcats.pop()
	#targmag=6.96	# this is the estimated thickness of the seismogenic zone.
	#
	## and let's strip out the far-field events (to make a prettier plot):
	#i=0
	#while i<len(c1.cat):
	#	#if rw[1]>33.25 or rw[1]<31.25 or rw[2]>-114.5 or rw[2]<-116.5: continue
	#	while rw[1]>33.25 or rw[1]<31.25 or rw[2]>-114.5 or rw[2]<-116.5:
	#
	#
	N0=int(getNofm(m=targmag, mc=mc, mt=7.6, dmstar=1.0, dms=1.0))
	Lr = 10.0**(targmag/2.0 - 1.76)
	thisr=rfactor
	print "Rupture Len: %f" % (thisr*Lr)
	#
	thetadeg=40.
	theta = thetadeg*2.0*math.pi/360.0
	mywinlen = winlen(targmag, mc, mt=7.6, doInt=False)
	
	dr=.1
	dx=dr*math.cos(theta)
	dy=dr*math.sin(theta)
	i=-6
	#i=5
	plt.figure(9)
	plt.clf()
	c1.rbomoriQuadPlot(mc=mc, targmag=targmag, bigmag=8.0, catnum=0, fignum=9, plotevents=False, intlist=[int(N0*.5), int(N0*.666), N0, int(N0*1.5)], weighted=weighted_ratios)
	plt.title('Full catalog (for reference)')
	fnum=10
	center =emcEvent[2], emcEvent[1]
	#while i<=6:
	for i in [-3,0,4]:
		#continue
		thisCenter = [center[0]-i*dx, center[1]+i*dy]
		c1.subcats+=[['circ%f' % i, circularcat(incat=c1.getcat(), latlon=[thisCenter[1], thisCenter[0]], Rkm=thisr*Lr)]]
		print "Center %d: %f, %f" % (i, thisCenter[0], thisCenter[1])
		f=plt.figure(fnum)
		plt.clf()
		#
		print "catlen: %d" % len(c1.getcat(len(c1.subcats)))
		#
		if mywinlen>(len(c1.getcat(len(c1.subcats)))*1.1): 
			i+=1.0
			fnum+=1
			continue
		c1.rbomoriQuadPlot(mc=mc, targmag=targmag, bigmag=8.0, catnum=len(c1.subcats), fignum=fnum, plotevents=False, intlist=[int(N0*.5), int(N0*.666), N0, int(N0*1.5)], thislw=2.25, weighted=weighted_ratios)
		#
		#
		f=plt.figure(fnum)
		thisAxes = f.get_axes()
		if len(thisAxes)<4: 
			i+=1.0
			fnum+=1
			continue
		#
		ax=thisAxes[2]
		ax.plot([emcEvent[0]], [1.], 'r^', ms=12, zorder=17)
		#
		ax=thisAxes[3]
		#
		circXY = getLLcircle(thisCenter[0], thisCenter[1], r=thisr*Lr*1000., rtype=0, N=500)
		Xcirc, Ycirc=c1.catmap(scipy.array(circXY[0]), scipy.array(circXY[1]))
		ax.plot(Xcirc, Ycirc, 'r--', zorder=8, lw=3, alpha=.8)
		llcenter = c1.catmap(thisCenter[0], thisCenter[1])
		llepicen = c1.catmap(emcEvent[2], emcEvent[1])
		llfs     = c1.catmap(emcFS[2], emcFS[1])
		ax.plot([llcenter[0]], [llcenter[1]], 'ro', ms=10, zorder=9)
		ax.plot([llepicen[0]], [llepicen[1]], 'm*', ms=25, zorder=10, label='Mainshock')
		ax.plot([llfs[0]], [llfs[1]], 'c*', ms=20, zorder=10, label='Foreshock')
		ax.legend(loc=0, numpoints=1 )
		# earthquakes:
		befores=[]
		afters=[]
		#for rw in c1.getcat(len(c1.subcats)):
		for rw in c1.getcat(0):
			#
			# let's narrow the map for a prettier plot.
			#if rw[1]>33.25 or rw[1]<31.25 or rw[2]>-114.5 or rw[2]<-116.5: continue
			if rw[0]<emcEvent[0]: befores+=[[rw[2], rw[1]]]
			if rw[0]>=emcEvent[0]: afters+=[[rw[2], rw[1]]]
		befores=zip(*befores)
		afters=zip(*afters)
		X,Y = c1.catmap(befores[0], befores[1])
		ax.plot(X,Y, 'b.', ms=3, zorder=8, alpha=.7, label='befores')
		X,Y = c1.catmap(afters[0], afters[1])
		ax.plot(X,Y, 'g.', ms=3, zorder=6, alpha=.7, label='after')
		ax.legend(loc=0, numpoints=1)
		#
		# now, show the mainshock and foreshock on the magnitude timeseries:
		ax=thisAxes[0]
		ax.plot(emcEvent[0], emcEvent[3], 'm*', ms=12, label='Mainshock: m=7.2, 4 Apr. 2010')
		ax.plot(emcFS[0], emcFS[3], 'c*', ms=12, label='Foreshock: m=5.8, 30 Dec. 2009')
		ax.legend(loc='upper left', numpoints=1)
		#
		#
		i+=1.0
		fnum+=1
	#
	# now, one elliptical catalog:
	if 1==1:
		xoffset=0.0
		yoffset=0.0
		thisCenter = [-115.2953+xoffset, 32.2862+yoffset]
		thistheta=-48.0
		# ellipcat(incat, latlon, a=10.0, b=5.0, ellipTheta=0.0):
		a=60.0	# in meters...
		b=20.0
		b=a/3.0
		# let's try an equal area transform with L=L_r?
		r0 = .5*10.0**(7.2/2.0 - 1.76)
		eps = 2.67		# as per Shcherbakof, parkfield (.4/.15)
		epsSqrt = math.sqrt(eps)
		a=r0*epsSqrt
		b=r0/epsSqrt
		# 
		# also consider Lr=120, which i think is a common measured value of the rupture length.
		#		
		#c1.subcats+=[['ellip%d' % i, ellipcat(incat=c1.getcat(0), latlon=[thisCenter[1], thisCenter[0]], a=a, b=b, ellipTheta=thistheta)]]
		c1.subcats+=[['ellip%d' % i, ellipcat(incat=c1.getcat(0), latlon=[thisCenter[1], thisCenter[0]], a=a, b=b, ellipTheta=thistheta)]]
		#plt.figure(47)
		#plt.clf()
		#plt.plot(map(operator.itemgetter(2), c1.getcat(len(c1.subcats))), map(operator.itemgetter(1), c1.getcat(len(c1.subcats))), 'b.')
		#plt.plot(map(operator.itemgetter(2), elcat), map(operator.itemgetter(1), elcat), 'g.')
		
		f=plt.figure(fnum)
		plt.clf()
		#
		print "catlen: %d" % len(c1.getcat(len(c1.subcats)))
		mywinlen = winlen(targmag, mc, mt=7.6, doInt=False)
		if mywinlen>(len(c1.getcat(len(c1.subcats)))*1.1): 
			i+=1.0
			fnum+=1
			#continue
		#N0=int(getNofm(m=targmag, mc=mc, mt=7.6, dmstar=1.0, dms=1.0))
		c1.rbomoriQuadPlot(mc=mc, targmag=targmag, bigmag=8.0, catnum=len(c1.subcats), fignum=fnum, plotevents=False, intlist=[int(N0*.5), int(N0*.666), N0, int(N0*1.5)], thislw=2.25, deltaLat=.5, deltaLon=.5, weighted=weighted_ratios)
		#
		#
		f=plt.figure(fnum)
		thisAxes = f.get_axes()
		if len(thisAxes)<4: 
			i+=1.0
			fnum+=1
			#continue
		#
		ax=thisAxes[2]
		ax.plot([emcEvent[0]], [1.], 'r^', ms=12, zorder=17)
		arrow_width=30
		ax.arrow(emcEvent[0]-dtm.timedelta(days=30), 4., 0., -2.5, width=arrow_width, head_width=3.*arrow_width, head_length=.1, color='m', zorder=11, alpha=.9) 	#
		#
		ax=thisAxes[3]
		#
		#circXY = getLLcircle(thisCenter[0], thisCenter[1], r=thisr*Lr*1000., rtype=0, N=500)
		circXY =  getEllipser0(x0=thisCenter[0], y0=thisCenter[1], r0=None, epsilon=None, thetaEllipse=thistheta, N=1000, a=a, b=b)
		#
		Xcirc, Ycirc=c1.catmap(scipy.array(circXY[0]), scipy.array(circXY[1]))
		ax.plot(Xcirc, Ycirc, 'b--', zorder=8, lw=3, alpha=.8)
		ax.fill(Xcirc, Ycirc, 'b', zorder=4, alpha=.2)
		llcenter = c1.catmap(thisCenter[0], thisCenter[1])
		llepicen = c1.catmap(emcEvent[2], emcEvent[1])
		llfs     = c1.catmap(emcFS[2], emcFS[1])
		#ax.plot([llcenter[0]], [llcenter[1]], 'ro', ms=10, zorder=9)
		ax.plot([llepicen[0]], [llepicen[1]], 'm*', ms=25, zorder=10, label='Mainshock')
		ax.plot([llfs[0]], [llfs[1]], 'c*', ms=20, zorder=10, label='Foreshock')
		ax.legend(loc=0, numpoints=1 )
		# earthquakes:
		befores=[]
		afters=[]
		#for rw in c1.getcat(len(c1.subcats)):
		ins=[[],[]]
		for rw in c1.getcat(0):
			#
			if rw[0]<emcEvent[0]: befores+=[[rw[2], rw[1]]]
			if rw[0]>=emcEvent[0]: afters+=[[rw[2], rw[1]]]
		#for rw in c1.getcat(len(c1.subcats)):
		#	ins[0]+=[rw[2]]
		#	ins[1]+=[rw[1]]
			
		befores=zip(*befores)
		afters=zip(*afters)
		X,Y = c1.catmap(befores[0], befores[1])
		ax.plot(X,Y, 'b.', ms=3, zorder=8, alpha=.7, label='befores')
		X,Y = c1.catmap(afters[0], afters[1])
		ax.plot(X,Y, 'g.', ms=3, zorder=6, alpha=.7, label='after')
		#xin, yin = c1.catmap(ins[0], ins[1])
		#ax.plot(xin, yin, 'r.', ms=3, zorder=9, alpha=.5)
		#
		ax.legend(loc=0, numpoints=1)
		#
		# foreshock, mainshock mags.
		ax=thisAxes[0]
		ax.plot(emcEvent[0], emcEvent[3], 'm*', ms=12, label='Mainshock: m=7.2, 4 Apr. 2010')
		ax.plot(emcFS[0], emcFS[3], 'c*', ms=12, label='Foreshock: m=5.8, 30 Dec. 2009')
		ax.legend(loc='upper left', numpoints=1)
		#
		#
		i+=1.0
		fnum+=1	
	return c1

def parkfieldQuads2(targmag=5.96, mc=1.5, rfactor=.5, weighted_ratios=False):
	# array of elliptical Parkfield RBTS. (production plot comes from here).
	Nsequence = getNofm(targmag, mc)
	intervalList = [int(Nsequence*.5), int(Nsequence*.666), Nsequence, int(Nsequence*1.5)]
	pfEvent = [dtm.datetime(2004, 9, 28, 22, 40, 0, 0, tzinfo=mhp.pytz.timezone('UTC')), 35.8182,  -120.366,  5.97,  8.58]
	a=mhp.parkfieldhm(mc=mc, todt=mhp.dtm.datetime.now(mhp.pytz.timezone('UTC')), ndithers=10, nContours=nContours1, bigmag=5.0, refreshcat=True, dt0=mhp.dtm.datetime(1990,1,1, tzinfo=mhp.pytz.timezone('UTC')), lats=[35.415, 36.5], lons=[-120.9, -119.97])
	#
	c1=a.getcat()
	while len(c1.subcats)>0:c1.subcats.pop()
	#
	Lr = 10.0**(targmag/2.0 - 1.76)
	thisr=rfactor
	#
	thetadeg=42.
	theta = thetadeg*2.0*math.pi/360.0
	dr=.05
	dx=dr*math.cos(theta)
	dy=dr*math.sin(theta)
	i=-2.0
	fnum=10
	center =pfEvent[2], pfEvent[1]
	# epicenter and also see: center: 35.885113, -120.440314 for good plots.
	while i<=9:
		thisCenter = [center[0]-i*dx, center[1]+i*dy]
		if i==9: thisCenter = [-120.5, 35.9]
		#c1.subcats+=[['circ%f' % i, circularcat(incat=c1.getcat(), latlon=[thisCenter[1], thisCenter[0]], Rkm=thisr*Lr)]]
		print "fnum: %d, center: %f, %f" % (fnum, thisCenter[1], thisCenter[0])
		c1.addEllipCat('ellip%f' % i, c1.cat, thetadeg, thisCenter[1], thisCenter[0], .4, .15)
		f=plt.figure(fnum)
		c1.rbomoriQuadPlot(mc=mc, targmag=targmag, bigmag=6.0, catnum=len(c1.subcats), fignum=fnum, plotevents=False, intlist=intervalList, thislw=2.25, deltaLat=.25, deltaLon=.25, weighted=weighted_ratios)
		mywinlen = winlen(targmag, mc, mt=7.6, doInt=False)
		if mywinlen>(len(c1.getcat(len(c1.subcats)))*1.1): 
			i+=1.0
			fnum+=1
			continue
		#
		print "catlen: %d" % len(c1.getcat(len(c1.subcats)))
		#
		f=plt.figure(fnum)
		thisAxes = f.get_axes()
		if len(thisAxes)<4: 
			i+=1.0
			fnum+=1
			continue
		#
		ax=thisAxes[2]
		ax.plot([pfEvent[0]], [1.], 'r^', ms=12, zorder=17)
		arrow_width=40.
		ax.arrow(pfEvent[0]-dtm.timedelta(days=20), 4., 0., -1.5, width=arrow_width, head_width=3.*arrow_width, head_length=.3, color='m', zorder=11, alpha=.9) 	#
		#
		ax=thisAxes[3]
		#
		#circXY = getLLcircle(thisCenter[0], thisCenter[1], r=thisr*Lr*1000., rtype=0, N=500)
		#Xcirc, Ycirc=c1.catmap(scipy.array(circXY[0]), scipy.array(circXY[1]))
		#ax.plot(Xcirc, Ycirc, 'm--', zorder=8, lw=2, alpha=.8)
		#
		# catalog ellipse:
		r0 = .5*10.0**(targmag/2.0 - 1.76)
		eps = math.sqrt(2.67)		# as per Shcherbakof, parkfield (.4/.15)
		#a=r0*eps
		a=40.
		b=15.
		#b=r0/(eps**2.)
		circXY =  getEllipser0(x0=thisCenter[0], y0=thisCenter[1], r0=None, epsilon=None, thetaEllipse=-thetadeg, N=500, a=a, b=b)
		Xcirc, Ycirc=c1.catmap(scipy.array(circXY[0]), scipy.array(circXY[1]))
		ax.plot(Xcirc, Ycirc, 'b--', zorder=8, lw=3, alpha=.8)
		ax.fill(Xcirc, Ycirc, 'b', zorder=3, alpha=.25)
		#
		llcenter = c1.catmap(thisCenter[0], thisCenter[1])
		llepicen = c1.catmap(pfEvent[2], pfEvent[1])
		ax.plot([llcenter[0]], [llcenter[1]], 'ro', ms=10, alpha=.8, zorder=10)
		ax.plot([llepicen[0]], [llepicen[1]], 'm*', ms=25, alpha=.8, zorder=9, label='Mainshock')
		ax.legend(loc=0, numpoints=1 )
		# earthquakes:
		befores=[]
		afters=[]
		#for rw in c1.getcat(len(c1.subcats)):
		for rw in c1.getcat(0):
			if rw[0]<pfEvent[0]: befores+=[[rw[2], rw[1]]]
			if rw[0]>=pfEvent[0]: afters+=[[rw[2], rw[1]]]
		befores=zip(*befores)
		afters=zip(*afters)
		X,Y = c1.catmap(befores[0], befores[1])
		ax.plot(X,Y, 'b.', ms=3, zorder=5, alpha=.7, label='befores')
		X,Y = c1.catmap(afters[0], afters[1])
		ax.plot(X,Y, 'g.', ms=3, zorder=4, alpha=.7, label='afters')
		ax.legend(loc=0, numpoints=1)
		'''
		X=map(operator.itemgetter(2), c1.getcat(len(c1.subcats)))
		Y=map(operator.itemgetter(1), c1.getcat(len(c1.subcats)))
		X,Y=c1.catmap(X,Y)
		ax.plot(X,Y, 'r.', ms=3, zorder=9, alpha=.24)
		'''
		#
		i+=1.0
		fnum+=1
	return c1

def parkfieldQuadsCircular(targmag=5.96, mc=1.5, rfactor=.5, weighted_ratios=False):
	Nsequence = getNofm(targmag, mc)
	intervalList = [int(Nsequence*.5), int(Nsequence*.666), Nsequence, int(Nsequence*1.5)]
	pfEvent = [dtm.datetime(2004, 9, 28, 22, 40, 0, 0, tzinfo=mhp.pytz.timezone('UTC')), 35.8182,  -120.366,  5.97,  8.58]
	a=mhp.parkfieldhm(mc=mc, todt=mhp.dtm.datetime.now(mhp.pytz.timezone('UTC')), ndithers=10, nContours=nContours1, bigmag=5.0, refreshcat=True, dt0=mhp.dtm.datetime(1990,1,1, tzinfo=mhp.pytz.timezone('UTC')))
	c1=a.getcat()
	while len(c1.subcats)>0:c1.subcats.pop()
	#
	Lr = 10.0**(targmag/2.0 - 1.76)
	thisr=rfactor
	#
	thetadeg=42.
	theta = thetadeg*2.0*math.pi/360.0
	#dr=.05
	dr=.01
	dx=dr*math.cos(theta)
	dy=dr*math.sin(theta)
	#
	a=40.
	b=15.
	center =pfEvent[2], pfEvent[1]
	ellipseXY =  getEllipser0(x0=center[0], y0=center[1], r0=None, epsilon=None, thetaEllipse=-thetadeg, N=500, a=a, b=b)
	Xellipse, Yellipse=c1.catmap(scipy.array(ellipseXY[0]), scipy.array(ellipseXY[1]))
	#
	i=-2.0
	fnum=10
	i=7
	
	#while i<8:
	# (in this configuration, pane 22 is the sweet spot).
	while (center[1]+i*dy)<36.0:
		thisCenter = [center[0]-i*dx, center[1]+i*dy]
		#c1.subcats+=[['circ%f' % i, circularcat(incat=c1.getcat(), latlon=[thisCenter[1], thisCenter[0]], Rkm=thisr*Lr)]]
		#c1.addEllipCat('ellip%f' % i, c1.cat, thetadeg, thisCenter[1], thisCenter[0], .4, .15)
		#
		if fnum!=22:
			i+=1
			fnum+=1
			continue
		#
		print "center: ", thisCenter
		c1.subcats+=[['circ%f' % i, circularcat(incat=c1.getcat(0), latlon=[thisCenter[1], thisCenter[0]], Rkm=thisr*Lr)]]
		f=plt.figure(fnum)
		c1.rbomoriQuadPlot(mc=mc, targmag=targmag, bigmag=6.0, catnum=len(c1.subcats), fignum=fnum, plotevents=False, intlist=intervalList, thislw=2.25, deltaLat=.25, deltaLon=.25, weighted=weighted_ratios)
		mywinlen = winlen(targmag, mc, mt=7.6, doInt=False)
		if mywinlen>(len(c1.getcat(len(c1.subcats)))*1.1): 
			i+=1.0
			fnum+=1
			continue
		#
		print "catlen: %d" % len(c1.getcat(len(c1.subcats)))
		#
		f=plt.figure(fnum)
		thisAxes = f.get_axes()
		if len(thisAxes)<4: 
			i+=1.0
			fnum+=1
			continue
		#
		ax = thisAxes[0]
		ax.plot([pfEvent[0]], [pfEvent[3]], 'm*', ms=10, label='mainshock')
		ax.legend(loc='upper left', numpoints=1)
		#
		ax=thisAxes[2]
		#ax.plot([pfEvent[0]], [1.1], 'm*', ms=15, zorder=17)
		ax.plot([pfEvent[0], pfEvent[0]], [.3,7.], 'm--', lw=3, alpha=.8, zorder=11)
		#
		ax=thisAxes[3]
		#
		circXY = getLLcircle(thisCenter[0], thisCenter[1], r=thisr*Lr*1000., rtype=0, N=500)
		Xcirc, Ycirc=c1.catmap(scipy.array(circXY[0]), scipy.array(circXY[1]))
		ax.plot(Xcirc, Ycirc, 'r--', zorder=8, lw=3, alpha=.8)
		
		circXY = getLLcircle(thisCenter[0], thisCenter[1], r=Lr*1000., rtype=0, N=500)
		Xcirc, Ycirc=c1.catmap(scipy.array(circXY[0]), scipy.array(circXY[1]))
		ax.plot(Xcirc, Ycirc, 'r-.', zorder=8, lw=3, alpha=.8)
		#
		## catalog ellipse:
		#r0 = .5*10.0**(targmag/2.0 - 1.76)
		#eps = math.sqrt(2.67)		# as per Shcherbakof, parkfield (.4/.15)
		#a=r0*eps
		#b=r0/eps
		ax.plot(Xellipse, Yellipse, 'b--', zorder=8, lw=2, alpha=.8)
		ax.fill(Xellipse, Yellipse, 'b', zorder=3, alpha=.25)
		#
		llcenter = c1.catmap(thisCenter[0], thisCenter[1])
		llepicen = c1.catmap(pfEvent[2], pfEvent[1])
		ax.plot([llcenter[0]], [llcenter[1]], 'ro', ms=10, zorder=9)
		ax.plot([llepicen[0]], [llepicen[1]], 'm*', ms=24, zorder=10, label='Mainshock')
		ax.legend(loc=0, numpoints=1 )
		# earthquakes:
		befores=[]
		afters=[]
		#for rw in c1.getcat(len(c1.subcats)):
		for rw in c1.getcat(0):
			if rw[0]<pfEvent[0]: befores+=[[rw[2], rw[1]]]
			if rw[0]>=pfEvent[0]: afters+=[[rw[2], rw[1]]]
		befores=zip(*befores)
		afters=zip(*afters)
		X,Y = c1.catmap(befores[0], befores[1])
		ax.plot(X,Y, 'b.', ms=3, zorder=8, alpha=.7, label='befores')
		X,Y = c1.catmap(afters[0], afters[1])
		ax.plot(X,Y, 'g.', ms=3, zorder=6, alpha=.7, label='afters')
		ax.legend(loc='lower center', numpoints=1)
		'''
		X=map(operator.itemgetter(2), c1.getcat(len(c1.subcats)))
		Y=map(operator.itemgetter(1), c1.getcat(len(c1.subcats)))
		X,Y=c1.catmap(X,Y)
		ax.plot(X,Y, 'r.', ms=3, zorder=9, alpha=.24)
		'''
		#
		i+=1.0
		fnum+=1
	return c1

def parkfieldQuads(targmag=5.96, mc=1.5, weighted_ratios=False):
	# testing differenc circle radii. this approach was not used much since it did not produce any (apparently) significant results.
	#
	Nsequence = getNofm(targmag, mc)
	pfEvent = [dtm.datetime(2004, 9, 28, 17, 15, 24, 249999, tzinfo=mhp.pytz.timezone('UTC')), 35.8182,  -120.366,  5.97,  8.58]
	a=mhp.parkfieldhm(mc=mc, todt=mhp.dtm.datetime.now(mhp.pytz.timezone('UTC')), ndithers=10, nContours=nContours1, bigmag=5.0, refreshcat=True, dt0=mhp.dtm.datetime(1990,1,1, tzinfo=mhp.pytz.timezone('UTC')))
	c1=a.getcat()
	while len(c1.subcats)>0:c1.subcats.pop()
	#
	rfactors = [.5,.75, 1.0, 1.5, 2.0, 2.5, 3.0]
	drawCircles=[]
	Lr = 10.0**(targmag/2.0 - 1.76)
	#
	myIndex=1
	for thisr in rfactors:
		c1.subcats+=[['circ%f' % thisr, circularcat(incat=c1.getcat(), latlon=[pfEvent[1], pfEvent[2]], Rkm=thisr*Lr)]]
		#
		fnum=myIndex+9
		c1.rbomoriQuadPlot(mc=mc, targmag=targmag, bigmag=5.0, catnum=myIndex, fignum=fnum, plotevents=False, intlist=[250, 500, 770, 1000], thislw=2.25, weighted=weighted_ratios)
		#
		f=plt.figure(fnum)
		thisAxes = f.get_axes()
		#
		ax=thisAxes[2]
		ax.plot([pfEvent[0]], [1.0], 'r^', ms=15, alpha=11)
		#
		ax=thisAxes[3]
		#
		drawCircles+=[[len(c1.subcats), thisr]]
		circXY = getLLcircle(pfEvent[2], pfEvent[1], r=thisr*Lr*1000., rtype=0, N=500)
		Xcirc, Ycirc=c1.catmap(scipy.array(circXY[0]), scipy.array(circXY[1]))
		ax.plot(Xcirc, Ycirc, 'r--', zorder=8, lw=2, alpha=.8)
		#
		circXY = getLLcircle(pfEvent[2], pfEvent[1], r=1.0*Lr*1000., rtype=0, N=500)
		XcircLr, YcircLr=c1.catmap(scipy.array(circXY[0]), scipy.array(circXY[1]))
		ax.plot(XcircLr, YcircLr, 'k--', lw=2, alpha=.8, zorder=7)
		#
		xEvents = map(operator.itemgetter(2), c1.getcat(myIndex))
		yEvents = map(operator.itemgetter(1), c1.getcat(myIndex))
		X,Y = c1.catmap(xEvents, yEvents)
		ax.plot(X,Y, 'b,')
		#
		myIndex+=1
	return c1

def chilequads(targmag=8.8, mc=5.0, weighted_ratios=False):
	Nsequence = getNofm(targmag, mc)
	a=mhp.chilehm(mc=mc, todt=mhp.dtm.datetime.now(mhp.pytz.timezone('UTC')), ndithers=10, nContours=nContours1, bigmag=7.0, lons=[135., 148.5], lats=[30., 45.25], refreshcat=True, dt0=mhp.dtm.datetime(1990,1,1, tzinfo=mhp.pytz.timezone('UTC')))
	c1=a.getcat()
	
	print "nSubcats: ", len(c1.subcats)
	while len(c1.subcats): c1.subcats.pop()
	print "nSubcats: ", len(c1.subcats)
	event=c1.getMainEvent()
	LatLon = event[1], event[2]
	
	Lr = 10.0**(.5*targmag - 1.76)
	print "Lr = %f" % Lr
	for thisr in rfactors:
		c1.subcats+=[['circ%f' % thisr, circularcat(incat=c1.getcat(), latlon=LatLon, Rkm=thisr*Lr)]]
		drawCircles+=[[len(c1.subcats), thisr]]
	#
	print "nSubcats: ", len(c1.subcats)
	for k in xrange(len(c1.subcats)):
		c1.subcats[k][1].sort(key=lambda x: x[0])
	# ... and finish eventually...

def chichiquads(targmag=7.3, mc=2.5, lfactor=.6, catname='cats/taiwan1994.cat', weighted_ratios=False):
	# 7.1, or ANSS at 7.3, 7.6? maybe it's up to us to find the "correct" magnitude?
	# m=7.3 and lfactor=.6 produce a really nice plot. The corresponding m(L)=7.46 does not appear to produce
	# anything useful.
	plt.ion()
	#
	c1=mhp.eqp.eqcatalog()
	c1.loadCatFromFile(catname, minmag=mc)
	c1.cat.sort(key = lambda x:x[0])
	#c1.loadCatFromFile('cats/taiwan1994.cat', minmag=mc)
	#1999, 9,
	chichidate = dtm.datetime( 1999, 9, 20, 17, 47, 15, 850000, tzinfo=pytz.timezone('UTC'))
	chichi = [chichidate, 23.8525, 120.8155, 7.3, 8.0]
	#
	Lr = 10.0**(targmag/2.0 - 1.76)
	#
	# now, make a circular subcat around chichi.
	while len(c1.subcats)>0: c1.subcats.pop()
	#
	# mainshock centric:
	c1.subcats+=[['r0', circularcat(c1.getcat(0), latlon=[chichi[1], chichi[2]], Rkm=Lr*lfactor)]]
	#
	# and a subcat centered around the CM of the aftershocks (approximately):
	#pre CM  (x,y): 121.580116, 23.625129
	#post CM (x,y): 120.984787, 23.882699
	aftershocksCM=[23.883, 120.985]
	c1.subcats+=[['r1', circularcat(c1.getcat(0), latlon=aftershocksCM, Rkm=Lr*lfactor)]]
	#
	'''
	ratios = c1.rbomoriQuadPlot(mc=mc, targmag=targmag, bigmag=targmag-1.1, catnum=1, fignum=0, plotevents=False, intlist=None, thislw=2.25)
	thisAxes = plt.figure(0).get_axes()
	ax=thisAxes[3]
	plotEvents(c1, targmag, mc, ax, chichidate)
	circLR = getLLcircle(chichi[2], chichi[1], r=lfactor*Lr*1000., rtype=0, N=500)
	XcircLr, YcircLr=c1.catmap(scipy.array(circLR[0]), scipy.array(circLR[1]))
	c1.catmap.plot(XcircLr, YcircLr, '-')
	'''
	#
	# and aftershock CM:
	c1.rbomoriQuadPlot(mc=mc, targmag=targmag, bigmag=9.5, catnum=2, fignum=2, plotevents=False, intlist=None, mapLLlat=22.5, mapLLlon=119.75, mapURlat=25.0, mapURlon=122.25, thislw=2.25, weighted=weighted_ratios)
	thisAxes = plt.figure(2).get_axes()
	ax=thisAxes[3]
	# plot relevant catalog and "dressing":
	#plotEvents(c1, targmag, mc, ax, chichidate)
	# plot all events:
	Xall, Yall, Xcirc, Ycirc=[],[],[],[]
	#Xall, Yall = c1.catmap(map(operator.itemgetter(2), c1.getcat(0)), map(operator.itemgetter(1), c1.getcat(0)))
	#Xcirc, Ycirc = c1.catmap(map(operator.itemgetter(2), c1.getcat(2)), map(operator.itemgetter(1), c1.getcat(2)))
	# we need to disciriminate by  mc, so we have to go old-school and do a loop:
	thismc=3.5
	my_bigmag=6.2
	symb='d'
	bigmagsXY = []
	for rw in c1.getcat(0):
		if rw[3]>thismc:
			Xall+=[rw[2]]
			Yall+=[rw[1]]
	for rw in c1.getcat(2):
		if rw[3]>thismc:
			Xcirc+=[rw[2]]
			Ycirc+=[rw[1]]
		if rw[3]>=my_bigmag:
			if rw[0]>chichi[0]: symb='s'
			if rw[0]<chichi[0]: symb='d'
			if rw[0]==chichi[0]: continue	# forego the mainshock; we'll plot it separately.
			x,y=c1.catmap(rw[2], rw[1])
			c1.catmap.plot([x], [y], symb, ms=20*(rw[3]/chichi[3]), alpha=.85, zorder=7, label='m=%.1f, %d/%d/%d' % (rw[3], rw[0].day, rw[0].month, rw[0].year))
	Xall, Yall = c1.catmap(Xall, Yall)
	Xcirc, Ycirc = c1.catmap(Xcirc, Ycirc)
	#
	c1.catmap.plot(Xall, Yall, 'g.', ms=3, alpha=.6, zorder=2, ax=ax)
	c1.catmap.plot(Xcirc, Ycirc, 'b.', ms=3, alpha=.6, zorder=3, ax=ax)
	xchichi, ychichi = c1.catmap(chichi[2], chichi[1])
	c1.catmap.plot([xchichi], [ychichi], 'k*', zorder=12, alpha=.8, ax=ax, ms=23)
	c1.catmap.plot([xchichi], [ychichi], 'r*', zorder=13, alpha=.8, ax=ax, ms=18, label='m=7.3, %d/%d/%d' % (chichi[0].day, chichi[0].month, chichi[0].year))
	#
	circLR = getLLcircle(aftershocksCM[1], aftershocksCM[0], r=lfactor*Lr*1000., rtype=0, N=500)
	XcircLr, YcircLr=c1.catmap(scipy.array(circLR[0]), scipy.array(circLR[1]))
	c1.catmap.plot(XcircLr, YcircLr, 'r--', lw=2, zorder=11)
	xcm, ycm=c1.catmap(aftershocksCM[1], aftershocksCM[0])
	ax.legend(loc=0, numpoints=1)
	#c1.catmap.plot([xcm], [ycm], 'ro', ms=10, zorder=12, alpha=.8, ax=ax, label='Center')
	'''
	c1.catmap.llcrnrlon=120.
	c1.catmap.llcrnrlat=22.75
	c1.catmap.urcrnrlon=122.
	c1.catmap.urcrnrlat=24.75
	c1.catmap.set_axes_limits(ax=ax)
	plt.legend(loc=0, numpoints=1)
	'''
	#
	ax=thisAxes[2]
	ax.set_xlim([dtm.datetime(1998,6,1, tzinfo=pytz.timezone('UTC')), dtm.datetime(2008,2,1, tzinfo=pytz.timezone('UTC'))])
	ax.set_ylim(.45, 6.3)
	#ax.plot([chichidate, chichidate], [.3, 6.], 'm-', lw=3, ms=10, alpha=.8, zorder=11)
	ax.plot([chichidate], [1.1], 'md', ms=10, alpha=.8, zorder=11)
	arrow_width=15.
	ax.arrow(chichidate, 6., 0., -3., width=arrow_width, head_width=3.*arrow_width, head_length=.3, color='m', zorder=11, alpha=.9) 	#
	ax=thisAxes[0]
	ax.plot([chichidate], [chichi[3]], 'md', lw=3, ms=12, alpha=.8, zorder=11)
	ax.set_ylim([mc-.1, chichi[3]+ .1])
	#
	# now, add a zoomed in plot of the mainshoc RBTS
	# myaxes+=[plt.axes([.1, .68, .45, .3], sharex=myaxes[0])]
	axescolor  = None #'#f6f6f6'
	thisAxes+=[plt.axes([.22, .55, .25, .2], frameon=True, axisbg=axescolor)]
	#
	myax=thisAxes[4]
	#myax.set_ylim(auto=True)
	winlen = mhp.getNsample(m=targmag, mc=mc)
	rbavelen=max(1, int(winlen/10.))
	#
	if weighted_ratios==False:
		r2=c1.plotIntervalRatiosAx(winlen=winlen, cat=c1.getcat(2), hitThreshold=1.0, bigmag=9.0, thisAx=myax, ratios=None, delta_t=1, avlen=rbavelen, mainEv=chichi, logZ=None, rbLegLoc=0, reverse=False)
	if weighted_ratios==True:
		r2=c1.plotWeightedIntervalRatiosAx(winlen=winlen, cat=c1.getcat(2), hitThreshold=1.0, bigmag=9.0, thisAx=myax, ratios=None, delta_t=1, avlen=rbavelen, mainEv=chichi, logZ=None, rbLegLoc=0, reverse=False)
	#
	myax.plot([chichi[0]], [1.15], 'md', ms=7, alpha=.8, zorder=11)
	myax.plot([chichi[0]], [1.15], 'kd', ms=9, alpha=.7, zorder=10)
	#
	#miny=min(map(operator.itemgetter(5), newratios))
	#maxy=max(map(operator.itemgetter(5), newratios))
	miny=.6
	maxy=3.
	dtstart=dtm.datetime(1999,2,1, tzinfo=pytz.timezone('UTC'))
	dtend=dtm.datetime(2000,1,1, tzinfo=pytz.timezone('UTC'))
	myax.set_ybound(lower=.9*miny, upper=maxy*1.1)
	myax.set_xbound(lower=dtstart, upper=dtend)
	#
	tickDates  = [dtstart]
	tickLabels = ['%d/%d' % (dtstart.month, dtstart.year)]
	#
	while tickDates[-1]<=dtend:
		
		newmonth=1+(1+tickDates[-1].month)%12
		#if newmonth==0: newmonth=12
		dYear = tickDates[-1].year-tickDates[0].year + (tickDates[-1].month)/12
		#
		#print "new date: ", dt1.year+dYear,newmonth,1
		tickDates  += [dtm.datetime(dtstart.year+dYear,newmonth,1, tzinfo=pytz.timezone('UTC'))]
		tickLabels += ['%d/%d' % (tickDates[-1].month, tickDates[-1].year)]
	#myax.set_xticks(map(mpd.date2num, tickDates))
	myax.set_xticks(tickDates)
	myax.set_xticklabels(tickLabels)	
	#
	myax=thisAxes[1]
	ax.legend(loc='lower right', numpoints=0)
	#
	return c1
	#return thisAxes[2]

# PAGEOPH-rbints
def sumatraQuad(mc=5.0, targmag=9.1, rbavelen=None, bigmag=9.5, intlist=None, catname='cats/sumatra.cat', refreshcat=False, plotevents=False, mt=7.55, lons=[92.0, 106.0],lats=[-9.0, 10.0], lfactor=.5, weighted_ratios=False):
	#
	# mc=4.75 might be ok, but i think we get better results for 5.0
	#
	# the full catalog produces a pretty nice signature for both 9.1 and 9.3. the 9.3 signature has that nice, quasi (log)linear
	# steady slope into the negative with a blue blip at the mainshock, and then again red carrying into a major cluster of events in 
	# late Jan. for the PAGEOPH paper, let's just stick to m=9.1 and use the full cat + 2 circular cats (maybe) as examples.
	#
	# ... except that, when we "correct" the N sacaling to use an averaged dm_bath value, we only have enough data for m=9.1 in the extended catalog.
	# the 'sweet spot' for a circular catalog is probably out there, but i'm not finding it easily, so let's just use sumatra to expand the discussion to
	# broad catalogs and transition to the prospect of a hazard map.
	#
	plt.ion()
	# for m=9.1, N=223, for m=8.6, N=70 (for mc=4.75)
	# the (a) catalog:
	if catname==None:
		refreshcat=True
		catname='cats/indonesia2012.cat'
	if refreshcat==True:
		#cl1=yp.catfromANSS(lon=[92.0, 106.0],lat=[-9.0, 10.0], minMag=3.5, dates0=[yp.dtm.datetime(1990,1,1, tzinfo=pytz.timezone('UTC')), yp.dtm.datetime(2012,6,19, tzinfo=pytz.timezone('UTC'))], fout=catname)
		cl1=atp.catfromANSS(lon=lons,lat=lats, minMag=mc, dates0=[dtm.datetime(1990,1,1, tzinfo=pytz.timezone('UTC')), dtm.datetime.now(pytz.timezone('UTC'))], fout=catname)
	#
	c1=eqp.eqcatalog([])
	c1.mt=mt
	c1.mc=mc
	c1.loadCatFromFile(catname)
	c1.cat.sort(key = lambda x:x[0])
	dlambda=1.76
	Lr=10.0**(targmag/2.0 - dlambda)
	#x0,y0 = c1.catmap(95.854, 3.316)	# mainshock
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
	#winlen=int(10**(targmag - (2.0+mc)))
	#winlen = Nofm(m=targmag, mc=mc, mt=7.6, b1=1.0, b2=1.5)*10**(-2.0)	# introducing -(dm+dm_s)
	#winlen=int(10*(int(winlen/10)))
	if rbavelen==None: avlen=mhp.getNave(m=targmag, mt=mt, mc=mc)
	winlen=int(mhp.getNsample(m=targmag, mc=mc, mt=mt, doint=True))
	#print "rbts for winlen=%d" % winlen
	#
	arb=c1.rbomoriQuadPlot(mc=mc, winlen=winlen, rbavelen=avlen, bigmag=bigmag, intlist=intlist, catnum=0, fignum=0, plotevents=plotevents, logZ=None, thislw=2.25, weighted=weighted_ratios)
	Xcat=map(operator.itemgetter(2), c1.getcat(0))
	Ycat=map(operator.itemgetter(1), c1.getcat(0))
	Xcat, Ycat = c1.catmap(Xcat, Ycat)
	#
	for i in xrange(len(foreshocks)):
		#c1.subcats+=[['r1', circularcat(c1.getcat(0), latlon=latlonMS, Rkm=Lr*lfactor)]]
		#if i>0: break		# only the full catalog is giving us much of a result.
		#
		thiscatnum=i
		thisfignum=i
		thismag=foreshocks[i][3]
		thisLr = 10.0**(thismag/2.0 - dlambda)
		#
		#if i>0:
		fsLatLon=[foreshocks[i][1], foreshocks[i][2]]
		
		if i>0: c1.subcats+=[['r%d' %i, circularcat(c1.getcat(0), latlon=fsLatLon, Rkm=Lr*lfactor)]]
		print "catalog %d, len: %d" % (thiscatnum, len(c1.getcat(thiscatnum)))
		#arb=c1.rbomoriQuadPlot(mc=mc, winlen=winlen, rbavelen=avlen, bigmag=bigmag, intlist=intlist, catnum=thiscatnum, fignum=thisfignum, plotevents=plotevents, logZ=None, thislw=2.25)
		#
		#arb=c1.rbomoriQuadPlot(mc=mc, targmag=targmag, rbavelen=None, bigmag=bigmag, intlist=intlist, catnum=thiscatnum, fignum=thisfignum, plotevents=plotevents, logZ=None, thislw=2.25)
		arb=c1.rbomoriQuadPlot(mc=mc, targmag=thismag, rbavelen=None, bigmag=bigmag, intlist=intlist, catnum=thiscatnum, fignum=thisfignum, plotevents=plotevents, logZ=None, thislw=2.25, deltaLat=1.5, deltaLon=1.5, weighted=weighted_ratios)
		if arb==None: continue
		plt.figure(thisfignum)
		thisAxes = plt.figure(thisfignum).get_axes()
		#
		ax=thisAxes[3]
		#if i>0:
		#
		circLR = getLLcircle(fsLatLon[1], fsLatLon[0], r=lfactor*Lr*1000., rtype=0, N=500)
		XcircLr, YcircLr=c1.catmap(scipy.array(circLR[0]), scipy.array(circLR[1]))
		#c1.catmap.plot(XcircLr, YcircLr, 'g-', lw=2, alpha=.8)
		ax.plot(XcircLr, YcircLr, 'r--', lw=2, alpha=.8)
		#
		x,y=c1.catmap(mainshock[2], mainshock[1])
		ax.plot([x], [y], 'm*', ms=18, alpha=.8, zorder=13, label='m=9.1, 2004/12/26')
		ax.plot([x], [y], 'k*', ms=21, alpha=.8, zorder=12)
		for fs in foreshocks:
			if fs[0]==mainshock[0]: continue
			x,y=c1.catmap(fs[2], fs[1])
			if fs[0]>mainshock[0]: marker='s'
			if fs[0]< mainshock[0]: marker='d'
			#
			thislabel=None
			if fs[3]<9.0: thislabel='m=%.1f, %d/%d/%d' % (fs[3], fs[0].year, fs[0].month, fs[0].day)
			ax.plot([x], [y], marker, ms=20*fs[3]/targmag, alpha=.8, zorder=11, label=thislabel )
		#
		littleX=map(operator.itemgetter(2), c1.getcat(thiscatnum))
		littleY=map(operator.itemgetter(1), c1.getcat(thiscatnum))
		littleX, littleY = c1.catmap(littleX, littleY)
		#ax.plot(Xcat[0:mevIndex], Ycat[0:mevIndex], 'b.', ms=3, alpha=.6, zorder=10)
		#ax.plot(Xcat[mevIndex:], Ycat[mevIndex:], 'r.', ms=3, alpha=.6, zorder=10)
		ax.plot(Xcat, Ycat, 'g.', ms=3, alpha=.6, zorder=7)
		ax.plot(littleX, littleY, 'b.', ms=3, alpha=6, zorder=8)
		plt.legend(loc=0, numpoints=1)
		#
		ax=thisAxes[2]
		ax.plot([mainshock[0]], [1.1], 'cd', ms=12, alpha=.8, zorder=11)
		ax.plot([foreshocks[i][0]], [1.1], 'md', ms=12, alpha=.8, zorder=11)
		arrow_width=30.
		if i!=2: arrow_width=45
		if i==3: arrow_width=60
		ax.arrow(foreshocks[i][0]-dtm.timedelta(days=10), 7., 0., -2.8, width=arrow_width, head_width=3.*arrow_width, head_length=.5, color='m', zorder=11, alpha=.9) 	#
		ax=thisAxes[0]
		ax.plot([mainshock[0]], [mainshock[3]], 'cd', ms=10, alpha=.8, zorder=5)
		ax.plot([foreshocks[i][0]], [foreshocks[i][3]], 'md', ms=10, alpha=.8, zorder=5)
		ax.set_ylim([mc, mainshock[3]+ .1])
	#
	# now, annotate the relevant shocks a bit.
	thisax = plt.figure(2).get_axes()
	thisax[2].set_xlim(dtm.datetime(2001, 12, 1, tzinfo=pytz.timezone('UTC')), dtm.datetime.now())
	thisax[2].set_ylim(.32, 6.3)
	#ax=thisax[2]
	#ax.dx=mpd.date2num(dtm.datetime(2004,1,1, tzinfo=pytz.timezone('UTC'))), y=3.5, dx=365., dy=-2.25, length_includes_head=True, width=.003)
	axescolor=None
	thisax+=[plt.axes([.25, .55, .27, .2], frameon=True, axisbg=axescolor)]
	ax = thisax[4]
	winlen = mhp.getNsample(m=foreshocks[2][3], mc=mc)
	rbavelen=max(1, int(winlen/10.))
	#
	if weighted_ratios==False:
		r2=c1.plotIntervalRatiosAx(winlen=winlen, cat=c1.getcat(2), hitThreshold=1.0, bigmag=9.0, thisAx=ax, ratios=None, delta_t=1, avlen=rbavelen, mainEv=foreshocks[2], logZ=None, rbLegLoc=0, reverse=False)
	if weighted_ratios==True:
		r2=c1.plotWeightedIntervalRatiosAx(winlen=winlen, cat=c1.getcat(2), hitThreshold=1.0, bigmag=9.0, thisAx=ax, ratios=None, delta_t=1, avlen=rbavelen, mainEv=foreshocks[2], logZ=None, rbLegLoc=0, reverse=False)
	#
	ax.plot([foreshocks[2][0]], [1.15], 'md', ms=7, alpha=.8, zorder=11)
	ax.plot([foreshocks[2][0]], [1.15], 'kd', ms=9, alpha=.7, zorder=10)
	#
	#miny=min(map(operator.itemgetter(5), newratios))
	#maxy=max(map(operator.itemgetter(5), newratios))
	miny=.35
	maxy=6.
	dtstart=dtm.datetime(2004,6,1, tzinfo=pytz.timezone('UTC'))
	dtend=dtm.datetime(2005,6,1, tzinfo=pytz.timezone('UTC'))
	ax.set_ybound(lower=miny, upper=maxy)
	ax.set_xbound(lower=dtstart, upper=dtend)
	#
	tickDates  = [dtstart]
	tickLabels = ['%d/%d' % (dtstart.month, dtstart.year)]
	#
	delta_months=4
	while tickDates[-1]<=dtend:
		#
		newmonth_number = tickDates[-1].month + delta_months
		deltaYear = int(newmonth_number/12)
		#
		newmonth=1+newmonth_number%12
		#if newmonth==0: newmonth=12
		#
		#print "new date: ", dt1.year+dYear,newmonth,1
		#tickDates  += [dtm.datetime(dtstart.year+dYear,newmonth,1, tzinfo=pytz.timezone('UTC'))]
		tickDates += [dtm.datetime(tickDates[-1].year+deltaYear, newmonth, 1, tzinfo=pytz.timezone('UTC'))]
		
		#print tickDates[-1]
		tickLabels += ['%d/%d' % (tickDates[-1].month, tickDates[-1].year)]
		
	#myax.set_xticks(map(mpd.date2num, tickDates))
	ax.set_xticks(tickDates)
	ax.set_xticklabels(tickLabels)	
	#
	####################################################
	'''
	# i don't think this zoom really adds much to the figure.
	# next aftershock:
	shockIndex=3
	thisax = plt.figure(shockIndex).get_axes()
	thisax[2].set_xlim(dtm.datetime(1992, 11, 1, tzinfo=pytz.timezone('UTC')), dtm.datetime.now()+dtm.timedelta(days=30))
	thisax[2].set_ylim(.25, 8.6)
	#ax=thisax[2]
	#ax.arrow(x=mpd.date2num(dtm.datetime(2004,1,1, tzinfo=pytz.timezone('UTC'))), y=3.5, dx=365., dy=-2.25, length_includes_head=True, width=.003)
	axescolor=None
	thisax+=[plt.axes([.2, .55, .27, .2], frameon=True, axisbg=axescolor)]
	ax = thisax[4]
	winlen = mhp.getNsample(m=foreshocks[shockIndex][3], mc=mc)
	rbavelen=max(1, int(winlen/10.))
	r2=c1.plotIntervalRatiosAx(winlen=winlen, cat=c1.getcat(shockIndex), hitThreshold=1.0, bigmag=9.0, thisAx=ax, ratios=None, delta_t=1, avlen=rbavelen, mainEv=foreshocks[shockIndex], logZ=None, rbLegLoc=0, reverse=False)
	ax.plot([foreshocks[shockIndex][0]], [1.15], 'md', ms=7, alpha=.8, zorder=11)
	ax.plot([foreshocks[shockIndex][0]], [1.15], 'kd', ms=9, alpha=.7, zorder=10)
	#
	#miny=min(map(operator.itemgetter(5), newratios))
	#maxy=max(map(operator.itemgetter(5), newratios))
	miny=.35
	maxy=6.
	dtstart=dtm.datetime(2004,9,1, tzinfo=pytz.timezone('UTC'))
	dtend=dtm.datetime(2008,6,1, tzinfo=pytz.timezone('UTC'))
	ax.set_ybound(lower=miny, upper=maxy)
	ax.set_xbound(lower=dtstart, upper=dtend)
	#
	tickDates  = [dtstart]
	tickLabels = ['%d/%d' % (dtstart.month, dtstart.year)]
	#
	delta_months=4
	while tickDates[-1]<=dtend:
		#
		newmonth_number = tickDates[-1].month + delta_months
		deltaYear = int(newmonth_number/12)
		#
		newmonth=1+newmonth_number%12
		#if newmonth==0: newmonth=12
		#
		#print "new date: ", dt1.year+dYear,newmonth,1
		#tickDates  += [dtm.datetime(dtstart.year+dYear,newmonth,1, tzinfo=pytz.timezone('UTC'))]
		tickDates += [dtm.datetime(tickDates[-1].year+deltaYear, newmonth, 1, tzinfo=pytz.timezone('UTC'))]
		
		print tickDates[-1]
		tickLabels += ['%d/%d' % (tickDates[-1].month, tickDates[-1].year)]
		#
		#if len(tickDates)>10: break
	#myax.set_xticks(map(mpd.date2num, tickDates))
	ax.set_xticks(tickDates)
	ax.set_xticklabels(tickLabels)
	'''
	#
	#return arb
	return c1
	
def plotEvents(c1=None, m=None, mc=None, ax=None, startdt=None):
	#
	# this function specifically plots two sets of earthquake events -- earthquakes before and earthquakes after the
	# mainshock event (startdt).
	#
	if startdt==None: startdt=c1.getcat(1)[0][0]
	wlen = int(getNofm(m=m, mc=mc))
	X=[]
	Y=[]
	'''
	for rw in c1.getcat(0):
		if rw[0]>startdt:
			#shocks+=[rw]
			X+=[rw[2]]
			Y+=[rw[1]]
			#
		#
		if len(X)>wlen: break
	'''
	#
	# find start date event:
	ct=c1.getcat(0)
	for i in xrange(len(ct)):
		if ct[i][0]>startdt:
			mevIndex=i-1	# recognizing that this might be one event off...
			break
	#
	Xpre  = map(operator.itemgetter(2), ct[i-wlen:i])
	Ypre  = map(operator.itemgetter(1), ct[i-wlen:i])
	Xpost = map(operator.itemgetter(2), ct[i:i+wlen])
	Ypost = map(operator.itemgetter(1), ct[i:i+wlen])
	#
	# where are the unweighted CM of these groups?
	print "pre CM  (x,y): %f, %f" % (scipy.mean(Xpre), scipy.mean(Ypre))
	print "post CM (x,y): %f, %f" % (scipy.mean(Xpost), scipy.mean(Ypost))
	
	print "shock len: %d" % len(X)
	Xpre,Ypre = c1.catmap(Xpre,Ypre)
	Xpost,Ypost = c1.catmap(Xpost,Ypost)
	c1.catmap.plot(Xpre,Ypre, 'r.', ms=3, alpha=.7, zorder=11)
	c1.catmap.plot(Xpost,Ypost, 'b.', ms=3, alpha=.7, zorder=11)
	#
	


def tohokuQuads(targmag=9.0, mc=5.0, weighted_ratios=False):
	# mc=4.75 might be more desirable... we might even get away with 4.25 though... i think the mc is not consistent in the regions, so smaller mc gives lots of noise (moreso on th maps than the rbts).
	# ... and anyway, it looks like the m<4.5 events are totally FUBAR, at least spatially. they seem to produce decent RBTS, but really 
	# crappy hazmaps. there may be an error that mis-locates the smaller earthquakes, so we end up with the right number and maybe even 
	# temporal sequencing of small earthquakes, but they're in the wrong place. see HazMaps for mc=4.5 vs mc=5.0
	#
	# the Tohoku haz map remains disturbingly elusive since the ANSS magnitude corrections. we'll eventually need to re-tune, preferably
	# in a systematic way. it seems that mc=5.0 removes most of the noise.
	#
	# both revieweres asked about spatial constraints on RBTS. R2 specifically asked about including northern Japan. so here are some examples.
	#
	#a=mhp.japanhm(mc=4.5, todt=mhp.dtm.datetime.now(mhp.pytz.timezone('UTC')), ndithers=10, nContours=nContours1, bigmag=7.0, lons=[135., 146.25], lats=[30., 41.75], refreshcat=True, dt0=mhp.dtm.datetime(1990,1,1, tzinfo=mhp.pytz.timezone('UTC')))
	Nsequence = getNofm(targmag, mc)
	a=mhp.japanhm(mc=mc, todt=mhp.dtm.datetime.now(mhp.pytz.timezone('UTC')), winlen=110, ndithers=10, nContours=nContours1, bigmag=7.0, lons=[135., 148.5], lats=[30., 45.25], refreshcat=True, dt0=mhp.dtm.datetime(1990,1,1, tzinfo=mhp.pytz.timezone('UTC')))
	c1=a.getcat()
	print "nSubcats: ", len(c1.subcats)
	while len(c1.subcats): c1.subcats.pop()
	print "nSubcats: ", len(c1.subcats)
	#
	'''
	c1.addLatLonSubcat(fullcat=c1.getcat(0), lats=[31.0, 45.0], lons=[135.0, 146.0])
	c1.addLatLonSubcat(fullcat=c1.getcat(0), lats=[35.0, 45.0], lons=[135.0, 146.0])
	c1.addLatLonSubcat(fullcat=c1.getcat(0), lats=[34.0, 41.0], lons=[135.0, 146.0])
	c1.addLatLonSubcat(fullcat=c1.getcat(0), lats=[35.0, 40.5], lons=[135.0, 146.0])
	c1.addLatLonSubcat(fullcat=c1.getcat(0), lats=[35.0, 41.5], lons=[135.0, 146.0])
	c1.addLatLonSubcat(fullcat=c1.getcat(0), lats=[35.0, 40.5], lons=[140.0, 146.0])
	'''
	#
	# now, let's add some round catalogs. for fun, let's use Tahir nomenclature "normalized" distance.
	# from some experience, we've found that L_r ~ L(D=2).
	tohokuEvent=c1.getMainEvent()
	tohokuLatLon = tohokuEvent[1], tohokuEvent[2]
	Lr = 10.0**(.5*targmag - 1.76)
	print "Lr = %f" % Lr
	#rfactors = [.25, .25, .5, .6, .75, 1.0, 1.5, 2.0]
	rfactors = [.5]
	drawCircles=[]
	for thisr in rfactors:
		c1.subcats+=[['circ%f' % thisr, circularcat(incat=c1.getcat(), latlon=tohokuLatLon, Rkm=thisr*Lr)]]
		drawCircles+=[[len(c1.subcats), thisr]]
	#c1.subcats+=[['c05', circularcat(incat=c1.getcat(), latlon=tohokuLatLon, Rkm=.5*Lr)]]
	#c1.subcats+=[['c10', circularcat(incat=c1.getcat(), latlon=tohokuLatLon, Rkm=1.0*Lr)]]
	#c1.subcats+=[['c15', circularcat(incat=c1.getcat(), latlon=tohokuLatLon, Rkm=1.5*Lr)]]
	#c1.subcats+=[['c20', circularcat(incat=c1.getcat(), latlon=tohokuLatLon, Rkm=2.0*Lr)]]
	#c1.subcats+=[['c25', circularcat(incat=c1.getcat(), latlon=tohokuLatLon, Rkm=2.5*Lr)]]
	#
	print "nSubcats: ", len(c1.subcats)
	for k in xrange(len(c1.subcats)):
		c1.subcats[k][1].sort(key=lambda x: x[0])
	#
	fnum0=10
	print "make %d subcat plots." % (len(c1.subcats)+1)
	tohokuEvent = c1.getMainEvent(thiscat=c1.getcat(0))
	# get a circle of R=L_r/2
	circLR = getLLcircle(tohokuEvent[2], tohokuEvent[1], r=.5*Lr*1000., rtype=0, N=500)
	XcircLr, YcircLr=c1.catmap(scipy.array(circLR[0]), scipy.array(circLR[1]))
	intlist=map(int, (10.0**(4.5-mc))*scipy.array([250, 500, 770, 1000, 1250]))
	for i in xrange(len(c1.subcats)+1):
		thistargmag=targmag
		if i==2: thistargmag=8.4
		c1.rbomoriQuadPlot(mc=mc, targmag=thistargmag, bigmag=9.5, catnum=i, fignum=fnum0+i, plotevents=False, intlist=intlist, thislw=2.25, deltaLat=1.5, deltaLon=1.5, weighted=weighted_ratios)
		#
		# draw the catalog boundaries:
		llrange = c1.getLatLonRange(c1.getcat(i))
		#dlat=llrange[1][0]-llrange[0][0]
		#dlon=llrange[1][1]-llrange[0][1]
		xleft = llrange[0][1]
		xright = llrange[1][1]
		ydown = llrange[0][0]
		yup = llrange[1][0]
		#
		f=plt.figure(fnum0+i)
		thisAxes = f.get_axes()
		# plot Mainshock:
		ax=thisAxes[2]
		ax.plot([tohokuEvent[0]], [1.0], 'r^', ms=12, alpha=11)
		# draw catalog boundaries:
		ax=thisAxes[3]
		#
		if not(i>0 and c1.subcats[i-1][0][0:4]=='circ'):
			# it's not one of our circular catalogs...
			xsquare=[xleft, xright, xright, xleft, xleft]
			ysquare=[ydown, ydown, yup, yup, ydown]
			for j in xrange(len(xsquare)):
				x,y = c1.catmap(xsquare[j], ysquare[j])
				xsquare[j]=x
				ysquare[j]=y		
		#
			ax.plot(xsquare, ysquare, 'r-', lw=2.5, zorder=11)
		#
		# and plot the catalog. blue inside, green outside the bounds.
		Xin, Xout, Yin, Yout = [], [], [], []
		'''
		for ev in c1.getcat(0):
			#evcolor='g'
			x,y = c1.catmap(ev[2], ev[1])
			if ev[2]>xleft and ev[2]<xright and ev[1]<yup and ev[1]>ydown:
				# evcolor='b'
				Xin+=[x]
				Yin+=[y]
			else:
				Xout+=[x]
				Yout+=[y]
			#
			#x,y = c1.catmap(ev[2], ev[1])
		'''
		for ev in c1.getcat(0):
			x,y = c1.catmap(ev[2], ev[1])
			Xout+=[x]
			Yout+=[y]
		for ev in c1.getcat(i):
			x,y = c1.catmap(ev[2], ev[1])
			Xin+=[x]
			Yin+=[y]
			
		# plotting big lists instead of a bunch of single points seems to be less resource intensive.
		ax.plot(Xin, Yin, 'b.', ms=5, zorder=9, alpha=.7)
		ax.plot(Xout, Yout, 'g.', ms=5, zorder=8, alpha=.7)
		#
		# what about the preceeding events (the end of the precursory sequence -- it's hard to say how long it was)
		# what might be cool is to color code them -t "hot" from 0 (hottest -> -inf), +t "cool" with t->0 the "coolest"... later.
		#mevindex=c1.getMainEvent(thiscat=c1.getcat(i))[-1]
		mevindex = tohokuEvent[-1]
		#Nsequence = getNofm(targmag, mc)
		
		Nsequence2=int(Nsequence*1.5)
		Xprecursor = scipy.array(map(operator.itemgetter(2), c1.getcat(i)[mevindex-Nsequence2:mevindex]))
		Yprecursor = scipy.array(map(operator.itemgetter(1), c1.getcat(i)[mevindex-Nsequence2:mevindex]))
		Xp, Yp = c1.catmap(Xprecursor, Yprecursor)
		ax.plot(Xp, Yp, 'r.', zorder=11, ms=5, alpha=.4)
		#
		# and plot tohoku:
		xt, yt = c1.catmap(tohokuEvent[2], tohokuEvent[1])
		ax.plot([xt], [yt], 'r*', zorder=15, ms=18, alpha=.7)
		ax.plot([xt], [yt], 'k*', zorder=14, ms=21, alpha=.8)
		#
		ax.plot(XcircLr, YcircLr, 'k--', lw=2, alpha=.8, zorder=7)
		if i>0 and c1.subcats[i-1][0][0:4]=='circ':
			myrfact = float(c1.subcats[i-1][0][4:])
			print "plot circle R: ", myrfact, myrfact*Lr, " fnum=", fnum0+i
			circXY = getLLcircle(tohokuEvent[2], tohokuEvent[1], r=myrfact*Lr*1000., rtype=0, N=500)
			Xcirc, Ycirc=c1.catmap(scipy.array(circXY[0]), scipy.array(circXY[1]))
			ax.plot(Xcirc, Ycirc, 'r--', zorder=8, lw=2.5, alpha=.8)
	#
	'''
	for circ in drawCircles:
		fnum = circ[0]+10
		r=circ[1]
		# getcircle(lon0, lat0, r, cm, n=1000)
		circXY = getcircle(lon0=tohokuEvent[2], lat0=tohokuEvent[1], r=r, n=500, cm=c1.catmap)
		plt.figure(fnum)
		c1.catmap.plot(circXY[0], circXY[1], '--', zorder=8)
	'''
	#
	return c1

def getNofm(m, mc, mt=7.6, dmstar=1.0, dms=1.0):
	if mt==None: mt=7.6
	#
	targmag=m
	if targmag<mt:
		# "small" earthquake
		winlen=10**(targmag-dmstar-dms-mc)	# where 2.0 is dmstar + dmprime
	if targmag>=mt:
		dmsprime = (1.0*dms*(mt-mc) + 1.5*dms*(targmag-mt))/(targmag-mc)	# but really, we don't know how this works outside the SS limit
		#winlen = 10**(1.0*(mt-mc) + 1.5*(targmag-mt-1.0) - dms)
		winlen = 10**(1.0*(mt-mc) + 1.5*(targmag-mt) - 2.0*dmsprime)
	#
	#winlen=int(10*round(winlen/10))
	#print "winlen0: %d" % winlen
	winlen=int(round(winlen,-1))
	#print "winlen0: %d" % winlen
	if winlen<1: winlen=1
	#
	return winlen

def circularcat(incat, latlon, Rkm=10.):
	outcat=[]
	#g1=ggp.WGS84.Inverse(self.loc[1], self.loc[0], inloc[1], inloc[0])y
	#                    (y, x, y, x)
	#r=g1['s12']/1000.0
	#
	for rw in incat:
		g1 = ggp.WGS84.Inverse(rw[1], rw[2], latlon[0], latlon[1])
		r=g1['s12']/1000.0
		#print "rws: ", rw[1], rw[2], latlon[0], latlon[1], r, Rkm
		if r<=Rkm: outcat+=[rw]
	#
	return outcat
	
def ellipcat(incat, latlon, a=10.0, b=5.0, ellipTheta=0.0):
	outcat=[]
	ellipTheta = ellipTheta*math.pi/180.
	#a*=1000.0
	#b*=1000.0
	#g1=ggp.WGS84.Inverse(self.loc[1], self.loc[0], inloc[1], inloc[0])y
	#                    (y, x, y, x)
	#r=g1['s12']/1000.0
	centerlon=latlon[0]*math.pi/180.
	#
	for rw in incat:
		#
		#if rw==latlon: continue
		# distance from center to earthquake.
		g1 = ggp.WGS84.Inverse(rw[1], rw[2], latlon[0], latlon[1])
		r=g1['s12']/1000.0	# ... in km
		#print "rws: ", rw[1], rw[2], latlon[0], latlon[1], r, Rkm
		#
		# rough, recta-lin approximation for angle to this event.
		y=(rw[1]-latlon[0])	# dy in km. 
		x=(rw[2]-latlon[1])*math.cos(centerlon)	# for our purposes here, we don't need the 111.1 
		#
		if x==0: 
			if y>0.: theta =	90.0
			if y<0.: theta = -90.0
		else:
			theta = math.atan(y/x)
			
		#
		thetaPrime = theta-ellipTheta
		#
		# r in this direction:
		rprime = a*b/(math.sqrt((b*math.cos(thetaPrime))**2.0 + (a*math.sin(thetaPrime))**2.0  ))
		#
		if r<=rprime: outcat+=[rw]
	#
	return outcat

def getLLcircle(lon0, lat0, r, rtype=0, N=500):
	if rtype=='km': rtype=0
	if rtype=='deg': rtype=1
	#
	M1=ggp.ALL
	X,Y=[],[]
	rvals=[]
	theta=0.0
	dtheta=360./float(N)
	#
	while theta<=360.:
		g1=ggp.WGS84.GenDirect(lat0, lon0, theta, rtype, r, M1)
		# appears to return: (length(deg), lat, lon
		#rvals+=[g1[0:4]]
		X+=[g1[2]]
		Y+=[g1[1]]
		theta+=dtheta
		#
	#
	return [X,Y]

def getcircle(x0, y0, r, N=1000):
	# let's do this with proper circular symmetry.
	# and let's just be screwy about it: lon0,lat0 are lon, lat. r is in km.
	# cm is a basemap "catalog map" object.
	dtheta=2.0*math.pi/(float(N))
	#
	X=[]
	Y=[]
	theta=0.0
	while theta<=math.pi*2.0:
		dx = math.cos(theta)
		dy = math.sin(theta)
		X+=[x0+dx]
		Y+=[y0+dy]
		theta+=dtheta
	#
	# return as coords or map points? for now, just map points so we'll plot directly.
	#
	return [X,Y]

def getEllipser0(x0, y0, r0, epsilon=1.0, thetaEllipse=0.0, N=1000, a=None, b=None, rtype=0):
	# proper circular symmetry....
	# this is an "ab-ratio ellipse". we define the ratio of the major to minor axis
	# and "r0", the ratio of the circle with equivalent area. so, this is a transform of a 
	# circle with r0 to an ellipse with abrati=a/b such that both shapes have equal area.
	#
	# a, b, r0 in km...
	#
	# note: default gives us a circle...
	#
	# epsilon: given the r0 transformation: 
	# pi*r0^2 = pi*ab
	# ab=r0^2
	# a equiv= epsilon*r0
	# b=r0/epsilon		(by subbing into area)
	#
	if epsilon==None: epsilon=a/float(b)
	abratio=float(epsilon)
	#dtheta=2.0*math.pi/(float(N))
	dtheta = 360.0/float(N)
	#thetaEllipse = thetaEllipse*2.0*math.pi/360.
	#
	X=[]
	Y=[]
	abratio=math.sqrt(abratio)
	if a==None: a=epsilon*r0
	if b==None: b=r0/epsilon
	a*=1000.
	b*=1000.	# convert to meters.
	theta=0.0
	rvals=[]
	M1=ggp.ALL
	while theta<=(360. + dtheta):
		thetaPrime = (theta-thetaEllipse)*math.pi/180.
		#
		#xprime = R[0]*math.cos(thetaconv) - R[1]*math.sin(thetaconv)
		#yprime = R[0]*math.sin(thetaconv) + R[1]*math.cos(thetaconv)
		#
		# rpirme = ab/( (bcos(theta))^2 + (a*sin(theta))^2)
		rprime = a*b/(math.sqrt((b*math.cos(thetaPrime))**2.0 + (a*math.sin(thetaPrime))**2.0  ))
		g1=ggp.WGS84.GenDirect(y0, x0, theta, rtype, rprime, M1)
		# appears to return: (length(deg), lat, lon
		rvals+=[g1[0:4]]
		#		
		#dx=rprime*math.cos(theta)
		#dy=rprime*math.sin(theta)
		
		#
		#X+=[x0+dx]
		#Y+=[y0+dy]
		#
		theta+=dtheta
	#
	Rs=zip(*rvals)
	#print rvals[0:5]
	# return as coords or map points? for now, just map points so we'll plot directly.
	#
	#return [X,Y]
	return [Rs[2], Rs[1]]

def ellipsetest(r0=100000., abratio=2.0, thetaEllipse=0.0, rtype=0, N=500, lat0=37.7833, lon0=-122.4167):
	# rtype==1 (True): degrees
	# rtype==0 (False): meters
	#
	# create a "mask" object:
	M1=ggp.ALL
	X,Y=[],[]
	rvals=[]
	theta=0.0
	dtheta=360./float(N)
	#
	while theta<=360.:
		g1=ggp.WGS84.GenDirect(lat0, lon0, (theta-thetaEllipse), rtype, r0, M1)
		# appears to return: (length(deg), lat, lon
		rvals+=[g1[0:4]]
		theta+=dtheta
	#
	Rs=zip(*rvals)
	plt.ion()
	plt.figure(0)
	plt.clf()
	plt.plot(Rs[2], Rs[1], '-')
	
	#plt.figure(1)
	#plt.clf()
	#plt.plot(Rs[3], Rs[0], '-')
	#plt.plot(Rs[3], Rs[1], '-')
	#plt.plot(Rs[3], Rs[2], '-')
	#

def circletest(r=100000., rtype=0, N=500, lat0=37.7833, lon0=-122.4167):
	# rtype==1 (True): degrees
	# rtype==0 (False): meters
	#
	# create a "mask" object:
	M1=ggp.ALL
	X,Y=[],[]
	rvals=[]
	theta=0.0
	dtheta=360./float(N)
	#
	while theta<=360.:
		g1=ggp.WGS84.GenDirect(lat0, lon0, theta, rtype, r, M1)
		# appears to return: (length(deg), lat, lon
		rvals+=[g1[0:4]]
		theta+=dtheta
	#
	Rs=zip(*rvals)
	plt.ion()
	plt.figure(0)
	plt.clf()
	plt.plot(Rs[2], Rs[1], '-')
	
	plt.figure(1)
	plt.clf()
	plt.plot(Rs[3], Rs[0], '-')
	plt.plot(Rs[3], Rs[1], '-')
	plt.plot(Rs[3], Rs[2], '-')
	#


def tohokuHMsigmas(targmag=9.0, mc=4.5,fcdt = mhp.dtm.datetime(2011,3,5, tzinfo=mhp.pytz.timezone('UTC')), ndithers=10):
	# mc=4.75 might be more desirable... we might even get away with 4.25 though... i think the mc is not consistent in the regions, so smaller mc gives lots of noise (moreso on th maps than the rbts).
	#
	# the Tohoku haz map remains disturbingly elusive since the ANSS magnitude corrections. we'll eventually need to re-tune, preferably
	# in a systematic way.
	#
	# both revieweres asked about spatial constraints on RBTS. R2 specifically asked about including northern Japan. so here are some examples.
	#
	#a=mhp.japanhm(mc=4.5, todt=mhp.dtm.datetime.now(mhp.pytz.timezone('UTC')), ndithers=10, nContours=nContours1, bigmag=7.0, lons=[135., 146.25], lats=[30., 41.75], refreshcat=True, dt0=mhp.dtm.datetime(1990,1,1, tzinfo=mhp.pytz.timezone('UTC')))
	Nsequence = getNofm(targmag, mc)
	
	#fcdt = mhp.dtm.datetime(2011,3,5, tzinfo=mhp.pytz.timezone('UTC'))
	a=mhp.japanhm(mc=mc, todt=fcdt, ndithers=ndithers, nContours=nContours, bigmag=8.0, lons=[135., 148.5], lats=[30., 45.25], refreshcat=True, dt0=mhp.dtm.datetime(1990,1,1, tzinfo=mhp.pytz.timezone('UTC')))
	c1=a.getcat()
	#tohokuY, tohokuX=c1.getMainEvent()[1:3]
	tohokuX, tohokuY=c1.catmap(142.369, 38.322)
	print "Tohoku: " , tohokuX, tohokuY
	mysigma=1.55
	dsigma=.01
	fnum=0
	dtstr = 'dt%s%s%s' % (str(fcdt.year), '00'+str(fcdt.month)[-2:], '00'+str(fcdt.day)[-2:])
	while mysigma<1.75:
		a=mhp.japanhm(mc=mc, todt=fcdt, ndithers=ndithers, nContours=nContours, bigmag=8.0, lons=[135., 148.5], lats=[30., 45.25], refreshcat=False, dt0=mhp.dtm.datetime(1990,1,1, tzinfo=mhp.pytz.timezone('UTC')), sigma=mysigma)
		mysigma+=dsigma
		#fnum+=3
		plt.figure(3)
		plt.plot([tohokuX], [tohokuY], 'r*', ms=18, alpha=.7, zorder=10)
		plt.plot([tohokuX], [tohokuY], 'k*', ms=21, alpha=.7, zorder=9)
		savedir = '/home/myoder/Dropbox/Research/rbPrecursorsTecto/tohokuHMs/%s/dith%d' % (dtstr, ndithers)
		dirses=savedir.split('/')
		mydir=''
		for thisdir in dirses:
			mydir = mydir + '/' + thisdir
			if glob.glob(mydir)==[]: os.system('mkdir %s' % mydir)
		#if glob.glob(savedir)==[]: os.system('mkdir %s' % savedir)
		print 'saving: %s/tohokuhm-sigma-%f-dith%d.png' % (savedir, mysigma, ndithers)
		plt.savefig('%s/tohokuhm-sigma-%f-dith%d.png' % (savedir, mysigma, ndithers))
	#
	return a.getcat()
	
def winlen(m, mc, mt=7.6, doInt=True):
	if m<mt:
		# "small" earthquake
		winlen=10**(m-2.0-mc)	# where 2.0 is dmstar + dmprime
	if m>=mt:
		dms = (1.0*(mt-mc) + 1.5*(m-mt))/(m-mc)
		#winlen = 10**(1.0*(mt-mc) + 1.5*(targmag-mt-1.0) - dms)
		winlen = 10**(1.0*(mt-mc) + 1.5*(m-mt) - 2.0*dms)
	#
	if doInt:
		winlen=int(round(winlen,-1))
		#print "winlen0: %d" % winlen
		if winlen<1: winlen=1
	#
	return winlen
		
