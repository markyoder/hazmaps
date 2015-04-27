import math
import scipy
#
import pylab as plt
import matplotlib as mpl
import numpy
import scipy.optimize as spo
from matplotlib.patches import Ellipse
import matplotlib.dates as mpd
#
import string
import sys
import os
import random
import time
import pytz
import glob
#
# gamma function lives here:
#import scipy.special
from scipy.special import gamma
from matplotlib import axis as aa
import threading
#from threading import Thread
#
import datetime as dtm
import calendar
import operator
import urllib
#import MySQLdb

#import yodapy as yp 
#import rbIntervals as rbi
import hazmap as hmp
import BASScast as bcp
import eqcatalog as eqp
import ANSStools as atp

#
# some notes:
# basemap(lon, lat) returns (x,y) in metres; "inverse" keyword inverts the process (UTM->lat/lon).
# so, in the next version of these tools, we will specify bin sizes in km (or m, or otherwise linear distances).
# will will first have to convert the catalog from (lon, lat) -> (x, y)
# i think the 'laea' projection option (Lambert Azimuthal Equal Area) will do the trick:
# convert all lat/lon in catalog to laea coordinates.
# bin by km.
# calc. rb-ratios, etc.
# convert bin coordinates to lat/lon and plot. there may be some overlapping
#


#def __init__(self, catalog=None, catnum=0, mc=2.5, gridsize=.5, ndithers=3, winlen=128, avlen=1, bigmag=5.0, fignum=0)
def makeHazMapDith(catname='cats/japancat4b.cat', catnum=0, mc=4.5, gridsize=None, ndithers=3, winlen=20, sigma=1.68, avlen=None, bigmag=6.5, fignum=0, logZ=None):
	if gridsize==None and winlen==None:
		# can't do anything with that, can we?
		return None
	midlat=getMidLat(catname, mc)
	if gridsize==None:
		gridsize=getgridsize(mc, winlen, midlat, sigma)
	if winlen==None:
		winlen=getWinlen(mc, gridsize, midlat, sigma)
	if avlen==None:
		avlen = int(winlen/10)
		if avlen==0: avlen=1

	# generic ditherd haz-map maker. then, hazMapTo(your-favorite-date-of-catastrophe).
	catalog=eqp.eqcatalog()
	catalog.loadCatFromFile(catname)
	#
	return makeHazMapDith2(catalog, catnum, mc, gridsize, ndithers, winlen, sigma, avlen, bigmag, fignum, logZ)

def makeHazMapDith2(catalog=None, catnum=0, mc=4.5, gridsize=None, ndithers=3, winlen=16, sigma=1.68, avlen=None, bigmag=6.5, fignum=0, logZ=None):
	#
	if catalog==None: return None
	if type(catalog).__name__=='string': return makeHazMapDith(catalog, catnum, mc, gridsize, ndithers, winlen, sigma, avlen, bigmag, fignum, logZ=logZ)	
	if gridsize==None and winlen==None:
		# can't do anything with that, can we?
		return None
	#midlat=getMidLat(catname, mc)
	midlat = scipy.mean(map(operator.itemgetter(1), catalog.getcat(catnum)))
	
	if gridsize==None:
		gridsize=getgridsize(mc, winlen, midlat, sigma)
	if winlen==None:
		winlen=getWinlen(mc, gridsize, midlat, sigma)
	if avlen==None:
		avlen = int(winlen/10)
		if avlen==0: avlen=1

	#
	llrange=catalog.getLatLonRange()
	centerlat=llrange[0][0] + (llrange[1][0]-llrange[0][0])/2.
	print centerlat
	#mysigmaexp=1.68
	mysigmaexp = sigma
	#	
	if gridsize==None and winlen==None: winlen=20	
	if winlen==None:
		# we've been given gridsize; calculate the winlen.
		winlen=getWinlen(mc, gridsize, centerlat, mysigmaexp)
		print "set winlen: %d from gridsize: %f" % (winlen, gridsize)
	if gridsize==None:
		gridsize=getgridsize(mc=mc, winlen=winlen, lat=centerlat, sigmaExp=mysigmaexp)
		print "set gridsize: %f from winlen: %d" % (gridsize, winlen)
	
	#objHM=hmp.hazCrack(catalog, catnum, mc, gridsize, ndithers, winlen, avlen, bigmag, fignum)
	objHM=hmp.hazCrack(catalog=catalog, catnum=catnum, mc=mc, gridsize=gridsize, ndithers=ndithers, winlen=winlen, avlen=avlen, bigmag=bigmag, fignum=fignum, logZ=logZ)
	#
	return objHM

def chichiStandardHM(catname='cats/chichistandard.cat', catnum=0, mc=3.0, gridsize=None, ndithers=5, winlen=20, sigma=1.68, avlen=None, bigmag=6.3, fignum=0, dt0=None, todt=(dtm.datetime.now(pytz.timezone('UTC'))), epicen=None, thresh=None, refreshcat=True, lons=[118.0, 125.7], lats=[20.0, 26.2], nContours=None):
	if dt0==None: dt0=todt - dtm.timedelta(days=int(365.24*12))
	if nContours==None:
		nContours=[-.9, -.8, -.7, -.6, -.5, -.4, -.3, -.2, -.1, 0.0, .1, .2, .3, .4, .5, .6, .7]
	return generalHazmap(catname=catname, catnum=catnum, mc=mc, gridsize=gridsize, ndithers=ndithers, winlen=winlen, sigma=sigma, avlen=avlen, bigmag=bigmag, fignum=fignum, dt0=dt0, todt=todt, epicen=epicen, thresh=thresh, refreshcat=refreshcat, lats=lats, lons=lons, nContours=nContours)

def makeChiChimap(catname='cats/taiwan1994.cat', catnum=0, mc=2.5, gridsize=None, ndithers=3, winlen=32, sigma=1.68, avlen=1, bigmag=6.0, fignum=0):
	#
	
	taiwancat=eqp.eqcatalog()
	taiwancat.loadCatFromFile('cats/taiwan.cat')
	mytcat=eqp.eqcatalog(taiwancat.getTimeRangeCat(taiwancat.getcat(0), dtm.datetime(1985,1,1, 0, 0, 0, 0, pytz.timezone('UTC')), dtm.datetime(2002,1,1, 0, 0, 0, 0, pytz.timezone('UTC'))))
	#
	chichidate = dtm.datetime( 1999, 9, 20, 17, 47, 15, 850000, tzinfo=pytz.timezone('UTC'))
	chichi = [chichidate, 23.8525, 120.8155, 7.3, 8.0]
	chichiEpicen=[chichi[1], chichi[2]]
	#
	#objHM=makeHazMapDith(catname, catnum, mc, gridsize, ndithers, winlen, sigma, avlen, bigmag, fignum)
	objHM=makeHazMapDith2(mytcat, catnum, mc, gridsize, ndithers, winlen, sigma, avlen, bigmag, fignum)
	objHM.epicen=chichiEpicen
	#objHM.hazMapTo(dtm.datetime(1999, 9, 10, 0, 0, 0, 0, pytz.timezone('UTC')))
	objHM.hazMapTo(dtm.datetime(1999, 9, 10, 0, 0, 0, 0, pytz.timezone('UTC')))	# let's fix everything to UTC; be sure we stop
																											# in advance of the mainshocks.
	objHM.simpleBoxes()
	hmx, hmy=objHM.hazmaps[0].catalog.catmap(120.98, 23.77)
	plt.plot([hmx], [hmy], 'r*', ms=18)
	#
	
	z2=objHM.contourMap(fignum+2)
	plt.figure(fignum+2)
	plt.plot([hmx], [hmy], 'b*', ms=18)
	#
	z=objHM.contourPlot(fignum=1)
	plt.figure(1)
	plt.plot([120.98], [23.77], 'r*', ms=18)
	#
	
	return objHM

def makecalnevmap(catname='cats/calnev1980.cat', catnum=0, mc=2.5, gridsize=None, ndithers=3, winlen=64, sigma=1.68, avlen=1, bigmag=5.5, fignum=0, todt=dtm.datetime.now(pytz.timezone('UTC'))):
	#
	objHM=makeHazMapDith(catname, catnum, mc, gridsize, ndithers, winlen, sigma, avlen, bigmag, fignum)
	#
	# let's use 1980 forward. incomplete data and low mc introduce lots of noise that looks like r<1 (data are added at an artificially increasing rate
	# if mc is systematicaly decreasing. obviously, this is a messy way to do this:
	#dtm0=dtm.datetime(1980,1,1, 0, 0, 0, 0, pytz.timezone('UTC'))
	#for hm in objHM.hazmaps:
	#	while hm.catalog.getcat(0)[0][0]<dtm0: hm.catalog.getcat(0).pop(0)
	#
	#objHM.hazMapTo(dtm.datetime(1999, 9, 10, 0, 0, 0, 0, pytz.timezone('UTC')))
	objHM.hazMapTo(todt)
	#objHM.simpleBoxes()
	objHM.simpleBoxes(fignum=fignum)
	objHM.contourMap(fignum=fignum+2)
	for fn in [fignum, fignum+2, fignum+3]:
		objHM.hazmaps[0].catalog.catmap.drawstates()

	#
	hmx, hmy=objHM.hazmaps[0].catalog.catmap(-121.1	, 35.7)	# sansim
	#plt.plot([hmx], [hmy], 'b*', ms=18)
	for fn in [fignum, fignum+2, fignum+3]:
		plt.figure(fn)
		plt.plot([hmx], [hmy], 'b*', ms=18)
	
	return objHM

def makeparkfieldmap(catname='cats/parkfield.cat', catnum=0, mc=1.5, gridsize=None, ndithers=3, winlen=60, sigma=1.68, avlen=6, bigmag=5.5, fignum=0, todt=dtm.datetime(2004,9,20, 0, 0, 0, 0, pytz.timezone('UTC'))):
	# for reference, we can make the aftershock plots like this:
	#
	# get gridsize/winlen:
	# getgridsize(mc=4.5, winlen=16, lat=45., sigmaExp=1.68)
	# def getWinlen(mc=4.5, gridsize=1.0, lat=45., sigmaExp=1.68):
	if gridsize==None and winlen==None:
		# can't do anything with that, can we?
		return None
	midlat=getMidLat(catname, mc)
	if gridsize==None:
		gridsize=getgridsize(mc, winlen, midlat, sigma)
	if winlen==None:
		winlen=getWinlen(mc, gridsize, midlat, sigma)
	#
	# noise levels:
	noisethresh=((winlen/2.)-.577 - math.log(winlen))/((winlen/2.)+.577 + math.log(winlen))
	noisethresh2=2.*(.577+math.log(winlen))/winlen
	noisethresh=abs(math.log10(noisethresh))
	noisethresh2=abs(math.log10(noisethresh2))
	print "noise thresholds: %f, %f" % (noisethresh, noisethresh2)
	#
	# pfc=eqp.eqcatalog()
	# pfc.loadCatFromFile('cats/parkfield.cat')
	#
	# pfcEl=eqp.eqcatalog(pfc.ellipseCat(None, -40., 35.9, -120.5, .15, .4))
	#
	# plot an ellipse around the ellipse catalog:
	# el=Ellipse([-120.5, 35.9], .8, .3, -40., facecolor='b', alpha=.25)
	# Xel, Yel = pfc.catmap(el.get_verts()[:,0],el.get_verts()[:,1])
	# mhp.plt.fill(Xel, Yel, ec='r', fc='c', alpha=.25)
	#
	# but for the haz-map, we don't need to worry about that (though we might want to hilight the aftershock zone in a plot)
	pfc0=eqp.eqcatalog()
	#pfc0.loadCatFromFile('cats/parkfield.cat')
	pfc0.loadCatFromFile(catname)
	pfc=eqp.eqcatalog(pfc0.getLatLonSubcat(pfc0.getcat(0), lats=[35.65, 36.3], lons=[-121.0, -120.]))
	#
	#objHM=makeHazMapDith(catname, catnum, mc, gridsize, ndithers, winlen, sigma, avlen, bigmag, fignum)
	objHM=makeHazMapDith2(pfc, catnum, mc, gridsize, ndithers, winlen, sigma, avlen, bigmag, fignum)
	#objHM.hazMapTo(dtm.datetime(2004, 9, 20, 0, 0, 0, 0, pytz.timezone('UTC')))
	objHM.epicen=None
	objHM.hazMapTo(todt)
	#objHM.hazMapTo(dtm.datetime.now(pytz.timezone('UTC')))
	objHM.simpleBoxes(mapres='l', thresh=noisethresh)
	
	for fn in [fignum, fignum+2]:
		hmx, hmy=objHM.hazmaps[0].catalog.catmap(-120.37	, 35.82)
		plt.figure(fn)
		plt.plot([hmx], [hmy], 'b*', ms=18)
	#
	'''
	plt.figure(fignum+1)
	plt.clf()
	z=objHM.contourPlot(fignum=fignum+1)
	plt.figure(fignum+1)
	plt.plot([-120.37], [35.82], 'r*', ms=18)
	'''
	
	z2=objHM.contourMap(fignum=fignum+2)
	plt.figure(fignum+2)
	plt.plot([hmx], [hmy], 'b*', ms=18)
	
	return objHM

def makesansimmap(catname='cats/sansim.cat', catnum=0, mc=1.5, gridsize=None, ndithers=3, winlen=64, sigma=1.68, avlen=1, bigmag=5.5, fignum=0):
	objHM=makeHazMapDith(catname, catnum, mc, gridsize, ndithers, winlen, sigma, avlen, bigmag, fignum)
	objHM.hazMapTo(dtm.datetime(2003, 12, 10, 0, 0, 0, 0, pytz.timezone('UTC')))
	#objHM.hazMapTo(dtm.datetime.now(pytz.timezone('UTC')))
	objHM.simpleBoxes()
	
	return objHM
	
def parseNZcat(finname=None, foutname=None):
	# parse catalog data from geonet.ort.nz
	fin=open(finname)
	fout=open(foutname, 'w')
	#
	for rw in fin:
		if rw[0]=='#' or rw[0:4] == 'CUSP': continue
		rws=rw.split(',')
		#print rws
		fout.write('%s/%s/%s %s:%s:%s\t%s\t%s\t%s\t%s\n' % (rws[5], rws[6], rws[7], rws[8], rws[9], rws[10], rws[1], rws[2], rws[11], rws[12] ))
	#
	fin.close()
	fout.close()
	#
	return None	
	
def makechristchurch(catname='cats/christchurchzoom.cat', catnum=0, mc=3.0, gridsize=None, ndithers=5, winlen=16, sigma=1.68, avlen=1, bigmag=5.5, fignum=0):
	# 'cats/christchurch-gn.cat'
	# this map really does not give much of a result. the bottom line is there is not enough data; it might make a good plot to discuss
	# minumum detection threshold. basically, m_min = m_c + dm* + dm_s + log(winlen) (aka, the estimated mag for 1 full bin). for mc=2.5, winlen=16, this
	# is m_min=6.2. at winlen=16, we seem to get a fair bit of noise, and for longer winlens, we cannot resolve the christchurch event (m=6.3).
	#
	# but, the prams {mc=3.0, winlen=16, ndithers=5} produces a usable result, albeit close to the noise limit.
	# this plot is challenging because it is nested in an aftershock sequence (this event is an aftershock of the canterbury event).
	# we like to think this method is independent of the aftershock situation, but in the end the wave-forms step all over one another
	# and produce noise and aliasing.
	# 
	objHM=makeHazMapDith(catname, catnum, mc, gridsize, ndithers, winlen, sigma, avlen, bigmag, fignum)
	objHM.hazMapTo(dtm.datetime(2011, 2, 15, 0, 0, 0, 0, pytz.timezone('UTC')))
	#objHM.hazMapTo(dtm.datetime.now(pytz.timezone('UTC')))
	objHM.simpleBoxes()
	#
	# christchurch (2011 2 22) 6.3
	hmx, hmy=objHM.hazmaps[0].catalog.catmap(172.7, -43.58)	
	plt.plot([hmx], [hmy], 'b*', ms=18, zorder=10, alpha=.4)
	#
	# canterbury (2010 9 4) 7.1
	# (but there is simply no data for this.)
	hmx, hmy=objHM.hazmaps[0].catalog.catmap(172.18, -43.55)	
	plt.plot([hmx], [hmy], 'g*', ms=18, zorder=10, alpha=.4)
	
	return objHM
	
def makeBigBearMap(catname='cats/hmine.cat', catnum=0, mc=2.5, gridsize=None, ndithers=3, winlen=32, sigma=1.68, avlen=1, bigmag=5.5, fignum=0):
	# 1992/06/28
	# 1992/06/28 11:57:34.13	34.2	-116.437	7.3
	objHM=makeHazMapDith(catname, catnum, mc, gridsize, ndithers, winlen, sigma, avlen, bigmag, fignum)
	objHM.hazMapTo(dtm.datetime(1992, 6, 26, 0, 0, 0, 0, pytz.timezone('UTC')))
	objHM.simpleBoxes()
	#
	hmx, hmy=objHM.hazmaps[0].catalog.catmap(-116.437	, 34.2)
	plt.plot([hmx], [hmy], 'b*', ms=18)
	#
	return objHM
	
def getHMinecat():
	# def catfromANSS(lon=[135., 150.], lat=[30., 41.5], minMag=4.0, dates0=[dtm.date(2005,01,01), None], Nmax=999999, fout='cats/mycat.cat'):
	atp.catfromANSS(lon=[-117.25, -115.5], lat=[33.5, 35.5], minMag=1.5, dates0=[dtm.datetime(1975,1,1, tzinfo=pytz.timezone('UTC')), None], fout='cats/hmine.cat')
	return None

def makeLandersMap(catname='cats/hmine.cat', catnum=0, mc=2.5, gridsize=None, ndithers=3, winlen=20, sigma=1.68, avlen=2, bigmag=6.0, fignum=0, dateto=dtm.datetime(1992, 6, 25, tzinfo=pytz.timezone('UTC'))):
	# 1999/10/16 09:46:44.13
	# 1999/10/16 09:46:44.13	34.594	-116.271	7.1
	objHM=makeHazMapDith(catname, catnum, mc, gridsize, ndithers, winlen, sigma, avlen, bigmag, fignum)
	# (1992, 6, 28)
	objHM.hazMapTo(dateto)
	objHM.simpleBoxes(fignum=fignum)
	objHM.contourMap(fignum=fignum+2)
	#hmx, hmy=objHM.hazmaps[0].catalog.catmap(-116.433, 34.217)
	hmx, hmy=objHM.hazmaps[0].catalog.catmap(-116.26, 34.13)

	#
	# plot HM:
	#hmx, hmy=objHM.hazmaps[0].catalog.catmap(-116.27, 34.59)
	for fn in [fignum, fignum+2, fignum+3]:
		plt.figure(fn)
		plt.plot([hmx], [hmy], 'b*', ms=18)
	
	return objHM

	
def makeHmineMap(catname='cats/hmine.cat', catnum=0, mc=2.5, gridsize=None, ndithers=3, winlen=32, sigma=1.68, avlen=1, bigmag=6.0, fignum=0):
	# 1999/10/16 09:46:44.13
	# 1999/10/16 09:46:44.13	34.594	-116.271	7.1
	objHM=makeHazMapDith(catname, catnum, mc, gridsize, ndithers, winlen, sigma, avlen, bigmag, fignum)
	objHM.epicen=None
	objHM.hazMapTo(dtm.datetime(1999, 10, 16, 0, 0, 0, 0, pytz.timezone('UTC')))
	objHM.simpleBoxes()
	objHM.contourMap(fignum=fignum+2)
	hmx, hmy=objHM.hazmaps[0].catalog.catmap(-116.27, 34.59)
	#
	# plot HM:
	for fn in [fignum, fignum+2]:
		plt.figure(fn)
		plt.plot([hmx], [hmy], 'b*', ms=18)
	#
	
	return objHM
	
def makeJapanHazMapDith(catname='cats/japancat.cat', catnum=0, mc=4.5, gridsize=None, ndithers=3, winlen=20, sigma=1.68, avlen=None, bigmag=7.0, fignum=0, dateto=dtm.datetime(2011, 3, 10, 0, 0, 0, 0, pytz.timezone('UTC')), epicen=None, thresh=0, logZ=1.0, dorefresh=False):
	if dorefresh==True:
		# makeJapCat(lon=[135., 146.], lat=[30., 41.5], minMag=4.5, dates0=[dtm.date(1990,01,01), None], Nmax=999999, fout='cats/japcat0.cat')
		thisnow=dtm.datetime.now(pytz.timezone('UTC'))
		a=makeJapCat(fout=catname, minMag=mc-.5, dates0=[dtm.datetime(2000,01,01, tzinfo=pytz.timezone('UTC')), thisnow])
	if avlen==None:
		avlen=int(winlen/10)
		if avlen<1: avlen=1
		print "avlen: %d" % avlen
	objHM=makeHazMapDith(catname, catnum, mc, gridsize, ndithers, winlen, sigma, avlen, bigmag, fignum, logZ=None)
	objHM.bigmag=bigmag
	objHM.epicen=epicen
	#objHM.hazMapTo(dtm.datetime(2011, 3, 11, 0, 0, 0, 0, pytz.timezone('UTC')))
	objHM.hazMapTo(dateto)
	#
	objHM.simpleBoxes(fignum=fignum, thresh=thresh, plotevents=0, mapres='l')
	objHM.hazmaps[0].catalog.catmap.drawmapscale(lon=137., lat=39., lon0=141., lat0=36.0, length=220., labelstyle='simple', barstyle='fancy')
	#
	# mainshock: 2011/3/11 5:45:24.56
	# 142,37, 38.30, m9.0
	evll=[142.37, 38.30]
	x,y=objHM.hazmaps[0].catalog.catmap(evll[0], evll[1])
	plt.figure(fignum)
	#plt.plot([x], [y], 'r*', ms=20)
	#
	objHM.contourMap(fignum=fignum+2, nconts=24)
	plt.figure(fignum+2)
	objHM.hazmaps[0].catalog.catmap.drawmapscale(lon=137., lat=39., lon0=141., lat0=36.0, length=220., labelstyle='simple', barstyle='fancy')
	plt.figure(fignum+3)
	objHM.hazmaps[0].catalog.catmap.drawmapscale(lon=137., lat=39., lon0=141., lat0=36.0, length=220., labelstyle='simple', barstyle='fancy')
	#plt.plot([x], [y], 'r*', ms=20)
	
	for fn in [fignum, fignum+2, fignum+3]:
		plt.figure(fn)
		plt.plot([x], [y], 'r*', ms=20, zorder=42)
	#
	plt.title('%d-%d-%d, $d\\lambda=%.4f, nDith$=$%d$, $m_c=%.2f$, $N_{rb} = %d$, $N_{ave}=%d$\n\n' % (dateto.year, dateto.month, dateto.day, objHM.gridsize, objHM.ndithers, objHM.mc, objHM.winlen, objHM.avlen))
	#
	return objHM

def makeLaquilaHazmap(catname='cats/laquila.cat', catnum=0, mc=2.5, gridsize=None, ndithers=3, winlen=10, sigma=1.68, avlen=1, bigmag=5.0, fignum=0, dateto=dtm.datetime(2009, 4, 1, 0, 0, 0, 0, pytz.timezone('UTC')), epicen=[13.38, 42.35]):
	#mainshock: 2009, 4, 6
	objHM=makeHazMapDith(catname, catnum, mc, gridsize, ndithers, winlen, sigma, avlen, bigmag, fignum)
	objHM.bigmag=bigmag
	objHM.epicen=epicen
	#objHM.hazMapTo(dtm.datetime(2011, 3, 11, 0, 0, 0, 0, pytz.timezone('UTC')))
	objHM.hazMapTo(dateto)
	objHM.simpleBoxes(fignum=fignum, thresh=-0.0, plotevents=0)
	#
	# mainshock: 2011/3/11 5:45:24.56
	# 142,37, 38.30, m9.0
	evll= epicen # [142.37, 38.30]
	x,y=objHM.hazmaps[0].catalog.catmap(evll[0], evll[1])
	plt.figure(fignum)
	#plt.plot([x], [y], 'r*', ms=20)
	#
	objHM.contourMap(fignum=fignum+2)
	plt.figure(fignum+2)
	#plt.plot([x]d, [y], 'r*', ms=20)
	
	for fn in [fignum, fignum+2, fignum+3]:
		plt.figure(fn)
		plt.plot([x], [y], 'r*', ms=20)
	#
	return objHM
#
def almanormovie(catname='cats/almanor.cat', catnum=0, mc=2.25, gridsize=None, ndithers=10, winlen=20, sigma=1.68, avlen=None, bigmag=3.0, fignum=0, dts=[dtm.datetime(2012, 5, 20, 0, 0, 0, 0, pytz.timezone('UTC')), dtm.datetime.now(pytz.timezone('UTC')), dtm.timedelta(hours=12)], catdate=dtm.datetime(1990, 1, 1, 0, 0, 0, 0, pytz.timezone('UTC')), epicen=None, thresh=None, refreshcat=True, lats=[39.5, 40.5], lons=[-121.7, -120.4], logZ=None, nameroot='almanor', moviedir='movies/almanor/', fps=15):
	todt=dts[1]
	dt0=catdate
	# 2012/7/1 3:25:20.450000	32.2072	-115.3432	4.64	3.18

	quakes = []
	quakes += [[dts[0], dts[1], -121.06, 40.19, 4.7]]
	#
	hm=generalHazMovie(catname=catname, catnum=catnum, mc=mc, gridsize=gridsize, ndithers=ndithers, winlen=winlen, sigma=sigma, avlen=avlen, bigmag=bigmag, fignum=fignum, dts=dts, catdate=catdate, epicen=epicen, thresh=thresh, refreshcat=refreshcat, lats=lats, lons=lons, logZ=logZ, nameroot=nameroot, moviedir=moviedir, events=quakes, fps=fps)
	#
	#os.system('avconv -i %s/contour/%s-%05d.png -r 50 %s/%s.mp4', (moviedir, nameroot, moviedir,nameroot))
	#
	return hm
#
def hminemovie(catname='cats/hmine.cat', catnum=0, mc=3.0, gridsize=None, ndithers=10, winlen=20, sigma=1.68, avlen=None, bigmag=7.9, fignum=0, dts=[dtm.datetime(1997, 1, 1, 0, 0, 0, 0, pytz.timezone('UTC')), dtm.datetime(2002, 1, 1, 0, 0, 0, 0, pytz.timezone('UTC')), dtm.timedelta(days=1)], catdate=dtm.datetime(1990, 1, 1, 0, 0, 0, 0, pytz.timezone('UTC')), epicen=None, thresh=None, refreshcat=True, lats=[33.5, 34.75], lons=[-117.25, -115.5], logZ=None, nameroot='hmine', moviedir='movies/hmine/'):
	todt=dts[1]
	dt0=catdate
	# 2012/7/1 3:25:20.450000	32.2072	-115.3432	4.64	3.18

	quakes = []
	quakes += [[dts[0], dts[1], -116.267, 34.6, 7.1]]
	#
	hm=generalHazMovie(catname=catname, catnum=catnum, mc=mc, gridsize=gridsize, ndithers=ndithers, winlen=winlen, sigma=sigma, avlen=avlen, bigmag=bigmag, fignum=fignum, dts=dts, catdate=catdate, epicen=epicen, thresh=thresh, refreshcat=refreshcat, lats=lats, lons=lons, logZ=logZ, nameroot=nameroot, moviedir=moviedir, events=quakes)
	#
	#os.system('avconv -i %s/contour/%s-%05d.png -r 50 %s/%s.mp4', (moviedir, nameroot, moviedir,nameroot))
	#
	return hm
#
# 2012/7/1 3:25:20.450000	32.2072	-115.3432	4.64	3.18
def  brawleyswarm(catname='cats/brawley.cat', catnum=0, mc=2.25, gridsize=None, ndithers=10, winlen=20, sigma=1.68, avlen=None, bigmag=7.9, fignum=0, dts=[dtm.datetime(2011, 9, 1, 0, 0, 0, 0, pytz.timezone('UTC')), dtm.datetime.now(pytz.timezone('UTC')), dtm.timedelta(days=1)], catdate=dtm.datetime(1990, 1, 1, 0, 0, 0, 0, pytz.timezone('UTC')), epicen=None, thresh=None, refreshcat=True, lats=[32.85, 33.1], lons=[-115.7, -115.3], logZ=None, nameroot='brawley', moviedir='movies/brawley/'):
	todt=dts[1]
	dt0=catdate
	# 2012/7/1 3:25:20.450000	32.2072	-115.3432	4.64	3.18

	quakes = []
	#quakes += [[dts[0], dts[1], -115.3432, 32.2072, 4.64]]
	#
	hm=generalHazMovie(catname=catname, catnum=catnum, mc=mc, gridsize=gridsize, ndithers=ndithers, winlen=winlen, sigma=sigma, avlen=avlen, bigmag=bigmag, fignum=fignum, dts=dts, catdate=catdate, epicen=epicen, thresh=thresh, refreshcat=refreshcat, lats=lats, lons=lons, logZ=logZ, nameroot=nameroot, moviedir=moviedir, events=quakes)
	#
	#os.system('avconv -i %s/contour/%s-%05d.png -r 50 %s/%s.mp4', (moviedir, nameroot, moviedir,nameroot))
	#
	return hm
#
# lats=[35.3, 36.5], lons=[-121.1, -119.9]
def ParkfieldMovie(catname='cats/pfmovie.cat', catnum=0, mc=1.5, gridsize=None, ndithers=10, winlen=30, sigma=1.68, avlen=None, bigmag=4.9, fignum=0, dts=[dtm.datetime(2001, 6, 1, 0, 0, 0, 0, pytz.timezone('UTC')), dtm.datetime(2007, 6, 1, 0, 0, 0, 0, pytz.timezone('UTC')), dtm.timedelta(days=1)], catdate=dtm.datetime(1990, 1, 1, 0, 0, 0, 0, pytz.timezone('UTC')), epicen=None, thresh=None, refreshcat=True, lats=[35.3, 36.5], lons=[-121.1, -119.9], logZ=None, nameroot='parkfield', moviedir='movies/parkfield/', framespeed=50):
	todt=dts[1]
	dt0=catdate
	# 2012/7/1 3:25:20.450000	32.2072	-115.3432	4.64	3.18
	#
	quakes = []
	quakes += [[dts[0], dts[1], -120.37, 35.81, 6.0]]
	#
	hm=generalHazMovie(catname=catname, catnum=catnum, mc=mc, gridsize=gridsize, ndithers=ndithers, winlen=winlen, sigma=sigma, avlen=avlen, bigmag=bigmag, fignum=fignum, dts=dts, catdate=catdate, epicen=epicen, thresh=thresh, refreshcat=refreshcat, lats=lats, lons=lons, logZ=logZ, nameroot=nameroot, moviedir=moviedir, events=quakes)
	#
	os.system('avconv -i %s/contour/%s-\%05d.png -r %d %s/%saconts.mp4' % (moviedir, nameroot, framespeed, moviedir, movieroot))
	#
	return hm
#
def elmayorAftershocks(mshock=5.0, refreshcat=True,lons=[-116.5, -114.5], lats=[31.75, 32.75], todt=dtm.datetime(2010,4,6, tzinfo=pytz.timezone('UTC')),catname='cats/elmayorshocks.cat'):
	# 
	# first, get a current map and all the aftershocks m>4.5?
	a=generalHazmap(lons=lons, lats=lats, mc=2.25, todt=todt, refreshcat=refreshcat, catname=catname, ndithers=10)
	#
	c1=a.hazmaps[0].catalog
	c2=eqp.eqcatalog()
	c2.loadCatFromFile(catname)
	shocks=[]
	for rw in c2.getcat(0):
		if rw[3]>=mshock and rw[0]>dtm.datetime(2010,4,3, 22, 40, 43, tzinfo=pytz.timezone('UTC')) and rw[0]<dtm.datetime.now(pytz.timezone('UTC')) :
			shocks+=[rw]
		#
	#
	plt.figure(2)
	for rw in shocks:
		#print rw, str(rw[0])
		x,y=c1.catmap(rw[2], rw[1])
		plt.plot([x], [y], '*', ms=12, label='$%.2f: %s' % (rw[3], str(rw[0])), zorder=17 )
	#
	plt.legend(loc=0)
	return a
	
#
def elmayorminimovie(catname='cats/elmayormini.cat', catnum=0, mc=2.25, gridsize=None, ndithers=10, winlen=20, sigma=1.68, avlen=None, bigmag=7.9, fignum=0, dts=[dtm.datetime(2012, 1, 1, 0, 0, 0, 0, pytz.timezone('UTC')), dtm.datetime.now(pytz.timezone('UTC')), dtm.timedelta(days=1)], catdate=dtm.datetime(1990, 1, 1, 0, 0, 0, 0, pytz.timezone('UTC')), epicen=None, thresh=None, refreshcat=True, lats=[32.0,32.5], lons=[-115.75, -114.8], logZ=None, nameroot='elmayormini', moviedir='movies/elmayormini/', framespeed=50, frame0=0):
	#
	if len(dts)==1: dts+=[dtm.datetime.now(pytz.timezone('UTC'))]
	if len(dts)==2: dts+=[dtm.timedelta(days=1)]
	todt=dts[1]
	dt0=catdate
	# 2012/7/1 3:25:20.450000	32.2072	-115.3432	4.64	3.18

	quakes = []
	#quakes += [[dts[0], dts[1], -115.3432, 32.2072, 4.64]]
	#
	hm=generalHazMovie(catname=catname, catnum=catnum, mc=mc, gridsize=gridsize, ndithers=ndithers, winlen=winlen, sigma=sigma, avlen=avlen, bigmag=bigmag, fignum=fignum, dts=dts, catdate=catdate, epicen=epicen, thresh=thresh, refreshcat=refreshcat, lats=lats, lons=lons, logZ=logZ, nameroot=nameroot, moviedir=moviedir, events=quakes, fps=framespeed, frame0=frame0,movieAnnotator=elmayorminiAnnoator)
	#
	#os.system('avconv -i %s/contour/%s-%s.png -r %d %s/%s.mp4' % (moviedir, nameroot, '%05d', framespeed, moviedir, nameroot))
	#
	return hm
#
def elmayorminiAnnoator(thisdt=None, mycat=None):
		# for now, let's do an example:
		# we want to show a precursory earthquake...
		#eqdt1=dtm.datetime(2012,7,1,6,36,7, int(.3*10**6))
		#eqdt2=dtm.datetime(2012,7,1,3,25,20, int(.8*10**6))
		eqdt0=dtm.datetime(2012,7,1,3,25, tzinfo=pytz.timezone('UTC'))
		eqdt1=dtm.datetime.now(pytz.timezone('UTC'))
		animation_time=dtm.timedelta(days=100)
		#
		if thisdt>=(eqdt0-animation_time) and thisdt<=(eqdt0+animation_time):
			feqdt0=mpd.date2num(eqdt0)
			fthisdt=mpd.date2num(thisdt)
			#
			L=5623.4*1.5	# approx L(m=5) in meters.
			x1,y1=mycat.catmap(-115.313, 32.193)
			x2,y2=mycat.catmap(-115.3, 32.205)
			x0=.5*(x1+x2) + L/4.0
			y0=.5*(y1+y2) + L/4.0	# customizing this a bit for style...
			#
			myalpha=.05+0.0*(1.0-(abs(feqdt0-fthisdt)/float(animation_time.days)))
			circ=plt.Circle((x0, y0), L, color='b', alpha=myalpha, zorder=11)
			ring=plt.Circle((x0, y0), L, color='b', alpha=.65, zorder=11, fill=False, lw=2.0)
			#
			myfig=plt.gcf()
			myfig.gca().add_artist(circ)
			myfig.gca().add_artist(ring)
			plt.plot([x1], [y1], '*', ms=15, alpha=.75, zorder=12)
			plt.plot([x2], [y2], '*', ms=15, alpha=.75, zorder=13)
			plt.show()
		#
		#
		if thisdt>=(eqdt1-animation_time):
			# developing earthquake?
			feqdt1=mpd.date2num(eqdt1)
			fthisdt=mpd.date2num(thisdt)
			L=5623.4*1.5	# approx L(m=5) in meters.
			#
			myalpha=.05+0.0*(1.0-(abs(feqdt1-fthisdt)/float(animation_time.days)))
			x0=-(115.0 + 19./60. + 36.7/3600.0)
			y0 = (32.0 + 15.0/60. + 52.66/3600.)
			x0,y0=mycat.catmap(x0,y0)
			
			circ=plt.Circle((x0, y0), L, color='b', alpha=myalpha, zorder=11)
			ring=plt.Circle((x0, y0), L, color='b', alpha=.75, zorder=11, fill=False, lw=3.0)
			myfig=plt.gcf()
			myfig.gca().add_artist(circ)
			myfig.gca().add_artist(ring)
			plt.show()				
#
def northElMayorMovie(catname='cats/northElMayor.cat', catnum=0, mc=2.25, gridsize=None, ndithers=10, winlen=20, sigma=1.68, avlen=None, bigmag=7.9, fignum=0, dts=[dtm.datetime(2010, 1, 1, 0, 0, 0, 0, pytz.timezone('UTC')), dtm.datetime.now(pytz.timezone('UTC')), dtm.timedelta(days=1)], catdate=dtm.datetime(1990, 1, 1, 0, 0, 0, 0, pytz.timezone('UTC')), epicen=None, thresh=None, refreshcat=True, lats=[32.6, 33.2], lons=[-116.5, -115.8], logZ=None, nameroot='northElMayor', moviedir='movies/northElMayor/', framespeed=50):
	todt=dts[1]
	dt0=catdate
	# 2012/7/1 3:25:20.450000	32.2072	-115.3432	4.64	3.18

	quakes = []
	#quakes += [[dts[0], dts[1], -115.3432, 32.2072, 4.64]]
	#
	hm=generalHazMovie(catname=catname, catnum=catnum, mc=mc, gridsize=gridsize, ndithers=ndithers, winlen=winlen, sigma=sigma, avlen=avlen, bigmag=bigmag, fignum=fignum, dts=dts, catdate=catdate, epicen=epicen, thresh=thresh, refreshcat=refreshcat, lats=lats, lons=lons, logZ=logZ, nameroot=nameroot, moviedir=moviedir, events=quakes)
	#
	os.system('avconv -i %s/contour/%s-\%05d.png -r %d %s/%saconts.mp4' % (moviedir, nameroot, framespeed, moviedir, nameroot))
	#
	return hm
#
def TohokuMovie(catname='cats/tohoku.cat', catnum=0, mc=4.75, gridsize=None, ndithers=15, winlen=20, sigma=1.68, avlen=None, bigmag=7.3, fignum=0, dts=[dtm.datetime(2009, 3, 1, 0, 0, 0, 0, pytz.timezone('UTC')), dtm.datetime.now(pytz.timezone('UTC')), dtm.timedelta(days=1)], catdate=dtm.datetime(1990, 1, 1, 0, 0, 0, 0, pytz.timezone('UTC')), epicen=None, thresh=None, refreshcat=True, lons=[135., 146.], lats=[30., 41.5], logZ=None, nameroot='tohoku', moviedir='movies/tohoku/', framespeed=50, conts=None):
	todt=dts[1]
	dt0=catdate
	# 2012/7/1 3:25:20.450000	32.2072	-115.3432	4.64	3.18

	quakes = []
	#quakes += [[dts[0], dts[1], -115.3432, 32.2072, 4.64]]
	#
	hm=generalHazMovie(catname=catname, catnum=catnum, mc=mc, gridsize=gridsize, ndithers=ndithers, winlen=winlen, sigma=sigma, avlen=avlen, bigmag=bigmag, fignum=fignum, dts=dts, catdate=catdate, epicen=epicen, thresh=thresh, refreshcat=refreshcat, lats=lats, lons=lons, logZ=logZ, nameroot=nameroot, moviedir=moviedir, events=quakes, conts=conts)
	#
	#os.system('avconv -i %s/contour/%s-\%05d.png -r %d %s/%saconts.mp4' % (moviedir, nameroot, framespeed, moviedir, movieroot))
	#
	return hm
#
# 2012/7/1 3:25:20.450000	32.2072	-115.3432	4.64	3.18
def mexiswarm(catname='cats/mexiswarm2012.cat', catnum=0, mc=2.25, gridsize=None, ndithers=10, winlen=20, sigma=1.68, avlen=None, bigmag=7.9, fignum=0, dts=[dtm.datetime(2011, 7, 1, 0, 0, 0, 0, pytz.timezone('UTC')), dtm.datetime.now(pytz.timezone('UTC')), dtm.timedelta(days=1)], catdate=dtm.datetime(1990, 1, 1, 0, 0, 0, 0, pytz.timezone('UTC')), epicen=None, thresh=None, refreshcat=True, lats=[31.9, 32.55], lons=[-115.6, -115.0], logZ=None, nameroot='mexi20120701', moviedir='movies/mexiswarm20120701/'):
	todt=dts[1]
	dt0=catdate
	# 2012/7/1 3:25:20.450000	32.2072	-115.3432	4.64	3.18

	quakes = []
	quakes += [[dts[0], dts[1], -115.3432, 32.2072, 4.64]]
	#
	hm=generalHazMovie(catname=catname, catnum=catnum, mc=mc, gridsize=gridsize, ndithers=ndithers, winlen=winlen, sigma=sigma, avlen=avlen, bigmag=bigmag, fignum=fignum, dts=dts, catdate=catdate, epicen=epicen, thresh=thresh, refreshcat=refreshcat, lats=lats, lons=lons, logZ=logZ, nameroot=nameroot, moviedir=moviedir, events=quakes)
	#
	os.system('avconv -i movies/mexiswarm20120701/contour/mexi20120701-%05d.png -r 50 movies/mexiswarm20120701/mexiswarmconts.mp4')
	#
	return hm
	

def socalswarm(catname='cats/scswarm2011.cat', catnum=0, mc=2.5, gridsize=None, ndithers=10, winlen=20, sigma=1.68, avlen=None, bigmag=7.9, fignum=0, dts=[dtm.datetime(2011, 2, 1, 0, 0, 0, 0, pytz.timezone('UTC')), dtm.datetime.now(pytz.timezone('UTC')), dtm.timedelta(days=1)], catdate=dtm.datetime(1990, 1, 1, 0, 0, 0, 0, pytz.timezone('UTC')), epicen=None, thresh=None, refreshcat=False, lats=[36.25, 36.75], lons=[-121.4, -120.9], logZ=None, nameroot='socal20110827', moviedir='movies/scswarm20110827b/'):
	#
	todt=dts[1]
	dt0=catdate
	#
	# base map should be like:
	# z=mhp.generalHazmap(catname='cats/socalswarm.cat', mc=2.5, lats=[36.25,36.75], lons=[-121.4, -120.9], refreshcat=True)

	#
	#print dts
	#print catdate
	#print lats
	#print lons
	#
	# quakes of interest are:
	#2011/8/27 7:18:21.150000	36.5843	-121.1808	4.64	7.84
	#2011/8/27 7:21:59.970000	36.6012	-121.2025	3.61	6.93
	#2011/7/6 7:18:52.360000	36.6652	-121.2928	3.8	7.59
	quakes = []
	quakes += [[dts[0], dts[1], -121.1808, 36.5843, 4.64]]
	quakes += [[dts[0], dts[1], -121.2025, 36.6012, 3.61]]
	quakes += [[dts[0], dts[1], -121.2928, 36.6652, 3.8]]
	#
	hm=generalHazMovie(catname=catname, catnum=catnum, mc=mc, gridsize=gridsize, ndithers=ndithers, winlen=winlen, sigma=sigma, avlen=avlen, bigmag=bigmag, fignum=fignum, dts=dts, catdate=catdate, epicen=epicen, thresh=thresh, refreshcat=refreshcat, lats=lats, lons=lons, logZ=logZ, nameroot=nameroot, moviedir=moviedir, events=quakes)
	#
	os.system('avconv -i movies/scswarm20110827/contour/socal20110827-%05d.png -r 50 movies/scswarm20110827/scswarmconts.mp4')
	#
	return hm

# 2011-8-27 (norish)socal "swarm":
def generalHazMovie(catname='cats/scswarm2011.cat', catnum=0, mc=4.75, gridsize=None, ndithers=4, winlen=20, sigma=1.68, avlen=None, bigmag=7.9, fignum=0, dts=[dtm.datetime(2010, 1, 1, 0, 0, 0, 0, pytz.timezone('UTC')), dtm.datetime.now(pytz.timezone('UTC')), dtm.timedelta(days=1)], catdate=dtm.datetime(1990, 1, 1, 0, 0, 0, 0, pytz.timezone('UTC')), epicen=None, thresh=None, refreshcat=False, lats=[-9.0, 10.0], lons=[92.0, 106.0], logZ=None, nameroot='hmmovie', moviedir='movies/', events=[], conts=None, fps=50, frame0=0, movieAnnotator=None):	
	#
	if len(dts)==1: dts+=[dtm.datetime.now(pytz.timezone('UTC'))]
	if len(dts)==2: dts+=[dtm.timedelta(days=1)]
	#
	dt0=catdate		# starting date for catalog.
	startdate=dts[0]	# forecast start date (first map).
	enddate=dts[1]
	todt=enddate
	deltat=dts[2]
	#
	# conts: could be an integer (number of contours), or an array of floats (contours to draw).
	if conts==None: conts=[-.9, -.7, -.5, -.3, -.1, 0.0, .1, .2, .3, .4, .5, .6, .7] 
	#
	# events: [startshow, stopshow, lon, lat, mag]
	#
	# eventually, build in directory-error-control here...
	#
	# (catname='cats/indonesia2012.cat', catnum=0, mc=4.75, gridsize=None, ndithers=4, winlen=20, sigma=1.68, avlen=None, bigmag=7.9, fignum=0, todt=dtm.datetime.now(pytz.timezone('UTC')), epicen=None, thresh=None, refreshcat=False, lats=[-9.0, 10.0], lons=[92.0, 106.0], dt0=None, logZ=None)
	# note: this will (if refreshing) fetch the full catalog from catdate -> todt/enddate
	objHM=generalHazmap(catname=catname, catnum=catnum, mc=mc, gridsize=gridsize, ndithers=ndithers, winlen=winlen, sigma=sigma, avlen=avlen, bigmag=bigmag, fignum=fignum, todt=todt, epicen=epicen, thresh=thresh, refreshcat=refreshcat, lats=lats, lons=lons, dt0=dt0, logZ=logZ)
	#return objHM
	#
	#fignums=[fignum, fignum+2, fignum+3]
	#figdirs=[moviedir+'binary/', moviedir+'contour/']
	#fignums=[fignum, fignum+3]
	fignums=[fignum+3]
	#figdirs=[moviedir+'binary/', moviedir+'contour/']
	figdirs=[moviedir+'contour/']
	#
	#drs = glob.glob(moviedir)
	for dr in figdirs:
		# does this folder exist?
		figdirses=dr.split('/')
		currentdir=''
		for thisdir in figdirses:
			a=glob.glob('*')
			currentdir = currentdir + thisdir + '/'
			if currentdir not in a:
				os.system('mkdir %s' % currentdir)
	#
	thisdt=startdate
	#framenum=0
	framenum=frame0
	framenumstr=('00000' + str(framenum))[-5:]
	print "initial frame-str: %s" % framenumstr
	while thisdt<dts[1]:
		#
		objHM.hazMapTo(thisdt)
		#objHM.simpleBoxes(fignum=fignum, thresh=-0.0, plotevents=0)	
		objHM.contourMap(fignum=fignum+2, nconts=conts)
		#
		i=0
		#
		mycat=objHM.hazmaps[0].catalog
		for fn in fignums:
			myfig=plt.figure(fn)
			for ev in events:
				if thisdt>=ev[0] and thisdt<=ev[1]:
					#x,y=objHM.hazmaps[0].catalog.catmap(ev[2], ev[3])
					x,y=mycat.catmap(ev[2], ev[3])
					plt.plot([x], [y], 'r*', ms=20, alpha=.75, zorder=15)
					#
			#
			# what we should do her is drop in a call to some customizatization function (passed as a parameter).
			#
			if movieAnnotator!=None: movieAnnotator(thisdt, mycat)		# annotation stuff. we'll need to know the date
																							# and have a pointer to the catalog so we can
																							# draw stuff
			#
			#if fn==fignums[-1]: continue	# we only want one contour plot
			#
			# save plots:
			#print "saving as: %s%s-%s" % (figdirs[i], nameroot, framenumstr)
			fnameroot='%s%s-%s' % (figdirs[i], nameroot, framenumstr)
			plt.savefig('%s.png' % fnameroot)
			#os.system('convert %s.png %s.jpg' % (fnameroot, fnameroot))
			i+=1
		#
		thisdt+=deltat
		framenum+=1
		framenumstr='00000' + str(framenum)
		framenumstr=framenumstr[-5:]
	#
	os.system('rm %s/%s.mp4' % (moviedir,nameroot))
	avconvstr='avconv -i %s/contour/%s-%s.png -r %d %s/%s.mp4' % (moviedir, nameroot, '%05d', fps, moviedir,nameroot)
	os.system(avconvstr)
	return objHM


def japanmovie(catname='cats/japancat4b.cat', catnum=0, mc=4.5, gridsize=None, ndithers=3, winlen=20, sigma=1.68, avlen=2, bigmag=7.0, fignum=0, dts=[dtm.datetime(2010, 1, 1, 0, 0, 0, 0, pytz.timezone('UTC')), dtm.datetime.now(pytz.timezone('UTC')), dtm.timedelta(days=1)], epicen=None, moviedir='movies/japan1/', framenum=0):
	objHM=makeJapanHazMapDith(catname=catname, catnum=catnum, mc=mc, gridsize=gridsize, ndithers=ndithers, winlen=winlen, sigma=sigma, avlen=avlen, bigmag=bigmag, fignum=fignum, dateto=dts[0], epicen=epicen)
	evll=[142.37, 38.30]	# event lat/lon
	#
	fignums=[fignum, fignum+2, fignum+3]
	figdirs=[moviedir+'binary/', moviedir+'contour/']
	#
	deltat=dts[2]
	thisdt=dts[0]
	#framenum=0
	framenumstr='00000'
	while thisdt<dts[1]:
		#
		objHM.hazMapTo(thisdt)
		objHM.simpleBoxes(fignum=fignum, thresh=-0.0, plotevents=0)	
		objHM.contourMap(fignum=fignum+2)
		#
		x,y=objHM.hazmaps[0].catalog.catmap(evll[0], evll[1])	# note, we have to do this every time because we can phase shift with the catalog.
		#
		i=0
		if thisdt>dtm.datetime(2011,1,1, 0, 0, 0, 0, pytz.timezone('UTC')): deltat=dtm.timedelta(hours=12)
		if thisdt>dtm.datetime(2011,2,15, 0, 0, 0, 0, pytz.timezone('UTC')): deltat=dtm.timedelta(hours=6)
		if thisdt>dtm.datetime(2011,4,15, 0, 0, 0, 0, pytz.timezone('UTC')): deltat=dtm.timedelta(hours=12)
		if thisdt>dtm.datetime(2011,6,1, 0, 0, 0, 0, pytz.timezone('UTC')):deltat=dtm.timedelta(days=1)
		#
		for fn in fignums:
			plt.figure(fn)
			plt.plot([x], [y], 'r*', ms=20)
			#
			if fn==fignums[-1]: continue	# we only want one contour plot
			#
			# save plots:
			fnameroot='%sjapanhm-%s.png' % (figdirs[i], framenumstr)
			plt.savefig('%s.png' % fnameroot)
			os.system('convert %s.png %s.jpg' % (fnameroot, fnameroot))
			i+=1
		#
		thisdt+=deltat
		framenum+=1
		framenumstr='00000' + str(framenum)
		framenumstr=framenumstr[-5:]
	#
	#for fd in figdirs:
		# something like:
		# mencoder mf://*.jpg -mf w=800:h=600:fps=ca25:type=jpg -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o output.avi
		# as per ( http://www.mplayerhq.hu/DOCS/HTML/en/menc-feat-enc-images.html )

		#os.system('mencoder mf://*.jpg -mf w=800:h=600:fps=15:type=jpg -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o japan.avi')
		#os.system('mencoder mf://%s*.jpg -mf w=800:h=600:fps=15:type=jpg -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o %s/japan.avi' % (fd,fd))
		# (use avconv...)
	#
	return objHM

class movieThread(threading.Thread):
	def __init__(self, catname='cats/mexicat-short.cat', catnum=0, mc=2.5, gridsize=None, ndithers=5, winlen=20, sigma=1.68, avlen=None, bigmag=5.0, fignum=0, dts=[dtm.datetime(2008, 1, 1, 0, 0, 0, 0, pytz.timezone('UTC')), dtm.datetime(2010,12,31, 0, 0, 0, 0, pytz.timezone('UTC')), dtm.timedelta(days=1)], epicen=None, moviedir='movies/mexi1/', framenum=0, evll=[-115.3, 32.13], logZ=None, nameroot='mexi', dN=1):
		threading.Thread.__init__(self)
		self.catname=catname
		self.catnum=catnum
		self.mc=mc
		self.gridsize=gridsize
		self.ndithers=ndithers
		self.winlen=winlen
		self.sigma=sigma
		self.avlen=avlen
		self.bigmag=bigmag
		self.fignum=fignum
		self.dts=dts
		self.epicen=epicen
		self.moviedir=moviedir
		self.framenum=framenum
		self.evll=evll
		self.logZ=logZ
		self.nameroot=nameroot
		self.dN=dN
		
	def run(self):
		a=meximovie(catname=self.catname, 		catnum=self.catnum, mc=self.mc, gridsize=self.gridsize, ndithers=self.ndithers, winlen=self.winlen, sigma=self.sigma, avlen=self.avlen, bigmag=self.bigmag, fignum=self.fignum, dts=self.dts, epicen=self.epicen, moviedir=self.moviedir, framenum=self.framenum, evll=self.evll, logZ=self.logZ, nameroot=self.nameroot, dN=self.dN)

def threadedMovie(catname='cats/mexicat-short.cat', catnum=0, mc=2.5, gridsize=None, ndithers=5, winlen=20, sigma=1.68, avlen=None, bigmag=5.0, fignum=0, dts=[dtm.datetime(2008, 1, 1, 0, 0, 0, 0, pytz.timezone('UTC')), dtm.datetime(2010,12,31, 0, 0, 0, 0, pytz.timezone('UTC')), dtm.timedelta(days=1)], epicen=None, moviedir='movies/mexi1/', framenum=0, evll=[-115.3, 32.13], logZ=None, nameroot='mexi', dN=1):
	# and this may not be working quite right just yet. there appears to be a problem separating the tkinter/pyplot
	# instances. we get a "main thred not in main loop" error.
	# dN will be number of threads.
	#moviethreads=[]
	for i in xrange(dN):
		#moviethreads+=
		t=movieThread(catname=catname, catnum=catnum, mc=mc, gridsize=gridsize, ndithers=ndithers, winlen=winlen, sigma=sigma, avlen=avlen, bigmag=bigmag, fignum=fignum, dts=dts, epicen=epicen, moviedir=moviedir, framenum=framenum + i, evll=evll, logZ=logZ, nameroot=nameroot, dN=dN)
		t.start()
		

def meximovie(catname='cats/mexicat-short.cat', catnum=0, mc=2.5, gridsize=None, ndithers=5, winlen=20, sigma=1.68, avlen=None, bigmag=5.0, fignum=0, dts=[dtm.datetime(2008, 1, 1, 0, 0, 0, 0, pytz.timezone('UTC')), dtm.datetime(2010,12,31, 0, 0, 0, 0, pytz.timezone('UTC')), dtm.timedelta(days=1)], epicen=None, moviedir='movies/mexi1/', framenum=0, evll=[-115.3, 32.13], logZ=1.0, nameroot='mexi', dN=1):
	# dN is the incrementation for frame naming, aka 0001, 0002, ... as opposed to maybe threading the process, in which
	# each thread starts with a \phi offset and increments in +phi.
	objHM=makeJapanHazMapDith(catname=catname, catnum=catnum, mc=mc, gridsize=gridsize, ndithers=ndithers, winlen=winlen, sigma=sigma, avlen=avlen, bigmag=bigmag, fignum=fignum, dateto=dts[0], epicen=epicen, logZ=logZ)
	#evll=[-115.3, 32.13]	# event lat/lon
	#
	fignums=[fignum, fignum+2, fignum+3]
	figdirs=[moviedir+'binary/', moviedir+'contour/', moviedir+'contour2/' ]
	#
	deltat=dts[2]
	thisdt=dts[0]
	#framenum=0
	framenumstr='00000'
	while thisdt<dts[1]:
		objHM.hazMapTo(thisdt)
		# set a max/min site so contour colors are consistent:
		objHM.hazmaps[0].rbRatios[0]=0.8
		objHM.hazmaps[1].rbRatios[0]=-1.0
		objHM.simpleBoxes(fignum=fignum, thresh=None, plotevents=0)	
		objHM.contourMap(fignum=fignum+2)
		#
		x,y=objHM.hazmaps[0].catalog.catmap(evll[0], evll[1])	# note, we have to do this every time because we can phase shift with the catalog.
		#
		i=0
		if thisdt>dtm.datetime(2009,9,1, 0, 0, 0, 0, pytz.timezone('UTC')): deltat=dtm.timedelta(hours=12)
		if thisdt>dtm.datetime(2010,2,15, 0, 0, 0, 0, pytz.timezone('UTC')): deltat=dtm.timedelta(hours=4)
		if thisdt>dtm.datetime(2010,5,15, 0, 0, 0, 0, pytz.timezone('UTC')): deltat=dtm.timedelta(hours=12)
		if thisdt>dtm.datetime(2010,7,1, 0, 0, 0, 0, pytz.timezone('UTC')):deltat=dtm.timedelta(days=1)
		#
		for fn in fignums:
			plt.figure(fn)
			plt.plot([x], [y], 'r*', ms=20)
			plt.title('%d-%d-%d, $d\\lambda=%.4f$, $nDith=%d$, $mc=%.2f$, $rblen=%d$, $avlen=%d$\n\n\n' % (thisdt.year, thisdt.month, thisdt.day, objHM.gridsize, objHM.ndithers, objHM.mc, objHM.winlen, objHM.avlen))
			#
			#if fn==fignums[-1]: continue	# we only want one contour plot
			#
			# save plots:
			fnameroot='%s%s-%s' % (figdirs[i], nameroot, framenumstr)
			plt.savefig('%s.png' % fnameroot)
			#os.system('convert %s.png %s.jpg' % (fnameroot, fnameroot))
			i+=dN
		#
		thisdt+=deltat
		framenum+=1
		framenumstr='00000' + str(framenum)
		framenumstr=framenumstr[-5:]
	#
	# and we've found a better way to do this with av-lib or something... (new fencoder or something)
	# i've been using this command-line:
	# avconv -i contour2/salton-%05d.png -r 50 saltonconts1.mp4
	#for fd in figdirs:
		# something like:
		# mencoder mf://*.jpg -mf w=800:h=600:fps=ca25:type=jpg -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o output.avi
		# as per ( http://www.mplayerhq.hu/DOCS/HTML/en/menc-feat-enc-images.html )

		#os.system('mencoder mf://*.jpg -mf w=800:h=600:fps=15:type=jpg -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o japan.avi')
		#os.system('mencoder mf://%s*.jpg -mf w=800:h=600:fps=15:type=jpg -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o %s/mexi.avi' % (fd,fd))
		# and this should work (note, the encoding is implied by the suffix):
		# avconv -i contour2/salton-%05d.png.png -r 50 saltonconts1.mp4
	#
	return objHM

def makeJapCat(lon=[135., 146.], lat=[30., 41.5], minMag=4.5, dates0=[dtm.datetime(2000,1,1, tzinfo=pytz.timezone('UTC')), None], Nmax=999999, fout='cats/japcat0.cat'):
	# get a basic catalog. then, we'll do a poly-subcat. we need a consistent catalog.
	
	if dates0[1]==None:
		# i think this needs a "date" object, and datetime breaks.
		# so, make a Now() for date.
		nowdtm=dtm.datetime.now(pytz.timezone('UTC'))
		#dates0[1]=dtm.date(nowdtm.year, nowdtm.month, nowdtm.day)
		dates0[1]=nowdtm
	#	
	'''
	catlist=eqp.getANSSlist(lon, lat, minMag, dates0, Nmax, None)
	f=open(fout, 'w')
	f.write("#anss catalog\n")
	f.write("#lon=%s\tlat=%s\tm0=%f\tdates=%s\n" % (str(lon), str(lat), minMag, str(dates0)))
	for rw in catlist:
		f.write('%s\t%s\t%s\t%s\n' % (rw[0], rw[1], rw[2], rw[4]))
	f.close()
	'''
	catlist=atp.catfromANSS(lon, lat, minMag, dates0, Nmax, fout)

def makeMexicat(lats=[31.5, 33.0], lons=[-116., -114.75], mc=2.5, dts=[dtm.datetime(1990,1,1, 0, 0, 0, 0, pytz.timezone('UTC')), None], fout='cats/mexicat0.cat'):
	# catfromANSS(lon=[135., 150.], lat=[30., 41.5], minMag=4.0, dates0=[dtm.date(2005,01,01), None], Nmax=999999, fout='cats/mycat.cat')
	return atp.catfromANSS(lat=lats, lon=lons, minMag=mc, dates0=dts, fout=fout)

def makeSocalcat(lat=[31., 37.], lon=[-119., -114.], minMag=2.25, dates0=[dtm.datetime(1990,1,1, 0, 0, 0, 0, pytz.timezone('UTC')), dtm.datetime.now(pytz.timezone('UTC'))], fout='cats/socal.cat'):
	a=atp.catfromANSS(lat=lat, lon=lon, minMag=minMag, dates0=dates0, fout=fout)
	return a

def makelaquilacat(lat=[40., 45.], lon=[11.,18.], mc=1.5, dts=[dtm.datetime(2000,1,1, 0, 0, 0, 0, pytz.timezone('UTC')), None], fout='cats/laquila.cat'):
	a=atp.catfromANSS(lat=lat, lon=lon, minMag=mc, dates0=dts, fout=fout)
	return a


def makeMeximap(catname='cats/mexicat-short.cat', catnum=0, mc=2.5, gridsize=None, ndithers=4, winlen=20, sigma=1.68, avlen=None, bigmag=6.5, fignum=0, todt=dtm.datetime(2010, 4, 1, 0, 0, 0, 0, pytz.timezone('UTC')), thresh=0):
	if avlen==None and winlen==None: avlen=1
	if avlen==None:
		avlen=int(math.ceil(winlen/10.))
	#
	scm = makeSocalMap(catname=catname, catnum=catnum, mc=mc, gridsize=gridsize, ndithers=ndithers, winlen=winlen, sigma=sigma, avlen=avlen, bigmag=bigmag, fignum=fignum, todt=todt, epicen=[-115.3, 32.13], thresh=thresh)
	# el mayor (baja, mexicali) earthquake:
	plt.figure(fignum)
	#hmx, hmy=scm.hazmaps[0].catalog.catmap(-115.3, 32.13)
	#hmx, hmy=scm.hazmaps[0].catalog.catmap(epicen[0], epicen[1])
	for fnum in [0,2,3]:
		plt.figure(fignum+fnum)
		# print "plotting mev on figure %d" % (fnum+fignum)
		scm.plotEvent([-115.3, 32.13], fignum=(fignum+fnum), symbolstr='k*', msize=24)
		scm.plotEvent([-115.3, 32.13], fignum=(fignum+fnum), symbolstr='r*', msize=20)
		
		#plt.plot([hmx], [hmy], 'b*', ms=22, zorder=10, alpha=.8)
		#plt.plot([hmx], [hmy], 'r*', ms=18, zorder=10, alpha=.8)
		#scm.plotEvent([-115, 32.], fignum=(fignum+fnum), symbolstr='ok', msize=24)
	plt.figure(fignum)
	#
	# and get a sampling of before/after earthquakes:
	thiscat=scm.hazmaps[0].catalog.getcat(0)
	catmag=mc+.5	# just to thin it out...
	maxnums=[500,1000]	# [before, after]
	beforecat=[]
	aftercat=[]
	# find mainshock...
	mevindex=0
	for rw in thiscat:
		if mpd.date2num(rw[0])>mpd.date2num(dtm.datetime(2010,4,1, 0, 0, 0, 0, pytz.timezone('UTC'))) and rw[3]>=7.0: break
		mevindex+=1
	i=0
	while len(beforecat)<maxnums[0] and (mevindex-i)>=0:
		#beforecat+=[thiscat[mevindex-maxnums[0]+i]]
		if thiscat[mevindex-i][3]>=catmag: beforecat.insert(0, thiscat[mevindex-i])
		i+=1
	i=0
	while len(aftercat)<maxnums[1] and (mevindex+i)<len(thiscat):
		if thiscat[mevindex+i][3]>=catmag: aftercat+=[thiscat[mevindex+i]]
		i+=1
	#
	scm.plotEvents(events=beforecat, symbolstr='r.', msize=6, fignum=fignum+3, alpha=.35, zorder=5)
	scm.plotEvents(events=aftercat, symbolstr='bs', msize=3, fignum=fignum+3, alpha=.15, zorder=5)
	# nsquares=None, mc=None, deltas=2.2, beta=2.0, rblen=None, dithers=None)
	magres=scm.guessMag(nsquares=1, mc=mc, deltas=2.0, beta=1.0, rblen=winlen, dithers=ndithers**2)
	magresDith=scm.guessMag(nsquares=ndithers**2, mc=mc, deltas=2.0, beta=1.0, rblen=winlen, dithers=ndithers**2)
	print "mag-res: %f/%f (dithering=%d)" % (magres, magresDith, ndithers)
	
	return scm

def makeSocalMap(catname='cats/socal.cat', catnum=0, mc=2.5, gridsize=None, ndithers=4, winlen=20, sigma=1.68, avlen=2, bigmag=6.5, fignum=0, todt=dtm.datetime(2010, 4, 1, 0, 0, 0, 0, pytz.timezone('UTC')), epicen=None, thresh=0, mapres='i',logZ=None, alpha=.9):
	# one mc, gridsize combination:mc=2.5, gridsize=.125, winlen=20 (approximately)
	# choose an optimal combination:
	# get gridsize/winlen:
	# getgridsize(mc=4.5, winlen=16, lat=45., sigmaExp=1.68)
	# def getWinlen(mc=4.5, gridsize=1.0, lat=45., sigmaExp=1.68):
	#impor
	# get gridsize/winlen:
	# getgridsize(mc=4.5, winlen=16, lat=45., sigmaExp=1.68)
	# def getWinlen(mc=4.5, gridsize=1.0, lat=45., sigmaExp=1.68):
	if gridsize==None and winlen==None:
		# can't do anything with that, can we?
		return None
	midlat=getMidLat(catname, mc)
	if gridsize==None:
		gridsize=getgridsize(mc, winlen, midlat, sigma)
	if winlen==None:
		winlen=getWinlen(mc, gridsize, midlat, sigma)
	#
	#
	#objHM=makeHazMapDith(catname, catnum, mc, gridsize, ndithers, winlen, sigma, avlen, bigmag, fignum)
	objHM=makeHazMapDith(catname, catnum, mc, gridsize, ndithers, winlen, sigma, avlen, bigmag, fignum, logZ=None)
	objHM.epicen=epicen
	# set map resolution. note, we use the catalog object in the first hazard-map object (objHM contains
	# a collection of hazard map objects) to draw maps et al.
	#objHM.hazmaps[0].catalog.mapres=mapres
	#
	objHM.hazMapTo(todt)
	#
	#
	magres=objHM.guessMag(nsquares=1, mc=mc, deltas=2.0, beta=1.0, rblen=winlen, dithers=ndithers**2)
	magresDith=objHM.guessMag(nsquares=ndithers**2, mc=mc, deltas=2.0, beta=1.0, rblen=winlen, dithers=ndithers**2)
	print "mag-res: %f/%f (dithering=%d)" % (magres, magresDith, ndithers)
	#
	plt.figure(fignum)
	objHM.simpleBoxes(fignum=fignum, thresh=thresh, mapres=mapres)
	drawgarnish(cm=objHM.hazmaps[0].catalog.catmap, mapres=mapres)
	#drawgarnish(cm=objHM.hazmaps[0].catalog.catmap, mapres=mapres)
	#
	# el mayor (baja, mexicali) earthquake:
	#hmx, hmy=objHM.hazmaps[0].catalog.catmap(-115.3, 32.13)
	#plt.plot([hmx], [hmy], 'b*', ms=18, zorder=10, alpha=.8)
	
	z2=objHM.contourMap(fignum+2, mapres=mapres, alpha=alpha)
	plt.figure(fignum+1)
	drawgarnish(cm=objHM.hazmaps[0].catalog.catmap, mapres=mapres)
	plt.figure(fignum+2)
	drawgarnish(cm=objHM.hazmaps[0].catalog.catmap, mapres=mapres)
	
	#plt.plot([hmx], [hmy], 'b*', ms=18)
	#
	'''
	for fn in [0,2,3]:
		plt.figure(fn)
		print "plot rivers and stuff %d" % fn
		#objHM.hazmaps[0].catalog.catmap.drawcountries(zorder=1, linewidth=1)
		#objHM.hazmaps[0].catalog.catmap.drawstates(zorder=1, linewidth=1)
		#objHM.hazmaps[0].catalog.catmap.drawmapscale(lon=-115., lat=33., lon0=-115.5, lat0=32.5, length=75., labelstyle='simple', barstyle='simple')
		#objHM.hazmaps[0].catalog.catmap.drawrivers(linewidth=1, color='b')
		drawgarnish(cm=objHM.hazmaps[0].catalog.catmap, mapres=mapres)
	'''
	
	return objHM

def checklognorms():
	a=saltonstuff() # returning the ct catalog object.
	#
	r1=a.rb.getIntervalRatios(minmag=2.25, windowLen=400, cat0=a.getcat(0), logZ=None)
	r2=a.rb.getIntervalRatios(minmag=2.25, windowLen=400, cat0=a.getcat(0), logZ=1.0)
	# so r1: (is normalized) -> averaged
	#    r2: (is raw) -> average -> normalize
	#
	dts=[]
	for rw in r1:
		dts+=[mpd.date2num(rw[1])]
	#
	r1ave=eqp.logaverageOver(map(mhp.operator.itemgetter(2), r1), 40)
	r2ave=eqp.logaverageOver(map(mhp.operator.itemgetter(2), r2), 40)
	
	plt.figure(7)
	plt.ion()
	plt.plot(dts, r1, '-', 'norm')
	plt.plot(dts, r2, '-', label='not-norm')

def pfquad(mc=1.5, targmag=6.0, rbavelen=None, bigmag=5.0, intlist=[25, 50, 100, 200, 300, 400, 500], catname='cats/parkfield.cat', refreshcat=False, plotevents=True,lats=[35.3, 36.5], lons=[-121.1, -119.9], logZ=None, A=.4, B=.15, theta=50):
	plt.ion()
	# catalog.addEllipCat('PFshock (.4x.15)', catalog.cat, pftheta, 35.9, -120.5, majAx, 0.15)
	# for m=9.1, N=223, for m=8.6, N=70 (for mc=4.75)
	# the (a) catalog:
	
	if catname==None:
		refreshcat=True
		catname='cats/parkfield.cat'
	if refreshcat==True:
		#cl1=atp.catfromANSS(lon=[85.0, 105.0],lat=[-10.0, 10.0], minMag=3.5, dates0=[eqp.dtm.date(2000,1,1), eqp.dtm.date(2012,4,14)], fout='cats/indonesia201204.cat')
		#cl1=atp.catfromANSS(lon=[85.0, 105.0],lat=[-10.0, 10.0], minMag=3.5, dates0=[eqp.dtm.date(1990,1,1), eqp.dtm.date(2012,4,14)], fout=catname)
		dLat=.6
		dLon=.6
		#cl1=atp.catfromANSS(lon=[-120.5-dLon, 35.9-dLat],lat=[-120.5+dLon, 35.9+dLat], minMag=mc, dates0=[eqp.dtm.datetime(1990,1,1, tzinfo=pytz.timezone('UTC')), eqp.dtm.datetime.now(pytz.timezone('UTC'))], fout=catname)
		cl1=atp.catfromANSS(lon=lons,lat=lats, minMag=mc, dates0=[eqp.dtm.datetime(1990,1,1, tzinfo=pytz.timezone('UTC')), eqp.dtm.datetime.now(pytz.timezone('UTC'))], fout=catname)
	#
	c1=eqp.eqcatalog([])
	c1.loadCatFromFile(catname)
	#majAx=.4
	majAx=A
	minAx=B
	pftheta=theta
	c1.addEllipCat('PFshock (.4x.15)', c1.cat, pftheta, 35.9, -120.5, majAx, minAx)
	#
	#winlen=int(10**(targmag - (2.0+mc)))
	winlen = Nofm(targmag, mc, mt=7.6, b1=1.0, b2=1.5)*10**(-2.0)	# introducing -(dm+dm_s)
	winlen=int(10*(int(winlen/10)))
	#
	arb=c1.rbomoriQuadPlot(mc=mc, winlen=winlen, catnum=1, rbavelen=rbavelen, bigmag=bigmag, intlist=intlist, plotevents=plotevents, logZ=logZ)
	#
	#return arb
	return c1
	
	

# Indonesia:
def indoQuad(mc=4.75, targmag=9.1, rbavelen=None, bigmag=7.5, intlist=[25, 50, 100, 200, 300, 400, 500], catname='cats/indonesia2012.cat', refreshcat=False, plotevents=False, mt=7.55):
	plt.ion()
	# for m=9.1, N=223, for m=8.6, N=70 (for mc=4.75)
	# the (a) catalog:
	if catname==None:
		refreshcat=True
		catname='cats/indonesia2012.cat'
	if refreshcat==True:
		#cl1=atp.catfromANSS(lon=[92.0, 106.0],lat=[-9.0, 10.0], minMag=3.5, dates0=[eqp.dtm.datetime(1990,1,1, tzinfo=pytz.timezone('UTC')), eqp.dtm.datetime(2012,6,19, tzinfo=pytz.timezone('UTC'))], fout=catname)
		cl1=atp.catfromANSS(lon=[92.0, 106.0],lat=[-9.0, 10.0], minMag=3.5, dates0=[eqp.dtm.datetime(1990,1,1, tzinfo=pytz.timezone('UTC')), eqp.dtm.datetime.now(pytz.timezone('UTC'))], fout=catname)
	#
	c1=eqp.eqcatalog([])
	c1.loadCatFromFile(catname)
	#
	#winlen=int(10**(targmag - (2.0+mc)))
	#winlen = Nofm(m=targmag, mc=mc, mt=7.6, b1=1.0, b2=1.5)*10**(-2.0)	# introducing -(dm+dm_s)
	#winlen=int(10*(int(winlen/10)))
	if rbavelen==None: avlen=getNave(m=targmag, mt=mt, mc=mc)
	winlen=int(getNsample(m=targmag, mc=mc, mt=mt, doint=True))
	print "rbts for winlen=%d" % winlen
	#
	arb=c1.rbomoriQuadPlot(mc=mc, winlen=winlen, rbavelen=avlen, bigmag=bigmag, intlist=intlist, plotevents=plotevents, logZ=None)
	#arb=c1.rbomoriQuadPlot(mc=mc, winlen=winlen, rbavelen=rbavelen, bigmag=bigmag, intlist=intlist, plotevents=plotevents, logZ=None)
	#
	#return arb
	return c1

def Nofm(m, mc, mt=7.6, b1=1.0, b2=1.5):
	# returns the full GR number.
	if m<mt:
		return 10**(b1*(m-mc))
	if m>=mt:
		return 10**(b1*(mt-mc) + b2*(m-mt))
#def Ngr(m, mc, mt=7.6, b1=1.0, b2=1.5):
#	return Nofm(m=m, mc=mc, mt=mt, b1=b1, b2=b2)
getNgr=Nofm
#
def getNomori(m, mc, mt=7.6, dm=1.0, b1=1.0, b2=1.5):
	#dms = getdmSample(m=m, mc=mc, mt=mt, dm=dm, b1=b1, b2=b2)
	if m<mt:
		return 10**(b1*(m-mc-dm))
	if m>=mt:
		return 10**(b2*(m-mt-dm) + b1*(mt-mc))
#
def getNsample(m, mc, mt=7.6, dm=1.0, b1=1.0, b2=1.5, dms0=1.0, doint=True, dmMode=1):
	# dmMode: mode for calculating dm for large earthquakes (see code below).
	targmag=m
	if targmag<mt:
		# "small" earthquake
		winlen=10**(targmag-dm-dms0-mc)	# where 2.0 is dmstar + dmprime
	if targmag>=mt:
		# we want to use an average value for both dms and dm, so dm (like dms) is determined from
		# the full sequence, not simply with respect to the mainshock.
		dmfactor = (1.0*(mt-mc) + 1.5*(targmag-mt))/(targmag-mc)
		#
		#dms = dms0*(1.0*(mt-mc) + dms0*1.5*(targmag-mt))/(targmag-mc)
		dms = dmfactor*dms0
		if dmMode==1: thisdm = dmfactor*dm	# relative to full sequence.
		if dmMode==0: thisdm = 1.5*dm		# relative to largest magnitude
		#winlen = 10**(1.0*(mt-mc) + 1.5*(targmag-mt-1.0) - dms)
		#winlen = 10**(1.0*(mt-mc) + 1.5*(targmag-mt-dm) - dms)
		#
		winlen = 10**(1.0*(mt-mc) + 1.5*(targmag-mt) - dms - thisdm)		# as opposed to dm reletive to b=1.5.
	#
	#winlen=int(10*round(winlen/10))
	#print "winlen0: %d" % winlen
	if doint: winlen=int(round(winlen,-1))
	#print "winlen0: %d" % winlen
	if winlen<1: winlen=1
	#
	return winlen
			
def getNsample_bkp(m, mc, mt=7.6, dm=1.0, b1=1.0, b2=1.5, dms0=1.0, doint=False):
	# doint: convert to base10 integer (technically base 10**dm...
	Nom = getNomori(m=m, mc=mc, mt=mt, dm=dm, b1=b1, b2=b2)
	#
	if m<mt:
		dms=dms0
	if m>=mt:
		dms = (b1*(mt-mc) + b2*(m-mt))/(m-mc)
	#
	Nsamp = Nom*10**(-dms)
	#if doint: Nsamp=10*int(Nsamp/10)
	#if doint: Nsamp=10*int(round(Nsamp/10.))
	if doint: Nsamp = int(round(Nsamp,-1))
	#
	return Nsamp
#
def getNave(m=None, mc=None, mt=7.6, dm=1.0, b1=1.0, b2=1.5, dms0=1.0, dmave=1.0):
	if m==None or mc==None: return None
	#
	nsample = getNsample(m=m, mc=mc, mt=mt, dm=dm, b1=b1, b2=b2, dms0=dms0)
	nave=int(round(nsample*10**(-dmave)))
	if nave<1: nave=1
	#
	return nave
#
def mofN(N, mc, mt=7.6, b1=1.0, b2=1.5):
	logNt = b1*(mt-mc)
	if math.log10(N)<=logNt:
		# small earthquakes
		m = math.log10(N)/b1 + mc
		#
	if math.log10(N)>logNt:
		# large earthquakes
		#m = (math.log10(N) - b1*(mt-mc))/b2 + mt
		m = (math.log10(N) + mc*b1)/b2 + mt * (1.0 - b1/b2)
		#
	#
	return m

def plotrbimap(datafile='data/sumatrarbimap.dat', fignum=1, mc=4.75, catfile='cats/indonesia2012.cat'):
	# mc is prammed for sumatra
	Xi=[]		# ftime
	Yj=[]		# winlen
	Zij=[]	# 2d
	thisZ = []
	#
	Xs = []
	Nvals = []
	# Ys = []
	#thisN=None
	setlen=0
	firstN=None
	prevN=None
	#
	f=open(datafile)
	rowindex=0
	for rw in f:
		if rw[0]=='#': continue
		rws=rw.split('\t')
		#print rws
		N=int(rws[0])
		fdt=float(rws[3])
		r=float(rws[4])	# maybe - rws[4]
		#
		if firstN==None: firstN=N
		if N==firstN: 
			setlen+=1
			Xs+=[fdt]	# this will be the longest set and so will have all the dates
		#
		if N!=prevN:
			#print "switching: ", prevN, N
			if len(Zij)>0:
				Nvals+=[prevN]
				# early-fill short sequences
				while len(Zij[-1])<setlen:
					Zij[-1].insert(0, 0.0)
					#
			#	#
			Zij+=[[]]
			prevN=N
			thislen=0
		
		thislen+=1
		#	
		#Xi+=[fdt]
		#Yj+=[N]
		Zij[-1]+=[r]
	while len(Zij[-1])<setlen:
		Zij[-1].insert(0, 1.0)
		
	Nvals+=[N]
	#
	# now, fix up the x,y vectors so they can be used for contouring:
	ftimes = []
	nvect = []
	'''
	while len(ftimes)<len(Zij[0]):
		ftimes+=[Xs]
	#
	for n in Nvals:
		nvect+=scipy.ones(len(Zij[0])).tolist()
	'''
	#X,Y=numpy.meshgrid(Xs, Nvals)
	#
	f.close()
	#
	# now, massage and contour-plot:
	#
	# first, convert fdates to year.fraction:
	yrs=[]
	for fdt in Xs:
		thisdt = mpd.num2date(fdt)
		thisyr = float(thisdt.year)
		ydays = thisdt.timetuple()[7]
		yrs += [thisyr + float(ydays)/365.24]
	#
	# and convert N to magnitude:
	magax=[]
	deltamstars = 2.0	# the delta m factors (dm*, dm_s)
	for N in Nvals:
		magax+=[mofN(N, mc, mt=7.6, b1=1.0, b2=1.5) + deltamstars]
	#
	plt.figure(fignum)
	plt.clf()
	plt.ion()
	#plt.contourf(Xs, Nvals, Zij, 25, alpha=.75)
	plt.contourf(yrs, magax, Zij, 25, alpha=.75)
	plt.colorbar()
	#
	# and plot earthquakes:
	eqyrs = []
	eqmag = []
	eqNs = []	# we might plot in N-space.
	f=open(catfile, 'r')
	for rw in f:
		if rw[0]=='#': continue
		rws=rw.split('\t')
		thismag = float(rws[3])
		if thismag<mc+3.0: continue
		#
		#print rws[0], rws
		thisfdt =  mpd.datestr2num(rws[0])
		thisdt = mpd.num2date(thisfdt)
		thisyr = float(thisdt.year)
		ydays = thisdt.timetuple()[7]
		#
		eqyrs += [thisyr + float(ydays)/365.24]
		eqmag += [thismag]
	#
	f.close()
	#
	plt.plot(eqyrs, eqmag, "*", ms=15)
	plt.xlabel('Date')
	plt.ylabel('Magnitude (potential)')
		
	#
	#
	return [Zij, Xs, Nvals]	
	#
	# raw data will need to be filled so that each sequence is the same length.
	# just make one long list to start with. we'll want ftime, winLen, +/-r

def makeSaltonrbimapdata(catname='cats/salton.cat'):
	c1=eqp.eqcatalog([])
	c1.loadCatFromFile(catname)
	c1.rb=eqp.rbi.intervalRecordBreaker(None)
	#
	a=makerbiMapData(mags=[5.5, 8.2], mc=2.5, dN=5, datafile='data/saltonrbidata.dat', refresh=False, ct=c1)
	return a
	
def makerbiMapData(mags=[4.75, 9.1], mc=4.75, dN=10, datafile='data/sumatrarbimap.dat', refresh=False, ct=None):
	#
	if ct==None: ct=indoQuad(mc=mc, targmag=8.0)
	#mc=ct.mc
	N1 = Nofm(m=mags[0], mc=mc)/100
	N2 = Nofm(m=mags[1], mc=mc)/100
	# for now, 10-integers:
	N1=10*int(N1/10)
	N2=10*int(N2/10)
	N=N1
	if N<10: N=10
	print "N, N1, N2", N, N1, N2
	f=open(datafile, 'w')
	f.write('#rbi map\n#n, dt, ftd, log(r)\n')
	f.close()
	while N<=N2:
		print "N=%d" % N
		theserbi=ct.rb.getIntervalRatios(minmag=mc, windowLen=N, cat0=ct.getcat(0), deltaipos=1, logZ=None)
		f=open(datafile, 'a')
		#
		for rw in theserbi:
			f.write('%d\t%d\t%s\t%f\t%f\n' % (N, rw[0], rw[1], mpd.date2num(rw[1]), -math.log10(rw[2]) ) )
			#
		#
		f.close()
		N+=dN
	return ct
		
def indoHazMapMovie(catname='cats/indonesia2012.cat', catnum=0, mc=4.75, gridsize=None, ndithers=4, winlen=20, sigma=1.68, avlen=None, bigmag=7.9, fignum=0, dts=[dtm.datetime(2000,1,1, 0, 0, 0, 0, pytz.timezone('UTC')), dtm.datetime(2012,6,14, 0, 0, 0, 0, pytz.timezone('UTC'))], epicen=None, thresh=None, fileindex=0):
	#
	plt.ion()
	thisdt=dtm.timedelta(days=1)
	#
	# refresh cat just in case.
	objHM=indoHazMap(catname=catname, catnum=catnum, mc=mc, gridsize=gridsize, ndithers=ndithers, winlen=winlen, sigma=sigma, avlen=avlen, bigmag=bigmag, fignum=fignum, todt=dts[0], epicen=epicen, thresh=thresh, refreshcat=True)
	#
	# set common scale.
	#maxr = math.log10(objHM.winlen)/objHM.logZ
	#minr = -math.log10(objHM.winlen)/objHM.logZ		# typically, 1.0, -1.0.
																	# typically, we see r-> -1 after big earthquakes, but r->1.0
																	# is farily rare, so let's typically set it back a bit.
	#maxr = .8*maxr
	maxr = .75
	minr = -.9
	print "logZ: ", objHM.logZ, "winlen: %d/%f" % (objHM.winlen, math.log10(objHM.winlen))
	#
	thisdate=dts[0]
	#fileindex=0
	while thisdate<=dts[1]:
		objHM.hazMapTo(thisdate)
		# and force two pixels to maxr, minr
		objHM.hazmaps[0].rbRatios[0]=maxr
		objHM.hazmaps[1].rbRatios[0]=minr
		#
		plt.figure(fignum)
		objHM.simpleBoxes(fignum=fignum, thresh=thresh)
		#	
		z2=objHM.contourMap(fignum+2)
		#
		indx=[0,2,3]
		moviefolders = {indx[0]:'movies/sumatra/boxy/', indx[1]:'movies/sumatra/contscat/', indx[2]:'movies/sumatra/conts/'}
		fnameroots = {indx[0]:'sumatrahmboxy', indx[1]:'sumatrahmcc', indx[2]:'sumatrahmconts'}	# we'll be skipping the middle one anyway.
		#
		for i in indx:
			plt.figure(i)
			plt.legend(loc='best', numpoints=1)
			fcdate='%d-%d-%d' % (thisdate.year, thisdate.month, thisdate.day)
			plt.title('%s, $d\\lambda = %.4f$, $nDith=%d$, \n$m_c=%.2f$, $rblen=%d$, $avlen=%d$\n\n\n' % (fcdate, objHM.gridsize, objHM.ndithers, objHM.mc, objHM.winlen, objHM.avlen))
			#
			if i==indx[1]: continue # we just want boxy and naked conts.
			fileindexstr = '000000000000' + str(fileindex)
			fnameout = moviefolders[i] + fnameroots[i] + fileindexstr[-6:]
			plt.savefig(fnameout + '.png')
			#os.system('convert %s.png %s.jpg' % (fnameout, fnameout))
			#
		fileindex+=1

		thisdate+=thisdt
	#
	# and print the movie with something like:
	# avconv -i sumatrahmboxy%06d.png -r 50 ../sumatraboxy.mp4
	# aka: avconv -i contour2/salton-%05d.png.png -r 50 saltonconts1.mp4

# define nepal dataset; then call with makeHazMapDith(**nepal_prams)
nepal_epi_lon = 84.698
nepal_epi_lat = 28.175
nepal_dlon = 5.
nepal_dlat = 5.
nepal_prams = {'todt':dtm.datetime(2015, 4, 25, tzinfo=pytz.timezone('UTC')), 'gridsize':.1, 'mc':5.0, 'catnum':0, 'catname':'cats/nepal_2015.cat', 'epicen':None, 'thresh':None, 'ndithers':10, 'winlen':20, 'sigma':1.68, 'avlen':None, 'bigmag':6.5, 'fignum':0, 'lons':[nepal_epi_lon-nepal_dlon, nepal_epi_lon+nepal_dlon], 'lats':[nepal_epi_lat-nepal_dlat, nepal_epi_lat+nepal_dlat], 'dt0':None, 'logZ':None, 'nContours':15, 'refreshcat':True}



def generalHazmap(catname='cats/indonesia2012.cat', catnum=0, mc=4.75, gridsize=None, ndithers=4, winlen=20, sigma=1.68, avlen=None, bigmag=7.9, fignum=0, todt=dtm.datetime.now(pytz.timezone('UTC')), epicen=None, thresh=None, refreshcat=False, lats=[-9.0, 10.0], lons=[92.0, 106.0], dt0=None, logZ=None, nContours=14):
	#
	# todates: m9: 2004 12 24
	#			  m86: 2005 3 28
	#			  m85, 2007 9 12
	#          m86,m82: 2012 04 11
	if dt0==None: dt0=eqp.dtm.datetime(1990,1,1, tzinfo=pytz.timezone('UTC'))
	for D in [dt0, todt]:
		if D.tzinfo==None:
			D=dtm.datetime(*D.timetuple()[:-2], tzinfo=pytz.timezone('UTC'))
		#
	if catname==None:
		refreshcat=True
		catname='cats/indonesia2012.cat'
	#(2012, 4, 5)
	#
	# one mc, gridsize combination:mc=2.5, gridsize=.125, winlen=20 (approximately)
	# choose an optimal combination:
	# get gridsize/winlen:
	# getgridsize(mc=4.5, winlen=16, lat=45., sigmaExp=1.68)
	# def getWinlen(mc=4.5, gridsize=1.0, lat=45., sigmaExp=1.68):
	#
	plt.ion()
	#
	if refreshcat:
		#cl1=atp.catfromANSS(lon=[85.0, 105.0],lat=[-10.0, 10.0], minMag=3.5, dates0=[eqp.dtm.date(1990,1,1), eqp.dtm.date(2012,12,14)], fout=catname)
		cl1=atp.catfromANSS(lon=lons,lat=lats, minMag=mc, dates0=[dt0, todt], fout=catname)
	else:
		thiscat=eqp.eqcatalog()
		thiscat.loadCatFromFile(fname=catname, minmag=mc)
		cl1=thiscat.cat
	#
	if avlen==None:
		avlen=int(winlen/10)
		if avlen==0: avlen=1
		print "avelen: %d" % avlen
	# get gridsize/winlen:
	# getgridsize(mc=4.5, winlen=16, lat=45., sigmaExp=1.68)
	# def getWinlen(mc=4.5, gridsize=1.0, lat=45., sigmaExp=1.68):
	if gridsize==None and winlen==None:
		# can't do anything with that, can we?
		return None
	midlat=getMidLat(catname, mc)
	if gridsize==None:
		gridsize=getgridsize(mc, winlen, midlat, sigma)
	if winlen==None:
		winlen=getWinlen(mc, gridsize, midlat, sigma)
	#
	#objHM=makeHazMapDith(catname, catnum, mc, gridsize, ndithers, winlen, sigma, avlen, bigmag, fignum)
	objHM=makeHazMapDith(catname, catnum, mc, gridsize, ndithers, winlen, sigma, avlen, bigmag, fignum, logZ=logZ)
	objHM.mapres='i'
	#print "hazmap made in mscm"
	objHM.epicen=epicen
	objHM.hazMapTo(todt)
	fcindex=-1
	while cl1[fcindex][0] and (-fcindex)>len(cl1) > todt: fcindex-=1
	#fcdatetime = cl1[-1][0]
	fcdatetime = cl1[fcindex][0]
	#
	#fcdate='%d-%d-%d' % (todt.year, todt.month, todt.day)
	fcdate='%d-%d-%d' % (fcdatetime.year, fcdatetime.month, fcdatetime.day)
	print "hazMapTo() complete. ", todt, "  :  ",  fcdatetime
	#
	magres=objHM.guessMag(nsquares=1, mc=mc, deltas=2.0, beta=1.0, rblen=winlen, dithers=ndithers**2)
	magresDith=objHM.guessMag(nsquares=ndithers**2, mc=mc, deltas=2.0, beta=1.0, rblen=winlen, dithers=ndithers**2)
	print "mag-res: %f/%f (dithering=%d)" % (magres, magresDith, ndithers)
	#
	plt.figure(fignum)
	objHM.simpleBoxes(fignum=fignum, thresh=thresh)
	#
	# el mayor (baja, mexicali) earthquake:
	#hmx, hmy=objHM.hazmaps[0].catalog.catmap(-115.3, 32.13)
	#plt.plot([hmx], [hmy], 'b*', ms=18, zorder=10, alpha=.8)
	
	z2=objHM.contourMap(fignum=fignum+2, nconts=nContours)
	plt.figure(fignum+2)
	#plt.plot([hmx], [hmy], 'b*', ms=18)
	#
	bigevs = []
	for ev in objHM.hazmaps[0].catalog.cat:
		if ev[3]>bigmag:
			x,y=objHM.cm(ev[2], ev[1])
			for i in [0,2,3]:
				plt.figure(i)
				plt.plot([x], [y], '*', ms=2*ev[3], zorder=10, alpha=.6, label='m%.2f, %d-%d-%d' % (ev[3], ev[0].year, ev[0].month, ev[0].day))
	for i in [0,2,3]:
		plt.figure(i)
		plt.legend(loc='best', numpoints=1)
		plt.title('%s, $d\\lambda = %.4f$, $nDith=%d$, \n$m_c=%.2f$, $rblen=%d$, $avlen=%d$\n\n\n' % (fcdate, objHM.gridsize, objHM.ndithers, objHM.mc, objHM.winlen, objHM.avlen))
	#
	# note also:
	# objHM.hazmaps[0].catalog.catmap.drawmapscale(lon=-115., lat=33., lon0=-115.5, lat0=32.5, length=75., labelstyle='simple',
	#
	return objHM	
'''
def getNsample(targmag, mc, mt=7.6, dm=1.0, dms=1.0, doreduce=False):
	# mt=7.55 is probaby ok too.
	# dm=Bath const., dms=sampling exponent.
	# estimate winlen from a target magnitude.
	if targmag<mt:
		# "small" earthquake
		winlen=10**(targmag-dm-dms-mc)	# where 2.0 is dmstar + dmprime
	if targmag>=mt:
		dms2 = dms*(1.0*(mt-mc) + 1.5*(targmag-mt))/(targmag-mc)
		#winlen = 10**(1.0*(mt-mc) + 1.5*(targmag-mt-1.0) - dms)
		winlen = 10**(1.0*(mt-mc) + 1.5*(targmag-mt) - 2.0*dms2)
	#
	if doreduce==True:
		winlen=int(round(winlen,-1))
		#print "winlen0: %d" % winlen
		if winlen<1: winlen=1
	#
	return winlen
'''

def precursorAnalyzer(objCat, mc=None, dt0=None, winlen=None, targmag=None, fignum=0):
	# cut out the precursory signature from the RBTS and do things to it.
	if dt0==None: dt0=objCat.getMainEvent()[0]
	if mc==None: mc=objCat.mc
	if winlen==None:
		if targmag==None: targmag = objCat.getMainEvent()[3]
		winlen=getNsample(targmag, mc, doreduce=True)
	meanLen=int(winlen/10)
	#
	rbdata=objCat.getNRBratios(winlen=winlen)
	rbdataT=zip(*rbdata)
	Rlogs=map(math.log10, rbdataT[4])
	X  = rbdataT[0]
	T  = rbdataT[1]
	Tf = map(mpd.date2num, rbdataT[1])
	#
	print "targmag, winlen, meanlen: ", targmag, winlen, meanLen, len(Rlogs)
	#
	RlogsMean =  []
	RlogsMeanG = []
	trendlens  = [0]
	trendMaxs = []
	trendMins = []
	trendlen=0
	for i in xrange(meanLen, len(Rlogs)):
		RlogsMean  += [scipy.mean(Rlogs[(i-meanLen):i])]
		#
		if len(RlogsMean)>=2:
			if RlogsMean[-1]*RlogsMean[-2]<0.0:
				# last two values have opposite parity; new "burst".
				trendlen=0
				if RlogsMean[-2]>=0.: trendMaxs+=[trendlens[-2]]
				if RlogsMean[-2]<0.:  trendMins+=[trendlens[-2]]
			#if RlogsMean[-1]>=0.: trendSign=1
			#if RlogsMean[-1]<0.:  trendSign=-1
			trendSign=1
			trendlen+=(1*trendSign)
			trendlens+=[trendlen]
			#
			#if RlogsMean[-1]>=0: trend
	plt.figure(72)
	plt.ion()
	plt.clf()
	ax=plt.gca()
	ax.set_xscale('log')
	ax.set_yscale('log')
	trendMaxs.sort()
	trendMaxs.reverse()
	trendMins.sort()
	trendMins.reverse()
	print "max-min lens: ", len(trendMaxs), len(trendMins)
	plt.plot(trendMaxs, range(1, len(trendMaxs)+1), 'b.-', label='Omori')	
	plt.plot(trendMins, range(1, len(trendMins)+1), 'g.-', label='PSA')
	plt.legend(loc=0)
	#
	# let's gather some statistics about r-persistence:
	ngt, nlt, meangt, meanlt, Tgt, Tlt = 0,0,0.0,0.0,0.,0.		# number gt/lt, mean gt/lt, time gt/lt (time part might require
	rgtdt, rltdt=0., 0.														# more sophisticated analysis to catch "crossings".
	for i in xrange(1,len(RlogsMean)):		# for now, let's use the smoothed values and just skip the first value so we can do dt easily.
		r=RlogsMean[i]
		dt=Tf[i]-Tf[i-1]
		#
		if r>0.0:
			ngt+=1
			meangt+=r
			Tgt+=dt
			rgtdt=r*dt
		if r<0.0:
			nlt+=1
			meanlt+=r
			Tlt+=dt
			rltdt=r*dt
	print "raw stats (ngt, nlt, sumgt, sumlt, Tgt, Tlt, sum(rgtdt), sum(rltdt)): ", ngt, nlt, meangt, meanlt, Tgt, Tlt, rgtdt, rltdt
	print "normalized stats (meangt, meanlt, TmeanGT, TmeanLT): ", meangt/float(ngt), meanlt/float(nlt), rgtdt/(Tf[-1]-Tf[0]), rltdt/(Tf[-1]-Tf[0])
	#
	# get r>1 and r<1 clusters.
	#
	acor1=[]
	acorlen=1
	for i in xrange(acorlen, len(RlogsMean)):
		#acro1+=[(Rlogs[i]*Rlogs[i-acorlen])]
		acor1+=[scipy.prod(Rlogs[i-acorlen:i])]
	#
	# "signal 1:" (experimental acor*r(t))
	sig1=scipy.array(acor1)*scipy.array(RlogsMean[acorlen:])
	#
	print "lens: ", len(RlogsMean)
	rbursts=getrbursts(RlogsMean)
	#
	plt.figure(fignum)
	plt.ion()
	plt.clf()
	plt.plot(RlogsMean, '.-')
	plt.plot([0., len(RlogsMean)], [0., 0.], 'k-')
	'''
	plt.plot([Tf[0], Tf[-1]], [0., 0.], '--k')
	#plt.plot(Tf, Rlogs, 'r--')
	plt.plot(Tf[meanLen:], RlogsMean, 'b.-')
	#plt.plot(Tf[meanLen:], RlogsMeanG, 'g.-')
	#print len(Tf[meanLen+1:]), len(acro1)
	plt.plot(Tf[meanLen+acorlen:], acor1, 'm.-')
	plt.plot(Tf[meanLen+acorlen:], sig1, 'c.-')
	#
	#plt.plot(Tf[meanLen+acorlen:], RlogsMean, 'm.-')
	#plt.plot(Tf[meanLen+acorlen:], RlogsMean, 'c.-')
	#
	'''
	plt.figure(fignum+1)
	plt.clf()
	#
	plt.plot(trendlens, 'm--')
	for i in xrange(len(RlogsMean)):
		if RlogsMean[i]>=0.0: plt.plot([i], [trendlens[i]], 'b.')
		if RlogsMean[i]<0.0:  plt.plot([i], [trendlens[i]], 'r.')
	#
	plt.figure(fignum+2)
	plt.clf()
	#Y=scipy.array(RlogsMean)*scipy.array(map(math.sqrt, map(float, trendlens)))/math.sqrt(winlen)
	Y=scipy.array(RlogsMean)*scipy.array(map(float, trendlens))/float(winlen)
	plt.plot(Y, 'r.-')
	plt.plot([0., len(Y)], [0., 0.], 'k-')
	#
	plt.figure(fignum+3)
	plt.clf()
	plt.plot(Tf[meanLen:], Y, 'r-')
	plt.plot(Tf[meanLen:], RlogsMean, 'b-')
	plt.plot([Tf[0], Tf[-1]], [0., 0.], 'k-')
	#
	return [Tf, RlogsMean, trendlens]

def getrbursts(rs):
	# for now, just get clusters. we'll add time, etc. later.
	rsgt=[[]]
	rslt=[[]]
	parity = (rs[0]>=0)
	for r in rs:
		newparity = (r>=0)
		if newparity!=parity:
			if parity==1: rsgt+=[[]]
			if parity==0: rslt+=[[]]
		#
		parity=newparity
		if parity==1: rsgt[-1]+=[r]
		if parity==0: rslt[-1]+=[r]
	return [rsgt, rslt]
	

def nz2013hm(catname='cats/nz2013hm.cat', catnum=0, mc=4.0, gridsize=None, ndithers=10, winlen=20, sigma=1.68, avlen=None, bigmag=5.0, fignum=0, dt0=None, todt=(dtm.datetime.now(pytz.timezone('UTC'))), epicen=None, thresh=None, refreshcat=True, lons=[173.94, 174.94], lats=[41.21, 42.24]):
	if dt0==None: dt0=todt - dtm.timedelta(days=int(365.24*12))
	return generalHazmap(catname=catname, catnum=catnum, mc=mc, gridsize=gridsize, ndithers=ndithers, winlen=winlen, sigma=sigma, avlen=avlen, bigmag=bigmag, fignum=fignum, dt0=dt0, todt=todt, epicen=epicen, thresh=thresh, refreshcat=refreshcat, lats=lats, lons=lons)

def tohokuFSseries(datenum=0):
	llms=[142.373, 38.297]
	llfs=[142.842, 38.435]
	
	fcdates=[dtm.datetime(2011,3,9, tzinfo=pytz.timezone('UTC')), dtm.datetime(2011,3,10, tzinfo=pytz.timezone('UTC')), dtm.datetime(2011,3,11, tzinfo=pytz.timezone('UTC')), dtm.datetime(2011,3,15, tzinfo=pytz.timezone('UTC')), dtm.datetime(2011,4,15, tzinfo=pytz.timezone('UTC')), dtm.datetime.now(pytz.timezone('UTC'))]
	
	fignum=0
	a=japanhm(mc=5.0, lons=[135., 147.], todt=fcdates[datenum], fignum=fignum)
	c1=a.getcat()
	x0, y0 = c1.catmap(llms[0], llms[1])
	x1, y1 = c1.catmap(llfs[0], llfs[1])
	plt.figure(3)
	c1.catmap.plot([x0], [y0], 'r*', ms=15, alpha=.8, zorder=11)
	c1.catmap.plot([x0], [y0], 'k*', ms=17, alpha=.8, zorder=10)
	
	c1.catmap.plot([x1], [y1], 'c*', ms=13, alpha=.8, zorder=10)

	'''
	fignum=4
	a=japanhm(mc=5.0, lons=[135., 147.], todt=dtm.datetime(2011,3,10, tzinfo=pytz.timezone('UTC')), fignum=fignum)
	c1=a.getcat()
	plt.figure(fignum+3)
	x0, y0=c1.catmap(llms[0], llms[1])
	x1, y1=c1.catmap(llfs[0], llfs[0])
	c1.catmap.plot([x0], [y0], 'r*', ms=15, alpha=.8, zorder=11)
	c1.catmap.plot([x0], [y0], 'k*', ms=17, alpha=.8, zorder=10)
	c1.catmap.plot([x1], [y1], 'c*', ms=12, alpha=.8, zorder=11)

	fignum=8
	a=japanhm(mc=5.0, lons=[135., 147.], todt=dtm.datetime(2011,3,11, tzinfo=pytz.timezone('UTC')), fignum=fignum)
	c1=a.getcat()
	x0, y0=c1.catmap(llms[0], llms[1])
	x1, y1=c1.catmap(llfs[0], llfs[0])
	plt.figure(fignum+3)
	c1.catmap.plot([x0], [y0], 'r*', ms=15, alpha=.8, zorder=11)
	c1.catmap.plot([x0], [y0], 'k*', ms=17, alpha=.8, zorder=10)
	c1.catmap.plot([x1], [y1], 'c*', ms=12, alpha=.8, zorder=11)

	fignum=12
	a=japanhm(mc=5.0, lons=[135., 147.], todt=dtm.datetime(2011,3,15, tzinfo=pytz.timezone('UTC')), fignum=fignum)
	c1=a.getcat()
	x0, y0=c1.catmap(llms[0], llms[1])
	x1, y1=c1.catmap(llfs[0], llfs[0])
	plt.figure(fignum+3)
	c1.catmap.plot([x0], [y0], 'r*', ms=15, alpha=.8, zorder=11)
	c1.catmap.plot([x0], [y0], 'k*', ms=17, alpha=.8, zorder=10)
	c1.catmap.plot([x1], [y1], 'c*', ms=12, alpha=.8, zorder=11)
	'''


def japanhm(catname='cats/japanhm.cat', catnum=0, mc=4.75, gridsize=None, ndithers=5, winlen=20, sigma=1.68, avlen=None, bigmag=7.3, fignum=0, dt0=None, todt=(dtm.datetime.now(pytz.timezone('UTC'))), epicen=None, thresh=None, refreshcat=True, lons=[135., 146.], lats=[30., 41.5], nContours=None):
	if dt0==None: dt0=todt - dtm.timedelta(days=int(365.24*12))
	if nContours==None:
		nContours=[-.9, -.8, -.7, -.6, -.5, -.4, -.3, -.2, -.1, 0.0, .1, .2, .3, .4, .5, .6, .7]
	return generalHazmap(catname=catname, catnum=catnum, mc=mc, gridsize=gridsize, ndithers=ndithers, winlen=winlen, sigma=sigma, avlen=avlen, bigmag=bigmag, fignum=fignum, dt0=dt0, todt=todt, epicen=epicen, thresh=thresh, refreshcat=refreshcat, lats=lats, lons=lons, nContours=nContours)

# sumatra
def indoHazMap(catname='cats/indonesia2012.cat', catnum=0, mc=4.75, gridsize=None, ndithers=4, winlen=20, sigma=1.68, avlen=None, bigmag=7.9, fignum=0, todt=(dtm.datetime.now(pytz.timezone('UTC'))), epicen=None, thresh=None, refreshcat=False, lats=[-9.0, 10.0], lons=[92.0, 106.0]):
	ghm = generalHazmap(catname=catname, catnum=catnum, mc=mc, gridsize=gridsize, ndithers=ndithers, winlen=winlen, sigma=sigma, avlen=avlen, bigmag=bigmag, fignum=fignum, todt=todt, epicen=epicen, thresh=thresh, refreshcat=refreshcat, lats=lats, lons=lons)
	#
	plt.figure(3)
	c1=ghm.hazmaps[0].catalog
	x0,y0 = c1.catmap(95.854, 3.316)	# mainshock
	m1a=plt.plot([x0],[y0], 'b*', ms=18, zorder=11, alpha=.7)
	m1b=plt.plot([x0],[y0], 'b*', ms=21, zorder=10, alpha=.9)
	#
	x1,y1=c1.catmap(96.085, 2.824)	# foreshock? m=7.4, ~2002
	m2=plt.plot([x1],[y1], 'g*', ms=18, zorder=10, alpha=.8)
	c1.catmap.drawmapscale(lon=94.0, lat=3.0, lon0=95.25, lat0=3.0, length=150, labelstyle='simple', barstyle='fancy')
	#
	return ghm
#
def parkfieldhm(catname='cats/parkfieldhm.cat', catnum=0, mc=1.5, gridsize=None, ndithers=10, winlen=20, sigma=1.68, avlen=None, bigmag=5.0, fignum=0, dt0=None, todt=(dtm.datetime.now(pytz.timezone('UTC'))), epicen=None, thresh=None, refreshcat=True, lats=[35.415, 36.215], lons=[-120.774, -119.974], nContours=None):
	if dt0==None: dt0=todt - dtm.timedelta(days=int(365.24*12))
	if nContours==None:
		nContours=[-.9, -.8, -.7, -.6, -.5, -.4, -.3, -.2, -.1, 0.0, .1, .2, .3, .4, .5, .6, .7]
		
	z=generalHazmap(catname=catname, catnum=catnum, mc=mc, gridsize=gridsize, ndithers=ndithers, winlen=winlen, sigma=sigma, avlen=avlen, bigmag=bigmag, fignum=fignum, todt=todt, epicen=epicen, thresh=thresh, refreshcat=refreshcat, lats=lats, lons=lons, nContours=nContours)
	plt.figure(3)
	z.cm.drawstates(zorder=11, color='k')
	z.cm.drawrivers(zorder=11, color='b')
	return z
	
#
def chilehm2014(catname='cats/chile.cat', catnum=0, mc=5.0, gridsize=None, ndithers=10, winlen=20, sigma=1.68, avlen=None, bigmag=7.8, fignum=0, todt=dtm.datetime.now(pytz.timezone('UTC')), epicen=None, thresh=None, refreshcat=False, lats=[-30, -10.], lons=[-80.0, -65.0]):
	'''
	# mainshock: 2014-04-01 16:46:46 UTC-07:00, (lon, lat) = [-70.817, -19.642]
	#
	# note: wikipedia lists two epicenter locations:
	#x1,y1=c2.catmap(-72.733, -35.909)
	#x2,y2=c2.catmap(-73.239, -36.290)
	# there do not appear to be sufficient data to draw a precursory HM, but the RBTS come out beautifully.
	#
	'''
	print "getting chilehm.", lats, lons
	#
	z=generalHazmap(catname=catname, catnum=catnum, mc=mc, gridsize=gridsize, ndithers=ndithers, winlen=winlen, sigma=sigma, avlen=avlen, bigmag=bigmag, fignum=fignum, todt=todt, epicen=epicen, thresh=thresh, refreshcat=refreshcat, lats=lats, lons=lons)
	plt.figure(3)
	z.cm.drawstates(zorder=11, color='k')
	z.cm.drawrivers(zorder=11, color='b')
	
	ms2014=[-70.817, -19.642]
	x,y = z.cm(ms2014[0], ms2014[1])
	plt.figure(2)
	plt.plot([x], [y], 'r*', zorder=11, alpha=.8, ms=18)
	
	return z

def chilehm(catname='cats/chile.cat', catnum=0, mc=4.5, gridsize=None, ndithers=4, winlen=20, sigma=1.68, avlen=None, bigmag=5.9, fignum=0, todt=dtm.datetime.now(pytz.timezone('UTC')), epicen=None, thresh=None, refreshcat=False, lats=[-40, -10.], lons=[-80.0, -65.0]):
	#
	# note: wikipedia lists two epicenter locations:
	#x1,y1=c2.catmap(-72.733, -35.909)
	#x2,y2=c2.catmap(-73.239, -36.290)
	# there do not appear to be sufficient data to draw a precursory HM, but the RBTS come out beautifully.
	#
	print "getting chilehm.", lats, lons
	#
	z=generalHazmap(catname=catname, catnum=catnum, mc=mc, gridsize=gridsize, ndithers=ndithers, winlen=winlen, sigma=sigma, avlen=avlen, bigmag=bigmag, fignum=fignum, todt=todt, epicen=epicen, thresh=thresh, refreshcat=refreshcat, lats=lats, lons=lons)
	plt.figure(3)
	z.cm.drawstates(zorder=11, color='k')
	z.cm.drawrivers(zorder=11, color='b')
	
	
	return z


def nevadahm(catname='cats/nevada.cat', catnum=0, mc=2.5, gridsize=None, ndithers=4, winlen=20, sigma=1.68, avlen=None, bigmag=5.9, fignum=0, todt=(dtm.datetime.now(pytz.timezone('UTC'))), epicen=None, thresh=None, refreshcat=False, lons=[-120.5, -113.25], lats=[34.25, 42.25]):
	#
	z=generalHazmap(catname=catname, catnum=catnum, mc=mc, gridsize=gridsize, ndithers=ndithers, winlen=winlen, sigma=sigma, avlen=avlen, bigmag=bigmag, fignum=fignum, todt=todt, epicen=epicen, thresh=thresh, refreshcat=refreshcat, lats=lats, lons=lons)
	plt.figure(3)
	z.cm.drawstates(zorder=11, color='k')
	z.cm.drawrivers(zorder=11, color='b')
	return z


def elmayorhm(catname='cats/elmayor.cat', catnum=0, mc=2.5, gridsize=None, ndithers=10, winlen=20, sigma=1.68, avlen=None, bigmag=6.0, fignum=0, todt=(dtm.datetime.now(pytz.timezone('UTC'))), epicen=None, thresh=None, refreshcat=True, lats=[31.0, 35.0], lons=[-118.0,-114.0], dt0=dtm.datetime(1990,1,1, tzinfo=pytz.timezone('UTC')), nContours=None):
	# lats=[32.128-1.2, 32.128+1.2], lons=[-115.303-1.2, -115.303+1.2]
	# lats=[31.85,34.55], lons=[-118.85, -114.15]
	# 32.128, -115.303
	# L ~ 66 km -- > .6 deg.
	z=generalHazmap(catname=catname, catnum=catnum, mc=mc, gridsize=gridsize, ndithers=ndithers, winlen=winlen, sigma=sigma, avlen=avlen, bigmag=bigmag, fignum=fignum, todt=todt, epicen=epicen, thresh=thresh, refreshcat=refreshcat, lats=lats, lons=lons, nContours=nContours)
	plt.figure(3)
	z.cm.drawstates(zorder=11, color='k')
	z.cm.drawrivers(zorder=11, color='b')
	return z

def saltonhm(catname='cats/salton.cat', catnum=0, mc=2.25, gridsize=None, ndithers=10, winlen=20, sigma=1.68, avlen=None, bigmag=6.0, fignum=0, todt=(dtm.datetime.now(pytz.timezone('UTC'))), epicen=None, thresh=None, refreshcat=False, lats=[31.85,34.55], lons=[-118.85, -114.15], dt0=dtm.datetime(1990,1,1, tzinfo=pytz.timezone('UTC'))):
	#
	z=generalHazmap(catname=catname, catnum=catnum, mc=mc, gridsize=gridsize, ndithers=ndithers, winlen=winlen, sigma=sigma, avlen=avlen, bigmag=bigmag, fignum=fignum, todt=todt, epicen=epicen, thresh=thresh, refreshcat=refreshcat, lats=lats, lons=lons)
	plt.figure(3)
	z.cm.drawstates(zorder=11, color='k')
	z.cm.drawrivers(zorder=11, color='b')
	return z
#
def cascadiahm(catname='cats/cascadia.cat', catnum=0, mc=3.75, gridsize=None, ndithers=10, winlen=20, sigma=1.68, avlen=None, bigmag=6.8, fignum=0, epicen=None, thresh=None, refreshcat=True, lats=[46.0, 55.0], lons=[-133.0, -121.0], dt0=dtm.datetime(1990,1,1, tzinfo=pytz.timezone('UTC')), todt=(dtm.datetime.now(pytz.timezone('UTC')))):
	#
	# as per teragrid execution model:
	# f=saltonstuff(mc=2.5, ndithers=10, winlen=20, sigma=1.68, thresh=0., lons=[-127.5, -117.5], lats=[36.75, 42.0], namebase='norcal', equakes=None, contres=totalres)
	#
	z=generalHazmap(catname=catname, catnum=catnum, mc=mc, gridsize=gridsize, ndithers=ndithers, winlen=winlen, sigma=sigma, avlen=avlen, bigmag=bigmag, fignum=fignum, todt=todt, epicen=epicen, thresh=thresh, refreshcat=refreshcat, lats=lats, lons=lons)
	plt.figure(3)
	z.cm.drawstates(zorder=11, color='k')
	z.cm.drawrivers(zorder=11, color='b')
	return z

def cascadiamovie(catname='cats/cascadia.cat', catnum=0, mc=4.0, gridsize=None, ndithers=10, winlen=20, sigma=1.68, avlen=None, bigmag=6.8, fignum=0, dts=[dtm.datetime(2010, 1, 1, 0, 0, 0, 0, pytz.timezone('UTC')), dtm.datetime.now(pytz.timezone('UTC')), dtm.timedelta(days=1)], catdate=dtm.datetime(1990, 1, 1, 0, 0, 0, 0, pytz.timezone('UTC')), epicen=None, thresh=None, refreshcat=True, lats=[46.0, 54.0], lons=[-132.0, -121.0], logZ=None, nameroot='cascadia', moviedir='movies/cascadia/'):
	#
	todt=dts[1]
	dt0=catdate
	quakes = []
	#
	hm=generalHazMovie(catname=catname, catnum=catnum, mc=mc, gridsize=gridsize, ndithers=ndithers, winlen=winlen, sigma=sigma, avlen=avlen, bigmag=bigmag, fignum=fignum, dts=dts, catdate=catdate, epicen=epicen, thresh=thresh, refreshcat=refreshcat, lats=lats, lons=lons, logZ=logZ, nameroot=nameroot, moviedir=moviedir, events=quakes)
	#
	#avconvstr='avconv -i %scontour/%s-%05d.png -r 50 %s/%sconts.mp4' % (moviedir, nameroot, moviedir, nameroor)
	#os.system(avconvstr)
	#
	return hm
#
def socalhm(catname='cats/socal.cat', catnum=0, mc=2.25, gridsize=None, ndithers=10, winlen=20, sigma=1.68, avlen=None, bigmag=6.0, fignum=0, epicen=None, thresh=None, refreshcat=True, lats=[31., 37.], lons=[-120., -114.], dt0=dtm.datetime(1990,1,1, tzinfo=pytz.timezone('UTC')), todt=(dtm.datetime.now(pytz.timezone('UTC')))):
	#
	# as per teragrid execution model:
	# f=saltonstuff(mc=2.5, ndithers=10, winlen=20, sigma=1.68, thresh=0., lons=[-127.5, -117.5], lats=[36.75, 42.0], namebase='norcal', equakes=None, contres=totalres)
	#
	z=generalHazmap(catname=catname, catnum=catnum, mc=mc, gridsize=gridsize, ndithers=ndithers, winlen=winlen, sigma=sigma, avlen=avlen, bigmag=bigmag, fignum=fignum, todt=todt, epicen=epicen, thresh=thresh, refreshcat=refreshcat, lats=lats, lons=lons)
	plt.figure(3)
	z.cm.drawstates(zorder=11, color='k')
	z.cm.drawrivers(zorder=11, color='b')
	return z
#
def norcalhm(catname='cats/norcal.cat', catnum=0, mc=2.5, gridsize=None, ndithers=10, winlen=20, sigma=1.68, avlen=None, bigmag=6.0, fignum=0, epicen=None, thresh=None, refreshcat=True, lats=[36.75,42.0], lons=[-129.5, -117.5], dt0=dtm.datetime(1990,1,1, tzinfo=pytz.timezone('UTC')), todt=(dtm.datetime.now(pytz.timezone('UTC')))):
	#
	# as per teragrid execution model:
	# f=saltonstuff(mc=2.5, ndithers=10, winlen=20, sigma=1.68, thresh=0., lons=[-127.5, -117.5], lats=[36.75, 42.0], namebase='norcal', equakes=None, contres=totalres)
	#
	z=generalHazmap(catname=catname, catnum=catnum, mc=mc, gridsize=gridsize, ndithers=ndithers, winlen=winlen, sigma=sigma, avlen=avlen, bigmag=bigmag, fignum=fignum, todt=todt, epicen=epicen, thresh=thresh, refreshcat=refreshcat, lats=lats, lons=lons)
	plt.figure(3)
	z.cm.drawstates(zorder=11, color='k')
	z.cm.drawrivers(zorder=11, color='b')
	return z
#
def saltonstuff(dorefresh=False, todt=dtm.datetime.now(pytz.timezone('UTC')), fromdt=dtm.datetime(1990,1,1, tzinfo=pytz.timezone('UTC')) , lats=[31.85,34.55], lons=[-118.85, -114.15], ndithers=10, winlen=20, doquad=True, mc=2.25, alpha=.9):
	# if timezone info must be specified, use the
	# tzi = datetime.tzinfo
	# dtm.datetime.now(pytz.timezone('UTC')
	#
	plt.ion()
	reload(bcp)
	bc1=bcp.BASScast()
	#
	thiswinlen=winlen
	thisavlen=int(thiswinlen/10)
	if thisavlen==0: thisavlen=1
	#
	if dorefresh == True: z=makeSocalcat(lat=lats, minMag=mc, lon=lons, fout='cats/salton.cat', dates0=[fromdt, todt])
	thisfcdt=todt
	#thisfcdt=dtm.datetime(2011,6,1, 0, 0, 0, 0, pytz.timezone('UTC'))
	q=makeSocalMap(catname='cats/salton.cat', winlen=thiswinlen, avlen=thisavlen, todt=thisfcdt, logZ=None, thresh=None, ndithers=ndithers, alpha=alpha)
	ct=q.hazmaps[0].catalog
	if doquad: qq=ct.rbomoriQuadPlot(catnum=0, mc=mc, winlen=400, rbavelen=40, bigmag=6.0, fignum=6, plotevents=True, intlist=[25, 50, 100, 200, 300, 400, 500], logZ=None)
	#
	# now, make some KMZ/L of the RB forecast; use the BASScast object to write the KML
	logZ = math.log10(float(thiswinlen))		# log-normalization factor (so -1 < r < 1 )
																# note: we're plotting r in log space, so the quantity in question
	logZinv = 1.0/math.log10(float(thiswinlen))															# is log(n1/n2)/log(N)
	print "logZ: %f, %f" % (logZ, logZinv)
	#
	rbarrays = q.contourArrays()
	#
	# so log-normalize in linear space:
	'''
	for i in xrange(len(rbarrays[2])):
		for ii in xrange(len(rbarrays[2][i])):
			if rbarrays[2][i][ii]!=None: rbarrays[2][i][ii] = rbarrays[2][i][ii]**(logZinv)
	'''		
	#
	contset=bc1.getContourSet(X_i = rbarrays[0], Y_i = rbarrays[1], Z_ij = rbarrays[2], contres = 10)	# nConts is some factor of this contyres
	bc1.writeKMLfile(cset=contset, fout='tmp/rbhmsalton.kml', colorbarname='tmp/rbhmsaltonscale.png', equakes=None)
	
	# time series (see specific time-series below):
	#logZ = math.log10(float(400))
	#logZinv = 1.0/logZ
	
	#rats1=ct.rb.getIntervalRatios(minmag=2.5, windowLen=400, cat0=ct.getcat(0), deltaipos=1, logZ=None)	# normalized
	#
	'''
	plt.figure(12)
	plt.clf()
	# normalize then mean (rats1 is normalized; take mean in this function).
	#a1=ct.rb.plotIntervalRatios(minmag=2.5, windowLen=400, cat0=ct.getcat(0), hitThreshold=1.0, bigmag=6.5, fignum=12, ratios=rats1, deltaipos=1, avlen=40, logZ=None)
	qqnorm=ct.rbomoriQuadPlot(catnum=0, mc=2.5, winlen=400, rbavelen=40, bigmag=6.0, fignum=12, plotevents=True, intlist=[25, 50, 100, 200, 300, 400, 500], logZ=None)
	plt.title('mean-normalized RB sequence')
	'''
#	
	#return ct
	return q

def makeKMLfromHM(objHM, outdir='kml', fnameroot='hmkml', logZ=None):
	if logZ==None: logZinv=1.0/(math.log10(float(objHM.winlen)))
	kmlfile = '%s/rbhm%s.kml' % (outdir, fnameroot)
	colorbarname='%s/rbhm%sscale.png' % (outdir, fnameroot)
	#
	bc1=bcp.BASScast()	# just a container for the KML functions
	rbarrays = objHM.contourArrays()
	contset=bc1.getContourSet(X_i = rbarrays[0], Y_i = rbarrays[1], Z_ij = rbarrays[2], contres = 16)
	bc1.writeKMLfile(cset=contset, fout=kmlfile, colorbarname=colorbarname, equakes=None)
	#
	kmzfile='%s/rbhm%s.kmz' % (outdir, fnameroot)
	os.system('zip -D -j %s %s %s' % (kmzfile, kmlfile, colorbarname))
	#
	return None

def drawgarnish(cm, mapres='i'):
	# cm is a catmap
	cm.drawcountries(zorder=1, linewidth=1)
	cm.drawstates(zorder=1, linewidth=1)
	#cm.drawlsmask(land_color='0.8', ocean_color='c', lsmask=None, lsmask_lons=None, lsmask_lats=None, lakes=True, resolution=mapres, grid=5)
	cm.drawmapscale(lon=-115., lat=33., lon0=-115.5, lat0=32.5, length=50, labelstyle='simple', barstyle='fancy')
	cm.drawrivers(linewidth=1, color='b')

def makeJapanHazMap(catname='cats/japancat4b.cat', catnum=0, mc=4.5, gridsize=1.0, phi=[0,0], rblen=32, avlen=1, bigmag=6.5, fignum=0):
	# for now, the hazard map object is in haitistuff.py. this should be moved to its own module, yodapy, or something like that.
	# getJapanQuadPlot(catname='cats/japancat4.cat', rbint=128, intervals=[64, 128, 256, 512], catnum=1, avelen=1):
	#rbcat=getjapanRB(catname, rblen, False)	# returns [catalog, rbObj]
	plt.ion()
	catalog=eqp.eqcatalog()
	catalog.loadCatFromFile(catname)
	#catalog=rbcat[0]
	#catalog.rb=rbcat[1]
	#
	#objHM=hpp.rbHazMap(catalog, catnum, mc, gridsize, phi, rblen, avlen, bigmag, fignum)
	objHM=hmp.rbHazMapLean(catalog, catnum, mc, gridsize, phi, rblen, avlen, bigmag)
	#
	return objHM			
	
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

def getMidLat(thiscat, mc=None):
	objcat=eqp.eqcatalog()
	objcat.loadCatFromFile(thiscat, mc)
	lats=map(operator.itemgetter(1), objcat.cat)
	midlat=numpy.mean(lats)
	#
	return midlat


def findClusters(gridlist=[], L=None):
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
	if L==None:
		# before giving up, guess a square:
		squareL=((float(gridlen))**.5)
		if squareL%1==0:
			# it could be a square...
			L=squareL
		else:
			# we don't know how to parse it.
			print "array width not defined."
			return None
	#
	clustIndex=0
	#
	i=0
	print "L=%d" % L
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
		
	for rw in clustMap:
		srcindex=newclusts[rw[0]]
		targetindex=newclusts[rw[1]]
		#
		#clustList[rw[1]][1]+=clustList[rw[0]][1]
		#clustList[rw[0]][1]=[]
		clustList[targetindex][1]+=clustList[srcindex][1]
		clustList[srcindex][1]=[]
		newclusts[srcindex] = targetindex
	
	#return [clustList, clustMap, clusterIDs]
	# return clusters.
	rlist=[]
	for rw in clustList:
		if len(rw[1])==0: continue
		rlist+=[rw[1]]
	return rlist

# a=mhp.generalHazMovie(catname='cats/isabel.cat', mc=2.25, ndithers=10, bigmag=4.0, dts=[mhp.dtm.datetime(2011,1,1, tzinfo=mhp.pytz.timezone('UCT')), dtm.datetime.now(mhp.pytz.timezone('UCT')),  dtm.timedelta(days=1)], nameroot='isabelmovie', lons=[-118.5, -118.25], lats=[35.1, 35.5], refreshcat=True, moviedir='movies/isabel')

'''
def doit():
   .....:     reload(hmp)
   .....:     reload(hmp.hmp)
   .....:     A=hmp.makeJapanHazMapDith(ndithers=10)
   .....:     A.simpleBoxes()


'''
	
