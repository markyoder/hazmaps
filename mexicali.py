import math
import scipy
import pylab
import matplotlib.pyplot as plt
import scipy.optimize as spo
#
from matplotlib.patches import Ellipse
import matplotlib.dates as mpd
#
import string
import sys
import os
import random
import time
#
# gamma function lives here:
from matplotlib import axis as aa
from threading import Thread
#
#
import datetime as dtm
import calendar
import operator
#import urllib
#import MySQLdb

import yodapy as ypp
import rbIntervals as rbi
import haitistuff as hpp

#def getMexiQuad(catname='cats/mexicali.cat', rbint=398, intervals=[64, 128, 256, 512], catnum=0, avelen=None, thislw=2.0):
def getMexiQuad(catname='cats/mexicali.cat', rbint=500, intervals=[64, 128, 256, 512], catnum=0, avelen=None, thislw=2.0, mc=2.5, bigmag=6.0, MEV=None, fignum=0):
	# rbint=500 is nominal
	plt.ion()
	if avelen==None: avelen=int(rbint/10)	# self similar averaging interval
	if avelen<1: avelen=1
	print "avelen: %d" % avelen
	rbcat=getRB(catname, rbint, False)	# returns [catalog, rbObj]
	thiscat=rbcat[0]
	thiscat.rb=rbcat[1]	# this is the modified form of the catalog object that quadPlot expects.
	#thiscat.addLatLonSubcat('centnortjapan', thiscat.getcat(0), [33., 42.], [138., 145.])
	#mc=2.5
	#bigmag=6.0
	# get main event?
	if MEV==None:
		mycat=thiscat.getcat(0)
		i=len(mycat)-1
		while i>=0 and MEV==None:
			if mycat[i][3]>7.0 and mycat[i][0]>dtm.datetime(2010,4,4) and mycat[i][0]<dtm.datetime(2010,4,5):
				MEV=mycat[i]
				print "MEV: %s" % MEV
			i-=1
	#
	#return rbcat
	#
	polys=[]
	#polys+=[[[138., 33.], [145., 33.], [145., 42.], [138., 42.], [138., 33.]]]
	#polys+=[[[138., 41.5], [138, 34.5], [140., 30.], [145., 30.], [147, 41.5]]]
	#polys+=[[[164.7070, 58.3095], [151.4355, 51.5087], [143.0859, 45.5833], [140.1855, 45.6448],  [138.2520, 42.5531], [135.0879, 37.0201], [137.1973, 32.1012], [138.1641, 21.5348], [151.7871, 21.2894], [164.7070, 58.3095]]]
	#
	# polycat(self, cat=None, verts=None):
	#thiscat.subcats+=[['jpoly1', thiscat.polycat(thiscat.getcat(0), polys[1]) ]]
	#thiscat.subcats+=[['ploycat2', thiscat.polycat(thiscat.getcat(0), polys[2]) ]]
	#
	#rbomoriQuadPlot(catalog=None, catnum=0, mc=2.5, winlen=501, rbthresh=1.0, bigmag=5.0, fignum=0, intlist=[32, 256, 512], rbavelen=1, thislw=1.0)
	#A=hpp.rbomoriQuadPlot(catalog=thiscat, catnum=catnum, mc=mc, winlen=rbint, rbthresh=1.0, bigmag=bigmag, fignum=2, intlist=intervals, rbavelen=avelen, thislw=thislw)
	# catnum=0, mc=2.5, winlen=501, rbthresh=1.0, bigmag=5.0, fignum=0, intlist=None, rbavelen=1, thislw=1.0
	A=thiscat.rbomoriQuadPlot(catnum=catnum, mc=mc, winlen=rbint, rbthresh=1.0, bigmag=bigmag, fignum=fignum, intlist=intervals, rbavelen=avelen, thislw=thislw, mainEV=MEV, plotevents=True)
	plt.figure(fignum)
	myax=thiscat.axlist[3]
	x,y=thiscat.catmap(MEV[2], MEV[1])
	plt.plot([x], [y], 'k*', ms=24,)
	plt.plot([x], [y], 'r*', ms=20, label='mainshock')
	handles, labels=myax.get_legend_handles_labels()
	handles.pop(-2)
	labels.pop(-2)
	myax.legend(handles, labels, loc='lower left', numpoints=1)
	#plt.legend(loc='lower left', numpoints=1)

	
	#A=hpp.rbomoriQuadPlot(thiscat, catnum, mc, rbint, 1.0, 6.5, 2, intervals, avelen, thislw)
	# axes list: A[0].axlist	# we add the axes list to the catalog before it is returned.
	#
	# draw the polygon?
	#return polys
	
	if catnum>0 and len(polys)>catnum:
		# assume we're using a polygon, but this should be better organized in the future (maybe a dictionary of polygons).
		thispoly=polys[catnum-1]
		pX=map(operator.itemgetter(0), thispoly)
		pY=map(operator.itemgetter(1), thispoly)
		polyX, polyY = A.catmap(pX, pY)
		#
		A.axlist[3].fill(polyX, polyY, alpha=.1, color='b')
		A.axlist[3].plot(polyX, polyY, 'b-', lw=3, alpha=.75)
	#
	#
	#return [polyX, polyY]
	return rbcat

def getRB(catname='cats/mexicali.cat', nits=256, doplot=True, mc=2.5, bigmag=6.0, myMEV=None):
	#
	catlist=[]
	# datetimeFromString
	dindex=0
	fin=open(catname)
	for rw in fin:
		if rw[0]=='#' or rw[0]=='\n': continue
		rws=rw.split('\t')
		#print rws
		strDtm=rws[0]
		if ":" in rws[1]:
			# date and time are separated by a tab, not space
			dindex=1
			strDtm=strDtm + ' ' + rws[1]
			#print "fixing."
		
		
		if float(rws[3+dindex])<mc: continue
		lat=float(rws[1+dindex])
		lon=float(rws[2+dindex])
		mag=float(rws[3+dindex])
		thisdt=ypp.datetimeFromString(strDtm)
		#
		catlist+=[[thisdt, lat, lon, mag]]
		#if myMEV==None:
		#	# make El Mayor the MEV:
		#	if thisdt>dtm.datetime(2010,4,4) and catlist[-1][3]>7.0: myMEV=catlist[-1][3]
		#print "row."
	
	rbc=ypp.eqcatalog(catlist)
	#return rbc
	rbj=rbi.intervalRecordBreaker(None)
	rbj.fullcat=catlist
	rbj.shockcat=rbj.fullcat
	#
	if doplot:
		rbc.plotCatMap()
		#
		# minmag=3.0, windowLen=10, cat0=None, hitThreshold=1.0, bigmag=5.0, thisAx=None, ratios=None, deltaipos=1, avlen=1, mainEV=None
		rbj.plotIntervalRatios(minmag=mc, windowLen=nits, cat0=rbj.fullcat, hitThreshold=1.0, bigmag=bigmag, mainEV=myMEV)
	
	return [rbc, rbj]	# [catalogObject, rbObject]

def getMexicaliCat(outFile='cats/mexicali.cat', startDate=dtm.datetime(1990, 01, 01), endDate=dtm.datetime.now(), minmag=2.0, lats=[31.0, 35.0], lons=[-118.0, -113.5] ):
	# catfromANSS(lon=[135., 150.], lat=[30., 41.5], minMag=4.0, dates0=[dtm.date(2005,01,01), None], Nmax=999999, fout='cats/mycat.cat'):
	CL=ypp.catfromANSS(lon=lons, lat=lats, minMag=minmag, dates0=[startDate, endDate], fout=outFile)
	
	return CL
