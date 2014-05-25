from math import *	
from scipy import *
import scipy
from pylab import *
from matplotlib import *
import numpy.fft as nft
import scipy.optimize as spo
#from matplotlib import pyplot as plt
import pylab as plt
from matplotlib import rc
from matplotlib.patches import Ellipse

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
from threading import Thread
#
#
import datetime as dtm
import calendar
import operator
import urllib
import MySQLdb

import yodapy as yp
import rbIntervals as rbi

def updateDB(catID=523):
	y.updateANSS2SQL(catID)



############################################
def freshHaitiRatios(wLens=[10, 15, 20, 30], bigShock=5.0, haitiCat='hati.cat', doShow=False):
	# full set of Haiti rb ratios from a fresh catalog...
	# update catalog:
	print "update MySQL catalog..."
	y.updateANSS2SQL()
	print "update local catalog, %s" % haitiCat
	getHaitiCatFromSQL(haitiCat)
	#
	return haitiRatios(wLens, bigShock, haitiCat, doShow)
	
def haitiRatios(wLens=[10, 15, 20, 30], bigShock=5.2, haitiCat='hati.cat', doShow=False):
	# full set of Haiti rb ratios
	# plot GR distributions:
	rbHaiti=getStandardHaitiRB()
	rbHaiti.GRshock(False, fname='images/haitiShock-GRdist.png')
	rbHaiti.GRfullcat(False, fname='images/haitiFull-GRdist.png')
	rbHaiti=None
	
	# calc ratios for wLens:
	print "calc ratios for wLens"
	for wLen in wLens:
		print "rb ratios for wLen=%d" % wLen
		getHaitiRBratio(wLen, bigShock, haitiCat, doShow, dtm.datetime.now())
	#
	return 0
	
#########################################


#######################################
def getHaitiCatFromSQL(outFile='haiti.cat', startDate=dtm.datetime(2000, 01, 01), endDate=dtm.datetime.now(), minmag=2.5, lats=[16.5, 20.5], lons=[-74.0, -70.0], catID=523 ):
	# haiti catalog is in flux, so let's use the catalog. we could go directly from ANSS if we like, but i think this will be easier.
	# note that the default prams are a bit broad, but we seem to get a signal anyway.
	#
	sqlSelect = "select eventDateTime, lat, lon, mag from Earthquakes where catalogID=%d and eventDateTime between '%s' and '%s' and lat between %f and %f and lon between %f and %f and mag>=%f order by eventDateTime asc" % (catID, str(startDate), str(endDate), lats[0], lats[1], lons[0], lons[1], minmag)
	sqlHost = 'localhost'
	sqlUser = 'myoder'
	sqlPassword = 'yoda'
	sqlPort = 3306
	sqlDB = 'QuakeData'
	myConn = MySQLdb.connect(host=sqlHost, user=sqlUser, passwd=sqlPassword, port=sqlPort, db=sqlDB)
	r1=myConn.cursor()
	r1.execute(sqlSelect)
	#	
	fout=open(outFile, 'w')
	fout.write("# haiti catalog from hatistuff.py\n#prams: outfile, startDate, endDate, minmag, lats, lons, catID\n")
	fout.write("#%s\t%s\t%s\t%s\t%s\t%s\t%s\d\n" % (outFile, startDate, endDate, minmag, lats, lons, catID))
	# catalog format (tab delimited):
	# 2004/09/27	01:53:42.78	35.5277	-120.8433	1.42
	catlen=0
	for rw in r1:
		fout.write("%d/%d/%d\t%d:%d:%d.%d\t%f\t%f\t%f\n" % (rw[0].year, rw[0].month, rw[0].day, rw[0].hour, rw[0].minute, rw[0].second, rw[0].microsecond, rw[1], rw[2], rw[3]))
		catlen+=1
	fout.close()
	r1.close()
	myConn.close()
	
	return catlen
def getStandardHaitiRB(hlat=18.46, hlon=-72.75, haitiTheta=-15, ra=1.25, rb=.25, haitiCat='hati.cat'):
	#reload(rbi)
	rbHaiti=rbi.intervalRecordBreaker(haitiCat, haitiTheta, hlat, hlon, ra, rb)
	rbHaiti.setAftershockCatalog(haitiCat, haitiTheta, hlat, hlon, ra, rb, None, dtm.datetime.now(), 0)
	return rbHaiti

def getHaitiCatPlot(doShow=False, doSave=True, rbh=None, saveName='images/catalogPlot.png'):
	# rbh is a Haiti specific intervalRecordBreaker() object.
	# note: this functionality has been integrated into the intervalRecordBreaker() object as obh.xyPlotCatalogs(doShow, doSave, saveName, ellipseCenter=[x,y])
	if rbh==None: rbh=getStandardHaitiRB()
	fcat=[]
	scat=[]
	for rw in rbh.shockCat:
		scat+=[rw[0:4]]
	for rw in rbh.fullCat:
		if rw not in scat: fcat+=[rw]
	#return [scat, fcat]
	
	plt.figure(0)	
	plt.clf()
	#
	ax=plt.gca()
	#el = Ellipse((-72.533,18.457), 1.25, .25, 15, facecolor='r', alpha=0.5)
	el = Ellipse((-72.75,18.46), 2.0*1.25, 2.0*.25, 15.0, facecolor='b', alpha=0.4)
	ax.add_artist(el)
	#
	plt.plot(map(operator.itemgetter(2), rbh.fullCat), map(operator.itemgetter(1), rbh.fullCat), '+')
	plt.plot(map(operator.itemgetter(2), rbh.shockCat), map(operator.itemgetter(1), rbh.shockCat), '.')
	#plt.plot(map(operator.itemgetter(2), fcat), map(operator.itemgetter(1), fcat), '+', label='Full Catalog')
	#plt.plot(map(operator.itemgetter(2), scat), map(operator.itemgetter(1), scat), '.', label='Aftershock zone')
	plt.plot([-72.533], [18.457], 'ro', label='M7.0 epicenter')
	plt.legend(loc='upper left', numpoints=1)
	if doSave: plt.savefig(saveName)
	if doShow: plt.show()
		
def getHaitiRBratio(wLen=10, bigShock=5.2, haitiCat='haiti.cat', doShow=True, maxDt=None):
	# plotting note: #lats=map(operator.itemgetter(1), rbi.shockCat)
	# catFname='parkcat.cat', theta=tTheta, clat=tLat, clon=tLon, ra=tA, rb=tB
	# haiti epicenter bits (so far, very approximate):
	if haitiCat==None: haitiCat='haiti.cat'
	
	hlat=18.46
	#hlon=-72.53	# actual epicenter
	hlon=-72.75		# approximate median position of events.
	haitiTheta=-15.0
	ra=1.25
	rb=.25
	
	#bigShock=5.7	# magnitude of a big aftershock.
	bigShockShocks=[]
	bigFullShocks=[]
	#
	#
	reload(rbi)
	rbHaiti=rbi.intervalRecordBreaker(haitiCat, haitiTheta, hlat, hlon, ra, rb)
	rbHaiti.setAftershockCatalog(haitiCat, haitiTheta, hlat, hlon, ra, rb, None, dtm.datetime.now(), 0)
	#
	#rbHaiti.GRshock(False, fname='images/haitiShock-GRdist.png')
	#rbHaiti.GRfullcat(False, fname='images/haitiFull-GRdist.png')
	#
	# get lists of bigshocks:
	rnum=0
	for rw in rbHaiti.shockCat:
		if rw[3]>=bigShock: bigShockShocks+=[[rnum] + rw]
		rnum+=1
	rnum=0
	for rw in rbHaiti.fullCat:
		if rw[3]>=bigShock: bigFullShocks+=[[rnum] + rw]
		rnum+=1
	#
	#rbHaitiSquare=rbi.intervalRecordBreaker(haitiCat, haitiTheta, hlat, hlon, ra, rb)	# a square catalog. we'll probably just impose a square catalog...
	# find mainshock(s):
	shockMainshockEvNum=0
	fullMainshockEvNum=0
	maxMag1=rbHaiti.shockCat[0][3]
	maxMag2=rbHaiti.fullCat[0][3]
	rwnum1=0
	rwnum2=0
	for rw in rbHaiti.shockCat:
		if rw[3]>maxMag1:
			maxMag1=rw[3]
			shockMainshockEvNum=rwnum1
		rwnum1+=1
	for rw in rbHaiti.fullCat:
		if rw[3]>maxMag2:
			maxMag2=rw[3]
			fullMainshockEvNum=rwnum2
		rwnum2+=1
	#
	# now, get ratios:
	#wLen=10
	hratios=rbHaiti.getIntervalRatios(3.5, wLen, rbHaiti.shockCat)
	if maxDt!=None:
		if maxDt>hratios[-1][1]: hratios+=[[hratios[-1][0]+1, maxDt, hratios[-1][2]]]
	fdts=[]
	for rw in hratios:
		 fdts+=[rw[1].toordinal() + float(rw[1].hour)/24 + rw[1].minute/(24*60) + rw[1].second/(24*3600) + rw[1].microsecond/(24*3600000000)]
	hratiosBoxy=rbHaiti.getIntervalRatios(3.5,wLen,rbHaiti.fullCat)
	if maxDt!=None:
		if maxDt>hratiosBoxy[-1][1]: hratiosBoxy+=[[hratiosBoxy[-1][0]+1, maxDt, hratiosBoxy[-1][2]]]
	fdtsBoxy=[]
	for rw in hratiosBoxy:
		 fdtsBoxy+=[rw[1].toordinal() + float(rw[1].hour)/24 + rw[1].minute/(24*60) + rw[1].second/(24*3600) + rw[1].microsecond/(24*3600000000)]
	#
	eventFloatDt=float(dtm.datetime(2010, 1,12).toordinal()) + 21.0/24.0 + 53.0/(24*60) + 10.0/(24*3600)
	plt.figure(0)
	plt.clf()
	plt.semilogy(fdts, map(operator.itemgetter(2), hratios), 'k.-')
	#plt.fill_between(fdts, map(operator.itemgetter(2), hratios), y2=1, color='b')
	plt.fill([fdts[0]]+ fdts + [fdts[-1], fdts[0]], [1] + yp.getValsAbove(map(operator.itemgetter(2), hratios), 1) + [1, 1], 'b')
	plt.fill([fdts[0]]+ fdts + [fdts[-1], fdts[0]], [1] + yp.getValsBelow(map(operator.itemgetter(2), hratios), 1) + [1, 1], 'r')
	# note: we don't set the y-log scale in these "fill()" commands. we can do that with axis.set_yscale('log') i think.
	# we achieve this by doing semilogy() plots below.
	plt.title("Haiti rupture area, time-time, wLen=%d" % wLen)
	plt.xlabel('time')
	plt.ylabel('$r=N_{rb-long} / N_{rb-short}$')
	#plt.axvline(x=eventFloatDt)
	nbigshocks=0
	for rw in bigShockShocks:
		if nbigshocks==0:
			plt.axvline(x=yp.datetimeToFloat(rw[1]), color='g', label='m > %f' % bigShock)
			nbigshocks+=1
		else:
			plt.axvline(x=yp.datetimeToFloat(rw[1]), color='g')
			
	plt.semilogy([eventFloatDt], [1], 'r^', ms=10)
	plt.axvline(x=eventFloatDt, color='r', lw=3, label='mainshock' )
	#	
	#plt.semilogy([y.datetimeToFloat(rbHaiti.shockCat[25][0])], [1], 'o')
	#plt.axvline(x=yp.datetimeToFloat(rbHaiti.shockCat[25][0]))
	plt.axhline(y=1, color='k')
	plt.legend(loc='upper left', numpoints=2)
	ax=plt.gca()
	fg=plt.gcf()
	ax.xaxis.set_major_formatter(dates.DateFormatter('%Y-%b-%d'))
	fg.autofmt_xdate()
	plt.savefig('images/rbHaitiRuptureTimeTime-Wlen%d.png' % wLen)
	#
	
	plt.figure(1)
	plt.clf()
	#plt.semilogy(fdtsBoxy, map(operator.itemgetter(2), hratiosBoxy), '.-')
	boxyShocks=[[],[]]
	for i in xrange(len(fdtsBoxy)):
		if fdtsBoxy[i]>=eventFloatDt:
			boxyShocks[0]+=[fdtsBoxy[i]]
			boxyShocks[1]+=[hratiosBoxy[i][2]]
			
	#plt.fill([fdtsBoxy[0]] + fdtsBoxy + [fdts[-1], fdtsBoxy[0]], [1]+ yp.getValsAbove(map(operator.itemgetter(2), hratiosBoxy), 1) + [1,1], 'b')
	#plt.fill([fdtsBoxy[0]] + fdtsBoxy + [fdts[-1], fdtsBoxy[0]], [1]+ yp.getValsBelow(map(operator.itemgetter(2), hratiosBoxy), 1) + [1,1], 'r')
	plt.fill([boxyShocks[0][0]] + boxyShocks[0] + [boxyShocks[0][-1], boxyShocks[0][0]], [1]+ yp.getValsAbove(boxyShocks[1], 1) + [1,1], 'b')
	plt.fill([boxyShocks[0][0]] + boxyShocks[0] + [boxyShocks[0][-1], boxyShocks[0][0]], [1]+ yp.getValsBelow(boxyShocks[1], 1) + [1,1], 'r')
	#
	#print "fdts: %s" % str([fdts[0]] + fdtsBoxy + [fdts[-1], fdts[0]])
	#print "vals: %s" % str([1]+ yp.getValsBelow(map(operator.itemgetter(2), hratiosBoxy), 1) + [1,1])
	#plt.semilogy(fdtsBoxy, yp.getValsBelow(map(operator.itemgetter(2), hratiosBoxy), 1), '.-')
	# mainshock:
	
	#plt.axvline(x=eventFloatDt)
	nbigshocks=0
	for rw in bigFullShocks:
		# print "shock events: %s (%f), %f" % (str(rw[1]), yp.datetimeToFloat(rw[1]), eventFloatDt)
		if nbigshocks==0:
			if yp.datetimeToFloat(rw[1])>eventFloatDt:
				plt.axvline(x=yp.datetimeToFloat(rw[1]), color='g', label='m > %f' % bigShock)
				nbigshocks+=1
		else:
			if yp.datetimeToFloat(rw[1])>eventFloatDt: plt.axvline(x=yp.datetimeToFloat(rw[1]), color='g')
		
	plt.semilogy([eventFloatDt], [1], 'r^', ms=10)
	plt.axvline(x=eventFloatDt, color='r', lw=3, label='mainshock' )
	#nbigshokcs=0	
	#for rw in bigShockShocks:
	#	if nbigshocks==0:
	#		plt.axvline(x=yp.datetimeToFloat(rw[1]), color='g', label='m > %f' % bigShock)
	#	else:
	#		plt.axvline(x=yp.datetimeToFloat(rw[1]), color='g')
	#	nbigshocks+=1
		
	plt.axhline(y=1, color='k')
	plt.title("4x4 box around Haiti, time-time Aftershocks, wLen=%d" % wLen)
	plt.xlabel('time')
	plt.ylabel('$r=N_{rb-long} / N_{rb-short}$')
	plt.legend(loc='upper left', numpoints=2)
	ax=plt.gca()
	fg=plt.gcf()
	ax.xaxis.set_major_formatter(dates.DateFormatter('%Y-%b-%d'))
	fg.autofmt_xdate()
	#plt.semilogy([eventFloatDt], [1], 'r^')
	plt.savefig('images/rbHaitiBoxyTimeTimeAftershocks-Wlen%d.png' % wLen)
	#
	plt.figure(2)
	plt.clf()
	plt.semilogy(fdtsBoxy, map(operator.itemgetter(2), hratiosBoxy), 'k.-')
	plt.fill([fdtsBoxy[0]] + fdtsBoxy + [fdts[-1], fdtsBoxy[0]], [1]+ yp.getValsAbove(map(operator.itemgetter(2), hratiosBoxy), 1) + [1,1], 'b')
	plt.fill([fdtsBoxy[0]] + fdtsBoxy + [fdts[-1], fdtsBoxy[0]], [1]+ yp.getValsBelow(map(operator.itemgetter(2), hratiosBoxy), 1) + [1,1], 'r')
	#print "fdts: %s" % str([fdts[0]] + fdtsBoxy + [fdts[-1], fdts[0]])
	#print "vals: %s" % str([1]+ yp.getValsBelow(map(operator.itemgetter(2), hratiosBoxy), 1) + [1,1])
	#plt.semilogy(fdtsBoxy, yp.getValsBelow(map(operator.itemgetter(2), hratiosBoxy), 1), '.-')
	
	#plt.axvline(x=eventFloatDt)
	nbigshocks=0
	for rw in bigFullShocks:
		if nbigshocks==0:
			plt.axvline(x=yp.datetimeToFloat(rw[1]), color='g', label='m > %f' % bigShock)
			nbigshocks+=1
		else:
			plt.axvline(x=yp.datetimeToFloat(rw[1]), color='g')
	
	plt.semilogy([eventFloatDt], [1], 'r^', ms=10)
	plt.axvline(x=eventFloatDt, color='r', lw=3, label='mainshock' )
	
	plt.axhline(y=1, color='k')
	plt.title("4x4 box around Haiti, time-time, wLen=%d" % wLen)
	plt.xlabel('time')
	plt.ylabel('$r=N_{rb-long} / N_{rb-short}$')
	plt.legend(loc='upper left', numpoints=2)
	ax=plt.gca()
	fg=plt.gcf()
	ax.xaxis.set_major_formatter(dates.DateFormatter('%Y-%b-%d'))
	fg.autofmt_xdate()
	#plt.semilogy([eventFloatDt], [1], 'r^')
	plt.savefig('images/rbHaitiBoxyTimeTime-Wlen%d.png' % wLen)
	
	plt.figure(3)
	plt.clf()
	plt.axhline(y=1, color='k')
	plt.semilogy(map(operator.itemgetter(0), hratiosBoxy), map(operator.itemgetter(2), hratiosBoxy), 'k.-')
	X=map(operator.itemgetter(0), hratiosBoxy)
	plt.fill([X[0]] + X + [X[-1], X[0]], [1] + yp.getValsAbove(map(operator.itemgetter(2), hratiosBoxy), 1) + [1,1], 'b')
	plt.fill([X[0]] + X + [X[-1], X[0]], [1] + yp.getValsBelow(map(operator.itemgetter(2), hratiosBoxy), 1) + [1,1], 'r')
	#
	#plt.axvline(x=fullMainshockEvNum)
	nbigshocks=0
	for rw in bigFullShocks:
		if nbigshocks==0:
			plt.axvline(x=rw[0], color='g', label='m > %f' % bigShock)
			nbigshocks+=1
		else:
			plt.axvline(x=rw[0], color='g')
	#
	plt.semilogy([fullMainshockEvNum], [1], 'r^', ms=10)
	plt.axvline(x=fullMainshockEvNum, color='r', lw=3, label='mainshock' )
	#	
	plt.title("4x4 box around Haiti, natural-time, wLen=%d" % wLen)
	plt.xlabel("Number of Events, n")
	plt.ylabel('$r=N_{rb-long} / N_{rb-short}$')
	plt.legend(loc='upper left', numpoints=2)
	# don't over-crowd labels...
	xtks0=map(operator.itemgetter(0), hratiosBoxy)
	xlbls0=map(operator.itemgetter(1), hratiosBoxy)
	nskip=len(xtks0)/15	#15 labels seem to fit pretty well.
	if nskip<1:nskip=1
	xlbls=[]
	for i in xrange(len(xlbls0)):
		lbl=''
		if i%nskip==0: lbl=str(xlbls0[i]) + '(%d)' % i
		xlbls+=[lbl]
	plt.xticks(xtks0, xlbls)
	#ax=plt.gca()
	fg=plt.gcf()
	#ax.xaxis.set_major_formatter(dates.DateFormatter('%Y-%b-%d'))
	fg.autofmt_xdate()
	plt.axis('auto')
	plt.savefig('images/rbHaitiBoxyNaturalTime-Wlen%d.png' % wLen)
	
	plt.figure(4)
	plt.clf()
	plt.axhline(y=1, color='k')
	plt.semilogy(map(operator.itemgetter(0), hratios), map(operator.itemgetter(2), hratios), 'k.-')
	X=map(operator.itemgetter(0), hratios)
	plt.fill([X[0]] + X+[X[-1], X[0]], [1] + yp.getValsAbove(map(operator.itemgetter(2), hratios), 1) + [1,1], 'b')
	plt.fill([X[0]] + X+[X[-1], X[0]], [1] + yp.getValsBelow(map(operator.itemgetter(2), hratios), 1) + [1,1], 'r')
	
	#plt.semilogy([25], [1], 'o')
	#plt.axvline(x=25)
	nbigshocks=0
	for rw in bigShockShocks:
		if nbigshocks==0:
			plt.axvline(x=rw[0], color='g', label='m > %f' % bigShock)
			nbigshocks+=1
		else:
			plt.axvline(x=rw[0], color='g')
	plt.semilogy([shockMainshockEvNum], [1], 'r^', ms=10)
	plt.axvline(x=shockMainshockEvNum, color='r', lw=3, label='mainshock')
		
	plt.title("Haiti rupture area, natural-time, wLen=%d" % wLen)
	plt.xlabel("Number of Events, n")
	plt.ylabel('$r=N_{rb-long} / N_{rb-short}$')
	# don't over-crowd labels...
	xtks0=map(operator.itemgetter(0), hratios)
	xlbls0=map(operator.itemgetter(1), hratios)
	nskip=len(xtks0)/9	#15 labels seem to fit pretty well.
	if nskip<1:nskip=1
	xlbls=[]
	for i in xrange(len(xlbls0)):
		lbl=''
		if i%nskip==0: lbl=str(xlbls0[i]) + '(%d)' % i
		xlbls+=[lbl]
	plt.xticks(xtks0, xlbls)
	#ax=plt.gca()
	fg=plt.gcf()
	#ax.xaxis.set_major_formatter(dates.DateFormatter('%Y-%b-%d'))
	fg.autofmt_xdate()
	plt.legend(loc='upper left', numpoints=2)
	plt.savefig('images/rbHaitiRuptureNaturalTime-Wlen%d.png' % wLen)
	
	if doShow: plt.show()
	
	# note: look at the shockCat record-breaking intervals plot; we pretty well forecast the 5.75 aftershock (event 25)
	
	return [rbHaiti, hratios]

##################################################################
# chile stuff:


def getChileCatFromSQL(outFile='chile.cat', startDate=dtm.datetime(1980, 01, 01), endDate=dtm.datetime.now(), minmag=5.0, lats=[-40, -30], lons=[-80.0, -65.0], catID=523 ):
	# haiti catalog is in flux, so let's use the catalog. we could go directly from ANSS if we like, but i think this will be easier.
	# note that the default prams are a bit broad, but we seem to get a signal anyway.
	#
	sqlSelect = "select eventDateTime, lat, lon, mag from Earthquakes where catalogID=%d and eventDateTime between '%s' and '%s' and lat between %f and %f and lon between %f and %f and mag>=%f order by eventDateTime asc" % (catID, str(startDate), str(endDate), lats[0], lats[1], lons[0], lons[1], minmag)
	sqlHost = 'localhost'
	sqlUser = 'myoder'
	sqlPassword = 'yoda'
	sqlPort = 3306
	sqlDB = 'QuakeData'
	myConn = MySQLdb.connect(host=sqlHost, user=sqlUser, passwd=sqlPassword, port=sqlPort, db=sqlDB)
	r1=myConn.cursor()
	r1.execute(sqlSelect)
	#	
	fout=open(outFile, 'w')
	fout.write("# chile catalog from hatistuff.py\n#prams: outfile, startDate, endDate, minmag, lats, lons, catID\n")
	fout.write("#%s\t%s\t%s\t%s\t%s\t%s\t%s\d\n" % (outFile, startDate, endDate, minmag, lats, lons, catID))
	# catalog format (tab delimited):
	# 2004/09/27	01:53:42.78	35.5277	-120.8433	1.42
	catlen=0
	for rw in r1:
		fout.write("%d/%d/%d\t%d:%d:%d.%d\t%f\t%f\t%f\n" % (rw[0].year, rw[0].month, rw[0].day, rw[0].hour, rw[0].minute, rw[0].second, rw[0].microsecond, rw[1], rw[2], rw[3]))
		catlen+=1
	fout.close()
	r1.close()
	myConn.close()
	
	return catlen

def getStandardchileRB(hlat=-35.846, hlon=-72.719, chileTheta=-55, ra=4.5, rb=3.75, chileCat='chile.cat'):
	reload(rbi)
	rbchile=rbi.intervalRecordBreaker(chileCat, chileTheta, hlat, hlon, ra, rb)
	rbchile.setAftershockCatalog(chileCat, chileTheta, hlat, hlon, ra, rb, None, dtm.datetime.now(), 0)
	return rbchile

def freshChileRatios(wLens=[15, 20, 30, 40, 50], bigShock=6.0, haitiCat='chile.cat', doShow=False, minMag=5.0):
	# full set of Haiti rb ratios from a fresh catalog...
	# update catalog:
	print "update MySQL catalog..."
	yp.updateANSS2SQL()
	print "update local catalog, %s" % haitiCat
	getChileCatFromSQL(chileCat)
	#
	return chileRatios(wLens, bigShock, haitiCat, doShow, minMag)
	
def chileRatios(wLens=[15, 20, 30, 40, 50], bigShock=6.0, chileCat='chile.cat', doShow=False, minMag=5.0):
	# full set of Haiti rb ratios
	# plot GR distributions:
	rbChile=getStandardchileRB()
	rbChile.GRshock(False, fname='images/chileShock-GRdist.png')
	rbChile.GRfullcat(False, fname='images/chileFull-GRdist.png')
	rbChile=None
	
	# calc ratios for wLens:
	print "calc ratios for wLens"
	for wLen in wLens:
		print "rb ratios for wLen=%d" % wLen
		getChileRBratios(wLen, bigShock, chileCat, doShow, dtm.datetime.now(), minMag)
	#
	return 0

def getChileCatPlot(doShow=False, doSave=True, rbh=None, saveName='images/chileCatalogPlot.png'):
	# rbh is a Chile specific (derived from haiti specific) intervalRecordBreaker() object.
	# note: this functionality has been integrated into the intervalRecordBreaker() object as obh.xyPlotCatalogs(doShow, doSave, saveName, ellipseCenter=[x,y])
	if rbh==None: rbh=getStandardchileRB()
	
	bangDate=dtm.datetime(2010, 2, 27, 6, 34, 14)
	
	fcat=[]
	scat=[]
	aftershocks=[]
	for rw in rbh.shockCat:
		scat+=[rw[0:4]]
	for rw in rbh.fullCat:
		if rw not in scat: fcat+=[rw]
		if rw[0]>=bangDate: aftershocks+=[rw]
	
	#print "aftershocks: %d" % len(aftershocks)
	#return [scat, fcat]
	
	plt.figure(0)	
	plt.clf()
	#
	ax=plt.gca()
	
	tLat=None # 35.9
	tLon=None # -120.5
	tTheta=None # 40.0		#47?	note: tTheta is the angle CCW of the x' (transformed) axis from the x axis.
	tA=None # .4		# ellipse axes
	tB=None # .15
	
	#el = Ellipse((-72.533,18.457), 1.25, .25, 15, facecolor='r', alpha=0.5)
	#el = Ellipse((-72.719,-35.846), 2.0*rbh.tA, 2.0*rbh.tB, rbh.tTheta, facecolor='b', alpha=0.4)
	el = Ellipse((rbh.tLon, rbh.tLat), 2.0*rbh.tA, 2.0*rbh.tB, -rbh.tTheta, facecolor='b', alpha=0.3)
	ax.add_artist(el)
	#
	#plt.plot(map(operator.itemgetter(2), rbh.fullCat), map(operator.itemgetter(1), rbh.fullCat), '.')
	plt.plot(map(operator.itemgetter(2), fcat), map(operator.itemgetter(1), fcat), 'b.')
	plt.plot(map(operator.itemgetter(2), scat), map(operator.itemgetter(1), scat), 'g.')
	#plt.plot(map(operator.itemgetter(2), rbh.shockCat), map(operator.itemgetter(1), rbh.shockCat), '.')
	plt.plot(map(operator.itemgetter(2), aftershocks), map(operator.itemgetter(1), aftershocks), 'r.')
	
	#plt.plot(map(operator.itemgetter(2), fcat), map(operator.itemgetter(1), fcat), '+', label='Full Catalog')
	#plt.plot(map(operator.itemgetter(2), scat), map(operator.itemgetter(1), scat), '.', label='Aftershock zone')
	plt.plot([rbh.tLon], [rbh.tLat], 'r*', ms=15, label='M8.8 epicenter')
	plt.legend(loc='upper left', numpoints=1)
	if doSave: plt.savefig(saveName)
	if doShow: plt.show()
	
	return rbh

def getChileRBratios(wLen=25, bigShock=5.2, chileCat='chile.cat', doShow=True, maxDt=None, minMag=5.0):
	# plotting note: #lats=map(operator.itemgetter(1), rbi.shockCat)
	# catFname='parkcat.cat', theta=tTheta, clat=tLat, clon=tLon, ra=tA, rb=tB
	# haiti epicenter bits (so far, very approximate):
	if chileCat==None: chileCat='chile.cat'
	
	#bigShock=5.7	# magnitude of a big aftershock.
	bigShockShocks=[]
	bigFullShocks=[]
	#
	#
	reload(rbi)
	rbchile=getStandardchileRB()
	hlat=rbchile.tLat
	#hlon=-72.53	# actual epicenter
	hlon=rbchile.tLon
	chileTheta=rbchile.tTheta
	ra=rbchile.tA
	rb=rbchile.tB
	print "prams: %f, %f, %f, %f, %f" % (hlat, hlon, chileTheta, ra, rb)
	# hlat=-35.846, hlon=-72.719, chileTheta=-55, ra=4.5, rb=3.75
	rbchile.setAftershockCatalog(chileCat, chileTheta, hlat, hlon, ra, rb, None, dtm.datetime.now(), 0)
	#
	rbchile.GRshock(False, fname='images/chileShock-GRdist.png')
	rbchile.GRfullcat(False, fname='images/chileFull-GRdist.png')
	#
	# get lists of bigshocks:
	rnum=0
	for rw in rbchile.shockCat:
		if rw[3]<minMag: continue
		if rw[3]>=bigShock: bigShockShocks+=[[rnum] + rw]
		rnum+=1
	rnum=0
	for rw in rbchile.fullCat:
		if rw[3]<minMag: continue
		if rw[3]>=bigShock: bigFullShocks+=[[rnum] + rw]
		rnum+=1
	#
	#rbchileSquare=rbi.intervalRecordBreaker(chileCat, chileTheta, hlat, hlon, ra, rb)	# a square catalog. we'll probably just impose a square catalog...
	# find mainshock(s):
	shockMainshockEvNum=0
	fullMainshockEvNum=0
	maxMag1=rbchile.shockCat[0][3]
	maxMag2=rbchile.fullCat[0][3]
	rwnum1=0
	rwnum2=0
	for rw in rbchile.shockCat:
		if rw[3]<minMag: continue
		if rw[3]>maxMag1:
			maxMag1=rw[3]
			shockMainshockEvNum=rwnum1
			eventFloatDt=yp.datetimeToFloat(rw[0])
		rwnum1+=1
	#print "event date: %f" % (eventFloatDt)
	for rw in rbchile.fullCat:
		if rw[3]<minMag: continue
		if rw[3]>maxMag2:
			maxMag2=rw[3]
			fullMainshockEvNum=rwnum2
		rwnum2+=1
	print "mainshock data: mag=%f/%f, evNum1=%d, evNum2=%d" % (maxMag1, maxMag2, shockMainshockEvNum, fullMainshockEvNum)
	#
	# now, get ratios:
	#wLen=10
		
	hratios=rbchile.getIntervalRatios(minMag, wLen, rbchile.shockCat)
	if maxDt!=None:
		if maxDt>hratios[-1][1]: hratios+=[[hratios[-1][0]+1, maxDt, hratios[-1][2]]]
	fdts=[]
	for rw in hratios:
		 fdts+=[rw[1].toordinal() + float(rw[1].hour)/24 + rw[1].minute/(24*60) + rw[1].second/(24*3600) + rw[1].microsecond/(24*3600000000)]
	hratiosBoxy=rbchile.getIntervalRatios(minMag,wLen,rbchile.fullCat)
	if maxDt!=None:
		if maxDt>hratiosBoxy[-1][1]: hratiosBoxy+=[[hratiosBoxy[-1][0]+1, maxDt, hratiosBoxy[-1][2]]]
	fdtsBoxy=[]
	for rw in hratiosBoxy:
		 fdtsBoxy+=[rw[1].toordinal() + float(rw[1].hour)/24 + rw[1].minute/(24*60) + rw[1].second/(24*3600) + rw[1].microsecond/(24*3600000000)]
	
	
	plt.figure(0)
	plt.clf()
	#plt.semilogy(fdts, map(operator.itemgetter(2), hratios), 'k.')
	#plt.fill_between(fdts, map(operator.itemgetter(2), hratios), y2=1, color='b')
	
	#plt.fill([fdts[0]]+ fdts + [fdts[-1], fdts[0]], [1] + yp.getValsAbove(map(operator.itemgetter(2), hratios), 1) + [1, 1], 'b')
	#plt.fill([fdts[0]]+ fdts + [fdts[-1], fdts[0]], [1] + yp.getValsBelow(map(operator.itemgetter(2), hratios), 1) + [1, 1], 'r')
	
	plt.axvline(x=eventFloatDt, color='c', lw=3, label='mainshock' )
	plt.fill_between(fdts, scipy.ones(len(fdts),int), map(operator.itemgetter(2), hratios), color='b', where=scipy.array([val>=1 for val in map(operator.itemgetter(2), hratios)]))
	plt.fill_between(fdts, scipy.ones(len(fdts),int), map(operator.itemgetter(2), hratios), color='r', where=scipy.array([val<=1 for val in map(operator.itemgetter(2), hratios)]))
	
	# note: we don't set the y-log scale in these "fill()" commands. we can do that with axis.set_yscale('log') i think.
	# we achieve this by doing semilogy() plots below.
	plt.title("Chile rupture area, time-time, wLen=%d" % wLen)
	plt.xlabel('time')
	plt.ylabel('$r=N_{rb-long} / N_{rb-short}$')
	plt.axvline(x=eventFloatDt)
	nbigshocks=0
	for rw in bigShockShocks:
		if nbigshocks==0:
			plt.axvline(x=yp.datetimeToFloat(rw[1]), color='g', label='m > %f' % bigShock)
			nbigshocks+=1
		else:
			plt.axvline(x=yp.datetimeToFloat(rw[1]), color='g')
			
	plt.semilogy([eventFloatDt], [1], 'r^', ms=10)
	#plt.axvline(x=eventFloatDt, color='r', lw=3, label='mainshock' )
	#	
	#plt.semilogy([y.datetimeToFloat(rbchile.shockCat[25][0])], [1], 'o')
	#plt.axvline(x=yp.datetimeToFloat(rbchile.shockCat[25][0]))
	plt.axhline(y=1, color='k')
	plt.legend(loc='lower right', numpoints=2)
	ax=plt.gca()
	fg=plt.gcf()
	ax.xaxis.set_major_formatter(dates.DateFormatter('%Y-%b-%d'))
	ax.set_ylim([.1,10])
	fg.autofmt_xdate()
	plt.savefig('images/rbchileRuptureTimeTime-Wlen%d.png' % wLen)
	#
	
	plt.figure(1)
	plt.clf()
	#plt.semilogy(fdtsBoxy, map(operator.itemgetter(2), hratiosBoxy), '.-')
	boxyShocks=[[],[]]
	for i in xrange(len(fdtsBoxy)):
		if fdtsBoxy[i]>=eventFloatDt:
			boxyShocks[0]+=[fdtsBoxy[i]]
			boxyShocks[1]+=[hratiosBoxy[i][2]]
			
	#plt.fill([fdtsBoxy[0]] + fdtsBoxy + [fdts[-1], fdtsBoxy[0]], [1]+ yp.getValsAbove(map(operator.itemgetter(2), hratiosBoxy), 1) + [1,1], 'b')
	#plt.fill([fdtsBoxy[0]] + fdtsBoxy + [fdts[-1], fdtsBoxy[0]], [1]+ yp.getValsBelow(map(operator.itemgetter(2), hratiosBoxy), 1) + [1,1], 'r')
	
	#plt.fill([boxyShocks[0][0]] + boxyShocks[0] + [boxyShocks[0][-1], boxyShocks[0][0]], [1]+ yp.getValsAbove(boxyShocks[1], 1) + [1,1], 'b')
	#plt.fill([boxyShocks[0][0]] + boxyShocks[0] + [boxyShocks[0][-1], boxyShocks[0][0]], [1]+ yp.getValsBelow(boxyShocks[1], 1) + [1,1], 'r')
	
	plt.fill_between(boxyShocks[0], scipy.ones(len(boxyShocks[0]),int), boxyShocks[1], color='b', where=scipy.array([val>=1 for val in boxyShocks[1]]))
	plt.fill_between(boxyShocks[0], scipy.ones(len(boxyShocks[0]),int), boxyShocks[1], color='r', where=scipy.array([val<=1 for val in boxyShocks[1]]))

	
	#plt.bar([boxyShocks[0][0]] + boxyShocks[0] + [boxyShocks[0][-1], boxyShocks[0][0]], [1]+ yp.getValsAbove(boxyShocks[1], 1) + [1,1], width=1, color='b', orientation='vertical')
	#plt.bar([boxyShocks[0][0]] + boxyShocks[0] + [boxyShocks[0][-1], boxyShocks[0][0]], [1]+ yp.getValsBelow(boxyShocks[1], 1) + [1,1], width=1, color='r', orientation='vertical')
	#
	#print "fdts: %s" % str([fdts[0]] + fdtsBoxy + [fdts[-1], fdts[0]])
	#print "vals: %s" % str([1]+ yp.getValsBelow(map(operator.itemgetter(2), hratiosBoxy), 1) + [1,1])
	#plt.semilogy(fdtsBoxy, yp.getValsBelow(map(operator.itemgetter(2), hratiosBoxy), 1), '.-')
	# mainshock:
	
	#plt.axvline(x=eventFloatDt)
	plt.axvline(x=eventFloatDt, color='c', lw=3, label='mainshock')
	nbigshocks=0
	for rw in bigFullShocks:
		# print "shock events: %s (%f), %f" % (str(rw[1]), yp.datetimeToFloat(rw[1]), eventFloatDt)
		if nbigshocks==0:
			if yp.datetimeToFloat(rw[1])>eventFloatDt:
				plt.axvline(x=yp.datetimeToFloat(rw[1]), color='g', label='m > %f' % bigShock)
				nbigshocks+=1
		else:
			if yp.datetimeToFloat(rw[1])>eventFloatDt: plt.axvline(x=yp.datetimeToFloat(rw[1]), color='g')
		
	plt.semilogy([eventFloatDt], [1], 'r^', ms=10)
	#plt.axvline(x=eventFloatDt, color='r', lw=3, label='mainshock' )
	#nbigshokcs=0	
	#for rw in bigShockShocks:
	#	if nbigshocks==0:
	#		plt.axvline(x=yp.datetimeToFloat(rw[1]), color='g', label='m > %f' % bigShock)
	#	else:
	#		plt.axvline(x=yp.datetimeToFloat(rw[1]), color='g')
	#	nbigshocks+=1
		
	plt.axhline(y=1, color='k')
	plt.title("Box around Chile, time-time Aftershocks, wLen=%d" % wLen)
	plt.xlabel('time')
	plt.ylabel('$r=N_{rb-long} / N_{rb-short}$')
	plt.legend(loc='lower right', numpoints=2)
	ax=plt.gca()
	fg=plt.gcf()
	ax.xaxis.set_major_formatter(dates.DateFormatter('%Y-%b-%d'))
	ax.set_ylim([.1,10])
	fg.autofmt_xdate()
	#plt.semilogy([eventFloatDt], [1], 'r^')
	plt.savefig('images/rbchileBoxyTimeTimeAftershocks-Wlen%d.png' % wLen)
	#
	plt.figure(2)
	plt.clf()
	#plt.semilogy(fdtsBoxy, map(operator.itemgetter(2), hratiosBoxy), 'k.')
	#plt.fill([fdtsBoxy[0]] + fdtsBoxy + [fdts[-1], fdtsBoxy[0]], [1]+ yp.getValsAbove(map(operator.itemgetter(2), hratiosBoxy), 1) + [1,1], 'b')
	#plt.fill([fdtsBoxy[0]] + fdtsBoxy + [fdts[-1], fdtsBoxy[0]], [1]+ yp.getValsBelow(map(operator.itemgetter(2), hratiosBoxy), 1) + [1,1], 'r')
	
	plt.fill_between(fdtsBoxy, scipy.ones(len(fdtsBoxy),int), map(operator.itemgetter(2), hratiosBoxy), color='b', where=scipy.array([val>=1 for val in map(operator.itemgetter(2), hratiosBoxy)]))
	plt.fill_between(fdtsBoxy, scipy.ones(len(fdtsBoxy),int), map(operator.itemgetter(2), hratiosBoxy), color='r', where=scipy.array([val<=1 for val in map(operator.itemgetter(2), hratiosBoxy)]))

	#print "fdts: %s" % str([fdts[0]] + fdtsBoxy + [fdts[-1], fdts[0]])
	#print "vals: %s" % str([1]+ yp.getValsBelow(map(operator.itemgetter(2), hratiosBoxy), 1) + [1,1])
	#plt.semilogy(fdtsBoxy, yp.getValsBelow(map(operator.itemgetter(2), hratiosBoxy), 1), '.-')
	
	#plt.axvline(x=eventFloatDt)
	plt.axvline(x=eventFloatDt, color='c', lw=3, label='mainshock' )
	nbigshocks=0
	for rw in bigFullShocks:
		if nbigshocks==0:
			plt.axvline(x=yp.datetimeToFloat(rw[1]), color='g', label='m > %f' % bigShock)
			nbigshocks+=1
		else:
			plt.axvline(x=yp.datetimeToFloat(rw[1]), color='g')
	
	plt.semilogy([eventFloatDt], [1], 'r^', ms=10)
	#plt.axvline(x=eventFloatDt, color='r', lw=3, label='mainshock' )
	
	plt.axhline(y=1, color='k')
	plt.title("Box around Chile, time-time, wLen=%d" % wLen)
	plt.xlabel('time')
	plt.ylabel('$r=N_{rb-long} / N_{rb-short}$')
	plt.legend(loc='upper left', numpoints=2)
	ax=plt.gca()
	fg=plt.gcf()
	ax.xaxis.set_major_formatter(dates.DateFormatter('%Y-%b-%d'))
	ax.set_ylim([.1,10])
	fg.autofmt_xdate()
	#plt.semilogy([eventFloatDt], [1], 'r^')
	plt.savefig('images/rbchileBoxyTimeTime-Wlen%d.png' % wLen)
	
	plt.figure(3)
	plt.clf()
	plt.axhline(y=1, color='k')
	#
	plt.axvline(x=fullMainshockEvNum, color='c', lw=3, label='mainshock' )
	#
	#plt.semilogy(map(operator.itemgetter(0), hratiosBoxy), map(operator.itemgetter(2), hratiosBoxy), 'k.')
	X=map(operator.itemgetter(0), hratiosBoxy)
	#plt.fill([X[0]] + X + [X[-1], X[0]], [1] + yp.getValsAbove(map(operator.itemgetter(2), hratiosBoxy), 1) + [1,1], 'b')
	#plt.fill([X[0]] + X + [X[-1], X[0]], [1] + yp.getValsBelow(map(operator.itemgetter(2), hratiosBoxy), 1) + [1,1], 'r')
	plt.fill_between(X, scipy.ones(len(X),int), map(operator.itemgetter(2), hratiosBoxy), color='b', where=scipy.array([val>=1 for val in map(operator.itemgetter(2), hratiosBoxy)]))
	plt.fill_between(X, scipy.ones(len(X),int), map(operator.itemgetter(2), hratiosBoxy), color='r', where=scipy.array([val<=1 for val in map(operator.itemgetter(2), hratiosBoxy)]))
	#
	#plt.axvline(x=fullMainshockEvNum)
	nbigshocks=0
	for rw in bigFullShocks:
		if nbigshocks==0:
			plt.axvline(x=rw[0], color='g', label='m > %f' % bigShock)
			nbigshocks+=1
		else:
			plt.axvline(x=rw[0], color='g')
	#
	plt.semilogy([fullMainshockEvNum], [1], 'r^', ms=10)
	#plt.axvline(x=fullMainshockEvNum, color='r', lw=3, label='mainshock' )
	#	
	plt.title("Box around Chile, natural-time, wLen=%d" % wLen)
	plt.xlabel("Number of Events, n")
	plt.ylabel('$r=N_{rb-long} / N_{rb-short}$')
	plt.legend(loc='upper left', numpoints=2)
	# don't over-crowd labels...
	xtks0=map(operator.itemgetter(0), hratiosBoxy)
	xlbls0=map(operator.itemgetter(1), hratiosBoxy)
	nskip=len(xtks0)/15	#15 labels seem to fit pretty well.
	if nskip<1:nskip=1
	xlbls=[]
	for i in xrange(len(xlbls0)):
		lbl=''
		if i%nskip==0: lbl=str(xlbls0[i]) + '(%d)' % i
		xlbls+=[lbl]
	plt.xticks(xtks0, xlbls)
	ax=plt.gca()
	fg=plt.gcf()
	ax.set_ylim([.1,10])
	#ax.xaxis.set_major_formatter(dates.DateFormatter('%Y-%b-%d'))
	fg.autofmt_xdate()
	plt.axis('auto')
	plt.savefig('images/rbchileBoxyNaturalTime-Wlen%d.png' % wLen)
	
	plt.figure(4)
	plt.clf()
	plt.axhline(y=1, color='k')
	#plt.semilogy(map(operator.itemgetter(0), hratios), map(operator.itemgetter(2), hratios), 'k.')
	X=map(operator.itemgetter(0), hratios)
	#plt.fill([X[0]] + X+[X[-1], X[0]], [1] + yp.getValsAbove(map(operator.itemgetter(2), hratios), 1) + [1,1], 'b')
	#plt.fill([X[0]] + X+[X[-1], X[0]], [1] + yp.getValsBelow(map(operator.itemgetter(2), hratios), 1) + [1,1], 'r')
	#plt.bar(X, map(operator.itemgetter(2), hratios), width=1, bottom=1, color='b', align='edge')
	#plt.bar([X[0]] + X+[X[-1], X[0]], [1] + yp.getValsBelow(map(operator.itemgetter(2), hratios), 1) + [1,1], width=1, color='r', orientation='vertical', align='edge')
	plt.axvline(x=shockMainshockEvNum, color='c', lw=3, label='mainshock')
	plt.fill_between(X, scipy.ones(len(X),int), map(operator.itemgetter(2), hratios), color='b', where=scipy.array([val>=1 for val in map(operator.itemgetter(2), hratios)]))
	plt.fill_between(X, scipy.ones(len(X),int), map(operator.itemgetter(2), hratios), color='r', where=scipy.array([val<1 for val in map(operator.itemgetter(2), hratios)]))
	#plt.semilogy([25], [1], 'o')
	#plt.axvline(x=25)
	nbigshocks=0
	for rw in bigShockShocks:
		if nbigshocks==0:
			plt.axvline(x=rw[0], color='g', label='m > %f' % bigShock)
			nbigshocks+=1
		else:
			plt.axvline(x=rw[0], color='g')
	plt.semilogy([shockMainshockEvNum], [1], 'r^', ms=10)
	#plt.axvline(x=shockMainshockEvNum, color='r', lw=3, label='mainshock')
		
	plt.title("Chile rupture area, natural-time, wLen=%d" % wLen)
	plt.xlabel("Number of Events, n")
	plt.ylabel('$r=N_{rb-long} / N_{rb-short}$')
	# don't over-crowd labels...
	xtks0=map(operator.itemgetter(0), hratios)
	xlbls0=map(operator.itemgetter(1), hratios)
	nskip=len(xtks0)/9	#15 labels seem to fit pretty well.
	if nskip<1:nskip=1
	xlbls=[]
	for i in xrange(len(xlbls0)):
		lbl=''
		if i%nskip==0: lbl=str(xlbls0[i]) + '(%d)' % i
		xlbls+=[lbl]
	plt.xticks(xtks0, xlbls)
	ax=plt.gca()
	ax.set_ylim([.1,10])
	fg=plt.gcf()
	#ax.xaxis.set_major_formatter(dates.DateFormatter('%Y-%b-%d'))
	fg.autofmt_xdate()
	plt.legend(loc='upper left', numpoints=2)
	plt.savefig('images/rbchileRuptureNaturalTime-Wlen%d.png' % wLen)
	
	if doShow: plt.show()
	
	# note: look at the shockCat record-breaking intervals plot; we pretty well forecast the 5.75 aftershock (event 25)
	
	return [rbchile, hratios]


##################################################################
# mexicali 2010:
def getMexicaliCat(mc=1.5):
	# hlat=32.258, hlon=-115.287, Theta=50.0, ra=2.0, rb=.5, thisCat='mexicali.cat'
	eqm=yp.eqcatalog([])
	eqm.rb=rbi.intervalRecordBreaker(None)
	xcenter=32.258
	ycenter=-115.287
	xsize=2.5
	ysize=2.5
	#startDt=dtm.datetime(2000,1,1)
	startDt=dtm.datetime(2000,1,1)
	eqm.setCatFromSQL(startDt, dtm.datetime.today(), lats=[xcenter-xsize, xcenter+xsize], lons=[ycenter-ysize, ycenter+ysize], minmag=mc, catalogName='Earthquakes', catalogID=523, ordering='asc')
	#
	#subcats: start with an xy limited cat:
	# addLatLonSubcat(self, subcatname='xysubcat', fullcat=None, lats=[], lons=[], llcols=[1,2]):
	# getLatLonSubcat(fullcat, lats, lons, llcols)
	eqm.addLatLonSubcat('xysubcat1', eqm.cat, [31.5, 33], [-116.25, -114.5], [1,2])
	eqm.addLatLonSubcat('xysubcat2', eqm.cat, [32, 32.5], [-115.5, -115])
	
	# eqm.addTimeRangeCat('postMS', eqm.getcat(0), dtm.datetime(2010, 01, 01), dtm.datetime.today())
	#eqm.addEllipCat('PFshock (.8 x .15)', eqm.cat, 40.0, 35.9, -120.5, 0.8, 0.15) # (parkfield reff.)
	eqm.addEllipCat('ellip1', eqm.getcat(0), 40.0, 32.258, -115.287, 1.0, .25)
	#
	return eqm

def getSocalCat(mc=1.5):
	# hlat=32.258, hlon=-115.287, Theta=50.0, ra=2.0, rb=.5, thisCat='mexicali.cat'
	eqm=yp.eqcatalog([])
	eqm.rb=rbi.intervalRecordBreaker(None)
	#xcenter=32.258
	#ycenter=-115.287
	#xsize=2.5
	#ysize=2.5
	lats=[31.5, 36.5]
	lons=[-122.5, -114.5]
	#
	#startDt=dtm.datetime(2000,1,1)
	startDt=dtm.datetime(1995,1,1)
	eqm.setCatFromSQL(startDt, dtm.datetime.today(), lats, lons, minmag=mc, catalogName='Earthquakes', catalogID=523, ordering='asc')
	#
	#subcats: start with an xy limited cat:
	# addLatLonSubcat(self, subcatname='xysubcat', fullcat=None, lats=[], lons=[], llcols=[1,2]):
	# getLatLonSubcat(fullcat, lats, lons, llcols)
	#
	# mexicali equakes:
	eqm.addLatLonSubcat('xysubcat1', eqm.cat, [31.5, 33], [-116.25, -114.5], [1,2])
	eqm.addLatLonSubcat('xysubcat2', eqm.cat, [32, 32.5], [-115.5, -115])
	
	# eqm.addTimeRangeCat('postMS', eqm.getcat(0), dtm.datetime(2010, 01, 01), dtm.datetime.today())
	#eqm.addEllipCat('PFshock (.8 x .15)', eqm.cat, 40.0, 35.9, -120.5, 0.8, 0.15) # (parkfield reff.)
	eqm.addEllipCat('ellip1', eqm.getcat(0), 40.0, 32.258, -115.287, 1.0, .25)
	#
	return eqm

def rbHazardMap(catalog=None, catnum=0, mc=2.5, gridsize=.5, phi=[0,0], winlen=128, bigmag=5.0, fignum=0):
	# make a hazard map from rb-sequence.
	# divide cat. into gridsize size square subcats, do RB over last winlen events in each square, draw on a map.
	#
	# phi=[xoffset, yoffset]. we'll accomplish this by simply adding phi to the LL x,y values. nominally this changes the area of the catalog a little bit,
	# but we don't throw away any data and we sill shift the bin centers.
	if catalog==None:
		#catalog=getMexicaliCat(mc)
		catalog=getSocalCat(mc)
		while len(catalog.subcats)>0: catalog.subcats.pop()
	avlen=10
	#
	thiscat0=catalog.getcat(catnum)
	mev=catalog.getMainEvent(thiscat0)
	thiscat=thiscat0[0:mev[4]-1]
	print "mainEvent: %s" % mev
	#
	llrange=catalog.getLatLonRange(thiscat)
	llrange[0]=[llrange[0][0]+phi[0], llrange[0][1]+phi[1]]
	deltaLat=abs(float(llrange[1][0])-llrange[0][0])
	deltaLon=abs(float(llrange[1][1])-llrange[0][1])
	Nlat=1+int(deltaLat/gridsize)
	Nlon=1+int(deltaLon/gridsize)	#number of cells in row/col
	ilat=0
	ilon=0	# index counters
	# add empty sub-cats:
	for i in xrange(Nlat*Nlon):
		catalog.subcats+=[[i, []]]
	#
	#return catalog
	lat0=llrange[0][0]
	lon0=llrange[0][1]
	#imax=len(thiscat)-1
	# now, skip through main-catalog, get cat-index for each event and append...
	#for icat in xrange(imax):
	#print "stats: %d, %d, (%f), %f, %f" % (Nlat, Nlon, gridsize, lat0, lon0)
	for rw in thiscat:
		#rw=thiscat[imax-icat]	# move backwards, from most recent to past.
	#	rw=thiscat[icat]
		#
		ilat=(rw[1]-lat0)/gridsize
		ilat=int(ilat)
		ilon=(rw[2]-lon0)/gridsize
		ilon=int(ilon)
		scindex=ilat*Nlon + ilon	# sub-catalog index
		#
		#print "lat, lon, ilat, ilon: [%d]: %f, %f, %d, %d (%d, %d, (%f), %f, %f)" % (scindex, rw[1], rw[2], ilat, ilon, Nlat, Nlon, gridsize, lat0, lon0)
		#
		#print "index: %d of %d/%d (%s)" % (scindex, len(catalog.subcats), imax, rw)
		#if len(catalog.subcats[scindex][1])<=(winlen+1): catalog.subcats[scindex][1].insert(0, rw)
		
		#catalog.subcats[scindex][1].insert(0, rw)
		catalog.subcats[scindex][1]+=[rw]
	#
	# now, calc. statistics from each square.
	# [index, totalInterval, r_nrb]
	gridStats=[]

	#return catalog
	for ct in catalog.subcats:
		if len(ct[1])<=(winlen+avlen): continue
		#print ct[0]
		#print ct[1][-1][0]
		#print ct[1][0][0]
		#
		thisratios=catalog.rb.getIntervalRatios(mc, winlen, ct[1])
		rs=map(operator.itemgetter(-1), thisratios)
		meanr=sum(rs[-avlen:])/float(avlen)
		#
		#gridStats+=[[ct[0], pylab.date2num(ct[1][-1][0])-pylab.date2num(ct[1][0][0]), catalog.rb.getIntervalRatios(mc, winlen, ct[1])[-1][-1]]]
		gridStats+=[[ct[0], pylab.date2num(ct[1][-1][0])-pylab.date2num(ct[1][0][0]), meanr]]
	#
	#plt.figure(5)
	reds=[]
	blues=[]
	print "make grid elements:"
	print llrange
	print lat0, lon0
	#plt.title('as of %s' % str(thiscat[-1][0]))
	plt.text(.4, .8, 'as of %s' % str(thiscat[-1][0]))
	print 'as of %s' % str(thiscat[-1][0])
	cm=catalog.plotCatMap(thiscat, True, False, None, None, 'best', True, 'g,', None, 1)
	for rw in gridStats:
		#
		clat=int(rw[0]/Nlon)*gridsize+lat0
		clon=int(rw[0]%Nlon)*gridsize+lon0
		#print rw, clat, clon, Nlon, gridsize
		#
		if rw[2]<=1:
			reds+=[[clat, clon, rw[1], rw[2]]]
			thispoly=getSquare([clon, clat], gridsize)
			rx, ry=cm(x=map(operator.itemgetter(0), thispoly), y=map(operator.itemgetter(1), thispoly))
			cm.plot(rx, ry, 'r-', lw=2)
			plt.fill(rx, ry, fc='r', ec='r', alpha=.4, lw=2)
			
		if rw[2]>1:
			blues+=[[clat, clon, rw[1], rw[2]]]
			thispoly=getSquare([clon, clat], gridsize)
			x=map(operator.itemgetter(0), thispoly)
			y=map(operator.itemgetter(1), thispoly)
			bx, by=cm(x,y)
			cm.plot(bx, by, 'b-', lw=2)
			plt.fill(bx, by, fc='c', ec='b', alpha=.4, lw=2)
	#
	#plt.plot(map(operator.itemgetter(1), reds), map(operator.itemgetter(0), reds), 'ro')
	#plt.plot(map(operator.itemgetter(1), blues), map(operator.itemgetter(0), blues), 'bo')
#	xred, yred=cm(map(operator.itemgetter(1), reds), map(operator.itemgetter(0), reds))
#	xblue, yblue=cm(map(operator.itemgetter(1), blues), map(operator.itemgetter(0), blues))
	#
#	cm.plot(xred, yred, 'yo')
#	cm.plot(xblue, yblue, 'co')
	#plt.show()
	#
	# and where did the earthquakes happen?
	bigeqs=[]
	for rw in thiscat0[mev[4]-1:]:
		if rw[3]>=bigmag:
			# bigeqs+=[rw]
			x,y=cm(rw[2], rw[1])
			cm.plot(x,y, '*', ms=15, label='m%f, %s' % (rw[3], str(rw[0])))
	x,y=cm(mev[2], mev[1])
	cm.plot(x,y, 'r*', ms=25,  label='ms:%f, %s' % (rw[3], str(rw[0])))		
	#
	plt.legend(loc='best', numpoints=1)
	
	
	#
	return catalog


class rbHazMap(object):
	# 
	#plt.ion()
	# note: this hazardmap has probably been bested by hazmaps/hazmap.py
	
	def __init__(self, catalog=None, catnum=0, mc=2.5, gridsize=.5, phi=[0,0], winlen=128, avlen=1, bigmag=5.0, fignum=0):
		# catalog: yp.eqcatalog object, phi: spatial phase shift(s) [lon, lat]
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
		print mev
		self.mev=mev
		self.thiscat=self.thiscat0[0:mev[4]-1]	# eliminate all events after mainshock...
		self.precatindex=len(self.thiscat)-1
		#self.thiscat=self.thiscat0
		#print "mainEvent: %s" % mev
		#
		self.llrange=self.catalog.getLatLonRange(self.thiscat)
		self.llrange[0]=[self.llrange[0][0]+phi[0], self.llrange[0][1]+phi[1]]
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
			#	print "stepping off pre-cat"
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
				#print 'm%f, %s' % (rw[3], str(rw[0]))
		#x,y=self.cm(self.mev[2], self.mev[1])
		#print x,y
		x,y=self.cm(mev[2], mev[1])
		#print x,y
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
		gridStats=[]
		if catalog==None: catalog=self.catalog
		reds=[]
		blues=[]
		#print 'as of %s' % str(self.thiscat[-1][0])
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
			#
			#meanr=sum(rs[-self.avlen:])/float(self.avlen)		# this will give the wrong result; we need to average logs...
			lrs=[]
			for i in xrange(len(rs)):
				lrs+=[math.log10(rs[i])]
			meanr=numpy.average(lrs[-self.avlen:])
			#
			#print meanr
			#
			#gridStats+=[[ct[0], pylab.date2num(ct[1][-1][0])-pylab.date2num(ct[1][0][0]), catalog.rb.getIntervalRatios(mc, winlen, ct[1])[-1][-1]]]
			#gridStats+=[[ct[0], pylab.date2num(ct[1][-1][0])-pylab.date2num(ct[1][0][0]), meanr]]
			rw=[ct[0], pylab.date2num(ct[1][-1][0])-pylab.date2num(ct[1][0][0]), meanr]
			#
			# boxes are already drawn we only have to update the color/alpha
			#clat=self.getLat(rw[0]) #int(rw[0]/self.Nlon)*self.gridsize+self.lat0
			#clon=self.getLon(rw[0]) #int(rw[0]%self.Nlon)*self.gridsize+self.lon0
			#
			if meanr<=1: boxcolor='r'
			if meanr>1: boxcolor='b'
			
			self.gridSquares[catindex].set_alpha(.4)
			self.gridSquares[catindex].set_color(boxcolor)
		#
		self.fig.canvas.draw()
	
			
def getSquare(xy=[], dx=.1):
	# xy = ll corner.
	return [xy, [xy[0], xy[1]+dx], [xy[0]+dx, xy[1]+dx], [xy[0]+dx, xy[1]], xy]	

def getParkfieldQuad(catalog=None,catnum=1, mc=1.5, winlen=256, rbthresh=1.0, bigmag=5.0, fignum=0, intlist=[32, 256, 512], avelen=1, thislw=1.0): 
	#print "winlen: %d" % winlen
	majAx=.4
	majAx=.6
	pftheta=45.0
	
	if catalog==None:
		catalog=yp.eqcatalog([])
		catalog.rb=rbi.intervalRecordBreaker(None)
		#catalog.setCatFromSQL(dtm.datetime(1995,1,1), dtm.datetime.now(), [34.4, 37.4], [-121.11, -119.4], 1.5, "Earthquakes", 523, 'asc')
		catalog.setCatFromSQL(dtm.datetime(1995,1,1), dtm.datetime.now(), [35.4, 36.4], [-121.1, -119.65], 1.5, "Earthquakes", 523, 'asc')
		while len(catalog.subcats)>0: x=catalog.subcats.pop()
		#catalog.addEllipCat('PFshock (.4x.15)', catalog.cat, 40.0, 35.9, -120.5, 0.4, 0.15)
		catalog.addEllipCat('PFshock (.4x.15)', catalog.cat, pftheta, 35.9, -120.5, majAx, 0.15)
		#catalog.addEllipCat('PFshock (.4x.15)', catalog.cat, pftheta, 36.1, -120.7, majAx, 0.15)
		#
	pfcat=rbomoriQuadPlot(catalog, catnum, mc, winlen, rbthresh, bigmag, fignum, intlist, avelen, thislw)
	#
	
	x,y = pfcat.catmap(map(operator.itemgetter(2), pfcat.getcat(0)), map(operator.itemgetter(1), pfcat.getcat(0)))
	pfcat.axlist[3].plot(x, y, 'g,', alpha=.5, zorder=1)
	#pfcat.plotCatMap(catalog=pfcat.getcat(0), legendLoc='best', myaxis=pfcat.axlist[3], doCLF=False)
	#x,y = pfcat.catmap(map(operator.itemgetter(2), pfcat.getcat(1)), map(operator.itemgetter(1), pfcat.getcat(1)))
	#pfcat.axlist[3].plot(x, y, 'b,', alpha=1.0, zorder=2)
	#
	# now, plot an ellipse over the aftershocks...
	# el = Ellipse([-120.5, 35.9], .8, .3, -40, facecolor='b', alpha=0.4)
	rbel=Ellipse([-120.5, 35.9], 2*majAx, .3, -pftheta, facecolor='b', alpha=0.4, lod=True)
	Xel, Yel = pfcat.catmap(rbel.get_verts()[:,0],rbel.get_verts()[:,1])
	pfcat.axlist[3].fill(Xel, Yel, color='b', alpha=.3, zorder=1)
	
	#
	plt.figure(0)
	pfcat.axlist[0].text(.2, .1, 'Magnitude m', rotation=90)
	pfcat.axlist[1].text(.2, .3, 'mean interval taut', rotation=90)
	pfcat.axlist[2].text(.2, .3, 'r(t)', rotation=90)
	

	return pfcat

def rbomoriQuadPlot(catalog=None, catnum=0, mc=2.5, winlen=501, rbthresh=1.0, bigmag=5.0, fignum=0, intlist=[32, 256, 512], rbavelen=1, thislw=1.0):
	# make an awesome quad plot: omor-times, rbRatios, mag-seismicity, map-catalog
	# catalog -> a yodapy.eqcatalog() object
	# thislw: linewidth
	#
	#rbavelen=1	# rb-averaging length (aka, <nrb>_thisnum
	print "winlen: %d" % winlen
	if catalog==None: catalog=getMexicaliCat()
	#intlist=[25, 256, 512]
	
	plt.figure(fignum)
	plt.clf()
	#ax0=plt.axes([.1,.1,.85, .35])
	# define axis (subplots) boundaries:
	ydelim=.03
	xdelim=.05
	xTS=[0.05, .5]
	yTS0=0.05
	dyTS=.3
			
	myaxes=[]
	nax=0
	# mags:
	x0=xTS[0]
	y0=yTS0+dyTS*nax
	#myaxes+=[plt.axes([xTS[0], y0, xTS[1], dyTS])]
	myaxes+=[plt.axes([.1, .03, .45, .3])]
	nax+=1
	# intervals:
	#x0=xTS[0]
	y0=yTS0+dyTS*nax
	#myaxes+=[plt.axes([xTS[0], y0, xTS[1], dyTS], sharex=myaxes[0])]
	myaxes+=[plt.axes([.1, .37, .45, .27], sharex=myaxes[0])]
	nax+=1
	# ratios:
	#x0=xTS[0]
	y0=yTS0+dyTS*nax
	#myaxes+=[plt.axes([xTS[0], y0, xTS[1], dyTS], sharex=myaxes[0])]
	myaxes+=[plt.axes([.1, .68, .45, .3], sharex=myaxes[0])]
	#
	# map:
	nax+=1
	xs=[xTS[1]+xdelim, .95]
	ys=[yTS0, 1.0]
	#myaxes+=[plt.axes([xs[0], xs[1], ys[0], ys[1]])]
	#myaxes+=[plt.axes([.6, .05, .35, .90], sharex=myaxes[0])]
	myaxes+=[plt.axes([.6, .05, .35, .90])]
	#
	# get RB ratios:
	try:
		catalog.rb
	except:
		catalog.rb=rbi.intervalRecordBreaker(None)
	#ratios=self.getIntervalRatios(minmag, windowLen, cat0, deltaipos, avlen)
#	ratios=catalog.rb.getIntervalRatios(mc, winlen, catalog.getcat(catnum), 1)
	#
	#plotIntervals(self, intervals=[10, 100, 1000], minmag=2.0, catalog=None, fignum=0, dtmlatlonmagCols=[0,1,2,3], plotDates=[None, None], thisAxes=None):
	#catalog.plotIntervals(intlist, mc, catalog.getcat(catnum), fignum, [0,1,2,3], [None, None], [myaxes[0], myaxes[1]])
	#format: plotInts(self, intervals=[10, 100, 1000], catalog=None, minmag=2.0, ax=None, dtmlatlonmagCols=[0,1,2,3], plotDates=[None, None]):
	catalog.plotInts(intervals=intlist, catalog=catalog.getcat(catnum), minmag=mc, ax=myaxes[1], thislw=thislw, legendPos='upper left')
	#myaxes[1].set_label('mean intervals $\\tau$')
	myaxes[1].set_ylabel('mean intervals $\\ < \\tau >$', size=14)
	#
	# format: #plotMags(self, catalog=None, minmag=2.0, ax=None, dtmlatlonmagCols=[0,1,2,3], plotDates=[None, None]):	
	catalog.plotMags(catalog.getcat(catnum), mc, myaxes[0])
	myaxes[0].set_ylabel('mag', size=14)
	#
	# plotIntervalRatiosAx(self, minmag=3.0, windowLen=10, cat0=None, hitThreshold=1.0, bigmag=5.0, thisAx=None, ratios=None, deltaipos=1, avlen=1):
	catalog.rb.plotIntervalRatiosAx(mc, winlen, catalog.getcat(catnum), rbthresh, bigmag, myaxes[2], None, 1, rbavelen)
	myaxes[2].set_ylabel('RB ratio $r(N=%d)$' % winlen, size=14)
	
	#plt.figure(fignum)
	myfsize=12
	myaxes[1].text(.1, .1, '$<\\tau>$', None, rotation='vertical', size=myfsize)
	myaxes[0].text(.1, .1, 'mags', None, rotation='vertical', size=myfsize)
	myaxes[2].text(.1, .1, '$r(1) = frac{n_{rb-large}}{n{rb-small}}$', None, rotation='vertical', size=myfsize)
	
	#
	X=catalog.plotCatMap(catalog.getcat(catnum), True, False, None, None, 'best', False, 'b,', myaxes[3], 0, .15)
	# and large events:
	bigEq=[]
	bigeqindex=0
	#print catalog.catmap
	for rw in catalog.getcat(catnum):
		if rw[3]>=bigmag and rw[0]>=catalog.getMainEvent(catalog.getcat(catnum))[0]:
			eqx,eqy = catalog.catmap(rw[2], rw[1])
			#catalog.catmap.plot(eqx, eqy, '*', label='m=%f, %s' % (rw[3], str(rw[0])))
			catalog.catmap.plot(eqx, eqy, '*', label='m=%.2f, %d' % (rw[3], bigeqindex), ms=15)
			bigeqindex+=1
			print "eq: %d, %f, %s" % (bigeqindex, rw[3], str(rw[0]))
	catalog.catmap.ax.legend(loc='best', numpoints=1)
#	#if len(bigEqs)>0:
	#	#XX=catalog.plotCatMap(bigEqs, True, False, None, 'best', False, '*', myaxes[3], fignum)
	#	XX=catalog.plotCatMap(bigEqs, True, False, None, None, 'best', False, '*', myaxes[3], fignum)
	# sharex=ax0
	#
	plt.show()
	
	# problem: the x-axis of RB-ratios is in float form, formatted as datetimes (basically using "plot" type functions). the intervals/magnitudes are date_plot
	# functions. in short, they have different x-axis variable types, so we can't share the rb-ratios x-axis with the other two. it probably makes sense to convert
	# the interval-ratio plots to the rb-ratio style (because the rb-rations uses the "fillbetween" function).
	
	catalog.axlist=myaxes
	return catalog
	

def mexicaliScript(eqm=None):
	eqm=0
	eqm="string"
	eqm=None
	#
	if eqm==None: eqm=getMexicaliCat()
	# addEllipCat(self, subcatname='newcat', fullcat=None, theta=0, clat=35.9, clon=-120.5, ra=1.0, rb=1.0):
	#eqm.addEllipCat('incipElip1', eqm.getcat(0), 60, 32.258, -115.287, .9, .35)
	eqm.addEllipCat('incipElip1', eqm.getcat(0), -60.0, 32.258, -115.287, .9/2.0, .35/2.0)	# note: this is opposite the direction and half the radii of plotting.
	elcatnum=len(eqm.subcats)
	#
	# and try a poly-catalog:
	polyverts=[[-115.4, 31.85], [-115.0, 32.3], [-115.6, 33.5], [-116.2, 32.75], [-115.3, 31.75]]
	eqm.subcats+=[['mexi-poly1', eqm.polycat(eqm.cat, polyverts)]]
	polycatnum=len(eqm.subcats)
	
	winlen=1500
	mc=1.5
	
	rbm=rbi.intervalRecordBreaker(None)
	# plotIntervalRatios(self, minmag=3.0, windowLen=10, cat0=None, hitThreshold=1.0, bigmag=5.0, fignum=0, ratios=None, deltaipos=1, avlen=10)
	#rbm.plotIntervalRatios(1.5, 2500, eqm.cat, 1.0, 5.0, 1)
#	rbm.plotIntervalRatios(1.5, 750, eqm.getcat(1), 1.0, 5.0, 3)
	rbm.plotIntervalRatios(mc, winlen, eqm.getcat(2), 1.0, 5.0, 4)
	
	rbm.plotIntervalRatios(mc, winlen, eqm.getcat(elcatnum), 1.0, 5.0, 5)
	eqm.plotIntervals([100, 250, 500, 750, 1000], mc, eqm.getcat(elcatnum), 6)
	#plotIntervals(self, intervals=[10, 100, 1000], minmag=2.0, catalog=None, fignum=0, dtmlatlonmagCols=[0,1,2,3], plotDates=[None, None]):
	rbm.plotIntervalRatios(mc, winlen, eqm.getcat(polycatnum), 1.0, 5.0, 7)
	eqm.plotIntervals([100, 250, 500, 750, 1000], mc, eqm.getcat(polycatnum), 8)

	# semi-manually get mainshock...	
	#mainshock=eqm.getMainEvent(eqm.getcat(2))
	irw=-1
	for rw in eqm.getcat(0):
		lat=rw[1]
		lon=rw[2]
		mag=rw[3]
		irw+=1
		# easy bits:
		if lat>32.4 or lat<32.1 or lon>-115 or lon<-115.5 or mag<7.0 or mag>7.5: continue
		# check date:
		if rw[0].year==2010 and rw[0].month==4 and rw[0].day==4 and rw[0].hour==22 and rw[0].minute==40 and (rw[0].second>20 and rw[0].second<59):
			# this is good enough, i think, to get the mainshock and allow the catalog to float a bit (aka, as corrected measurements are added).
			mainshock=rw+[irw]
			break
	imev=mainshock[-1]	# index of 
	# plotCatMap(self, catalog=None, doShow=True, doSave=False, saveName='catalogPlot.png', epicenter=None, legendLoc='upper left', doCLF=True, eqicon='b,')
	#
	mapfignum=1
	plt.figure(mapfignum)
	plt.clf()
	cm=eqm.plotCatMap(eqm.getcat(0)[0:imev-1000], mapfignum, False, None, None, 'upper left', True, 'b,')
	cm=eqm.plotCatMap(eqm.getcat(0)[imev:], mapfignum, False, None, None, 'upper left', False, 'g,')
	cm=eqm.plotCatMap(eqm.getcat(0)[imev-1000:imev], mapfignum, False, None, None, 'upper left', False, 'r,')
	#cm=eqm.plotCatMap(eqm.getcat(elcatnum), mapfignum, False, None, None, 'upper left', False, 'c,')
	
	# now, can we draw the ellipse we want to use?
	from matplotlib.patches import Ellipse
	#f=plt.figure(0)
	#ax1=f0.gca()
	el = Ellipse([-115.287, 32.258], .9, .35, 60, facecolor='b', alpha=0.4)
	Xel, Yel = cm(el.get_verts()[:,0],el.get_verts()[:,1])
	cm.plot(Xel, Yel, '-r', lw=2)
	cm.ax.fill(Xel, Yel, ec='r', fc='r', alpha=.4)
	#
	# draw the polygon:
	polyx, polyy = cm(map(operator.itemgetter(0), polyverts), map(operator.itemgetter(1), polyverts))
	cm.plot(polyx, polyy, 'g-')
	
	
	return eqm

#def mexicaliEventsMovie(fnum=0, eqm=None):	
def mexicaliEventsMovie(fnum=0, eqm=None, mapcat=None):
	doSave=False
	if eqm==None: eqm=mexicaliScript()
	if mapcat==None: mapcat=eqm.cat
	print "catlen: %d" % len(mapcat)
	#
	#mapcenter=[-115.287, 32.258]
	f0=plt.figure(fnum)
	plt.clf()
	plt.ion()	# interactive mode on...
	plt.show()
	#
	#set up map:
	dWidth=0.5
	llr=eqm.getLatLonRange(mapcat)	# latLonRange
	llr[0][0]-=dWidth
	llr[0][1]-=dWidth
	llr[1][0]+=dWidth
	llr[1][1]+=dWidth
	#
	bigmag=5.0
	bigmags=[[],[], []]
	for rw in mapcat:
		if rw[3]>=bigmag:
			bigmags[0]+=[rw[2]]
			bigmags[1]+=[rw[1]]
			bigmags[2]+=[rw[3]]
	# semi-manually get main event (we might have  broad catalog that includes a larger earthquake in the area):
	irw=-1
	for rw in mapcat:
		lat=rw[1]
		lon=rw[2]
		mag=rw[3]
		irw+=1
		# easy bits:
		if lat>32.4 or lat<32.1 or lon>-115 or lon<-115.5 or mag<7.0 or mag>7.5: continue
		# check date:
		if rw[0].year==2010 and rw[0].month==4 and rw[0].day==4 and rw[0].hour==22 and rw[0].minute==40 and (rw[0].second>20 and rw[0].second<59):
			# this is good enough, i think, to get the mainshock and allow the catalog to float a bit (aka, as corrected measurements are added).
			mainshock=rw+[irw]
			break
		
	######
	print "mainshock: %s" % (str(mainshock))	
	
	cntr=[float(llr[0][0])+(llr[1][0]-float(llr[0][0]))/2.0, float(llr[0][1])+(llr[1][1]-float(llr[0][1]))/2.0]
	catmap=yp.Basemap(llcrnrlon=llr[0][1], llcrnrlat=llr[0][0], urcrnrlon=llr[1][1], urcrnrlat=llr[1][0], resolution ='l', projection='tmerc', lon_0=cntr[1], lat_0=cntr[0])
	canvas=yp.FigureCanvas(f0)	# this comes from: from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
	catmap.ax=f0.add_axes([0,0,1,1])
	f0.set_figsize_inches((12/catmap.aspect,12.))
	#
	catmap.drawcoastlines(color='gray')
	catmap.drawcountries(color='gray')
	catmap.drawstates(color='gray')
	catmap.fillcontinents(color='beige')
	
	catmap.drawmeridians(range(int(llr[0][1]-2.0), int(llr[1][1]+2.0)), color='k', labels=[1,1,1,1])
	catmap.drawparallels(range(int(llr[0][0]-2.0), int(llr[1][0]+2.0)), color='k', labels=[1, 1, 1, 1])
	
	#xfull, yfull=catmap(map(operator.itemgetter(2), catalog), map(operator.itemgetter(1), catalog))
	#epx, epy=catmap(epicenter[0], epicenter[1])
	epx, epy = catmap(mainshock[2], mainshock[1])
	catmap.plot(epx, epy, 'r*')
	
	#catmap.plot(xfull, yfull, 'b,', label='Full Catalog')
	#catmap.plot(xfull, yfull, eqicon, label='Full Catalog')
	dN=100
	beforeIcon='b,'
	incipIcon='r,'
	afterIcon='g,'
	thisIcon=beforeIcon
	i=dN
	# fetch coords once so we don't mess with that during the graphical bit.
	X=map(operator.itemgetter(2), mapcat)
	Y=map(operator.itemgetter(1),mapcat)
	maxi=len(mapcat)-1
	imev=mainshock[4]
	logimax=int(math.ceil((math.log10(maxi))))
	zeroBaseString="0"
	for ii in xrange(logimax):
		zeroBaseString+="0"
		
	while i<maxi:
		# first, let's not over-flow, and let's demarcate at the event:
		if i>maxi: i=maxi
		if i>imev and i<imev+dN: i=imev
		if i<=imev: thisIcon=beforeIcon
		if i>imev-1000: thisIcon=incipIcon
		if i>imev: thisIcon=afterIcon
		#
		theseX, theseY=catmap(X[i-dN:i-1], Y[i-dN:i-1])
		catmap.plot(theseX, theseY, thisIcon)
		if float(int(i/dN))%100==0: catmap.plot(epx, epy, 'r*')
		#plt.title("date: %s" % str(mapcat[i][0]))
		#
		# make movie frames:
		
		indexString=zeroBaseString + str(i)
		indexString=indexString[-logimax-1:]
		foutname="movieFrames1/movieFrames1-%s.png" % indexString
		canvas.print_figure(foutname)
		os.system("convert %s %s%s" % (foutname, foutname[:-3], 'jpg'))
		####
		#
		i+=dN
	
	theseX, theseY=catmap(X[imev-1000:imev], Y[imev-1000:imev])
	catmap.plot(theseX, theseY, 'r,')
	catmap.plot(epx, epy, 'r*')
	bigMX, bigMY=catmap(bigmags[0], bigmags[1])
	catmap.plot(bigMX, bigMY, 'c*')
	#
	# if we are inclned to save:
	if doSave and saveName!=None: canvas.print_figure(saveName)
	
	return eqm

def mexicaliEventsMovie2(fnum=0, eqm=None):
	# in this case, a running window of N events changing dN at a time...
	# this is taken from a standard pylab plot script from james, but basemap appears to handle
	# some of these plotting featuers differently, so this doesn't really work.
	
	doSave=False
	if eqm==None: eqm=mexicaliScript()
	mapcat=eqm.cat
	#
	#mapcenter=[-115.287, 32.258]
	f0=plt.figure(fnum)
	plt.clf()
	plt.ion()	# interactive mode on...
	plt.show()
	#
	#set up map:
	dWidth=0.5
	llr=eqm.getLatLonRange(eqm.cat)	# latLonRange
	llr[0][0]-=dWidth
	llr[0][1]-=dWidth
	llr[1][0]+=dWidth
	llr[1][1]+=dWidth
	#
	# semi-manually get main event (we might have  broad catalog that includes a larger earthquake in the area):
	irw=-1
	for rw in eqm.cat:
		lat=rw[1]
		lon=rw[2]
		mag=rw[3]
		irw+=1
		# easy bits:
		if lat>32.4 or lat<32.1 or lon>-115 or lon<-115.5 or mag<7.0 or mag>7.5: continue
		# check date:
		if rw[0].year==2010 and rw[0].month==4 and rw[0].day==4 and rw[0].hour==22 and rw[0].minute==40 and (rw[0].second>20 and rw[0].second<59):
			# this is good enough, i think, to get the mainshock and allow the catalog to float a bit (aka, as corrected measurements are added).
			mainshock=rw+[irw]
			break
		
	######
	print "mainshock: %s" % (str(mainshock))	
	
	cntr=[float(llr[0][0])+(llr[1][0]-float(llr[0][0]))/2.0, float(llr[0][1])+(llr[1][1]-float(llr[0][1]))/2.0]
	catmap=yp.Basemap(llcrnrlon=llr[0][1], llcrnrlat=llr[0][0], urcrnrlon=llr[1][1], urcrnrlat=llr[1][0], resolution ='l', projection='tmerc', lon_0=cntr[1], lat_0=cntr[0])
	canvas=yp.FigureCanvas(f0)	# this comes from: from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
	catmap.ax=f0.add_axes([0,0,1,1])
	f0.set_figsize_inches((12/catmap.aspect,12.))
	#
	catmap.drawcoastlines(color='gray')
	catmap.drawcountries(color='gray')
	catmap.drawstates(color='gray')
	catmap.fillcontinents(color='beige')
	
	catmap.drawmeridians(range(int(llr[0][1]-2.0), int(llr[1][1]+2.0)), color='k', labels=[1,1,1,1])
	catmap.drawparallels(range(int(llr[0][0]-2.0), int(llr[1][0]+2.0)), color='k', labels=[1, 1, 1, 1])
	
	#xfull, yfull=catmap(map(operator.itemgetter(2), catalog), map(operator.itemgetter(1), catalog))
	#epx, epy=catmap(epicenter[0], epicenter[1])
	epx, epy = catmap(mainshock[2], mainshock[1], 'ro')
	catmap.plot(epx, epy, 'r*')
	
	#catmap.plot(xfull, yfull, 'b,', label='Full Catalog')
	#catmap.plot(xfull, yfull, eqicon, label='Full Catalog')
	winLen=500
	dN=100
	beforeIcon='b,'
	afterIcon='g,'
	thisIcon=beforeIcon
	i=winLen
	# fetch coords once so we don't mess with that during the graphical bit.
	X=map(operator.itemgetter(2), mapcat)
	Y=map(operator.itemgetter(1),mapcat)
	maxi=len(mapcat)-1
	imev=mainshock[4]
	#line1=catmap.plot([], [], thisIcon) 
	while i<len(mapcat):
		# first, let's not over-flow, and let's demarcate at the event:
		if i>maxi: i=maxi
		if i>imev and i<imev+dN: i=imev
		if i<=imev: thisIcon=beforeIcon
		if i>imev: thisIcon=afterIcon
		#
		#theseX, theseY=catmap(X[i-dN:i-1], Y[i-dN:i-1])
		theseX, theseY=catmap(X[i-winLen:i-1], Y[i-winLen:i-1])
		#plt.clf()
		#canvas.draw()
		catmap.plot(theseX, theseY, thisIcon)
		#line1.set_data(theseX, theseY, thisIcon)
		canvas.draw()
		#
		i+=dN
		
	# if we are inclned to save:
	if doSave and saveName!=None: canvas.print_figure(saveName)
	#
	return eqm
	

def getMexicaliCatFromSQL(outFile='mexicali.cat', startDate=dtm.datetime(1990, 01, 01), endDate=dtm.datetime.now(), minmag=2.0, lats=[31.0, 35.0], lons=[-118.0, -113.5], catID=523 ):
	yp.updateANSS2SQL()
	#
	sqlSelect = "select eventDateTime, lat, lon, mag from Earthquakes where catalogID=%d and eventDateTime between '%s' and '%s' and lat between %f and %f and lon between %f and %f and mag>=%f order by eventDateTime asc" % (catID, str(startDate), str(endDate), lats[0], lats[1], lons[0], lons[1], minmag)
	sqlHost = 'localhost'
	sqlUser = 'myoder'
	sqlPassword = 'yoda'
	sqlPort = 3306
	sqlDB = 'QuakeData'
	myConn = MySQLdb.connect(host=sqlHost, user=sqlUser, passwd=sqlPassword, port=sqlPort, db=sqlDB)
	r1=myConn.cursor()
	r1.execute(sqlSelect)
	#	
	fout=open(outFile, 'w')
	fout.write("# chile catalog from hatistuff.py\n#prams: outfile, startDate, endDate, minmag, lats, lons, catID\n")
	fout.write("#%s\t%s\t%s\t%s\t%s\t%s\t%s\d\n" % (outFile, startDate, endDate, minmag, lats, lons, catID))
	# catalog format (tab delimited):
	# 2004/09/27	01:53:42.78	35.5277	-120.8433	1.42
	catlen=0
	for rw in r1:
		fout.write("%d/%d/%d\t%d:%d:%d.%d\t%f\t%f\t%f\n" % (rw[0].year, rw[0].month, rw[0].day, rw[0].hour, rw[0].minute, rw[0].second, rw[0].microsecond, rw[1], rw[2], rw[3]))
		catlen+=1
	fout.close()
	r1.close()
	myConn.close()
	
	return catlen

def getStandardMexicaliRB(hlat=32.258, hlon=-115.287, Theta=50.0, ra=2.0, rb=.5, thisCat='mexicali.cat'):
	#reload(rbi)
	#thisrb=rbi.intervalRecordBreaker(thisCat, Theta, hlat, hlon, ra, rb, dtm.datetime(2010,3,3), dtm.datetime.now(), 0)
	#thisrb.setAftershockCatalog(thisCat, Theta, hlat, hlon, ra, rb, dtm.datetime(2010,3,3), dtm.datetime.now(), 0)
	
	#thisrb=rbi.intervalRecordBreaker(thisCat, Theta, hlat, hlon, ra, rb)
	#thisrb.setAftershockCatalog(thisCat, Theta, hlat, hlon, ra, rb, None, dtm.datetime.now(), 0)
	
	thisrb=rbi.intervalRecordBreaker(thisCat, Theta, hlat, hlon, ra, rb, dtm.datetime(2010, 4, 4, 22, 40, 42))
	thisrb.setAftershockCatalog(thisCat, Theta, hlat, hlon, ra, rb, dtm.datetime(2010, 4, 4, 22, 40, 42), dtm.datetime.now(), 0)
	
	return thisrb

def getStandardTaiwanRB(hlat=23.0, hlon=121.0, Theta=50.0, ra=2.0, rb=.5, thisCat='cats/taiwan.cat'):
#def getStandardTaiwanRB(hlat=23.0, hlon=121.0, Theta=50.0, ra=2.0, rb=.5, thisCat='shinesStuff/chighi2.dat'):
	# 'shinesStuff/chichiareadata.txt'
	thisrb=rbi.intervalRecordBreaker(thisCat, Theta, hlat, hlon, ra, rb, dtm.datetime(1998, 1, 1, 1, 1, 1))
	thisrb.setAftershockCatalog(thisCat, Theta, hlat, hlon, ra, rb, dtm.datetime(1998, 1, 1, 1, 1, 1), dtm.datetime.now(), 0)
	return thisrb

def freshMexicaliRatios(wLens=[15, 20, 30, 40, 50], bigShock=6.0, haitiCat='mexicali.cat', doShow=False, minMag=5.0):
	# full set of Haiti rb ratios from a fresh catalog...
	# update catalog:
	print "update MySQL catalog..."
	yp.updateANSS2SQL()
	print "update local catalog, %s" % haitiCat
	#getChileCatFromSQL(haitiCat)
	#
	return mexicaliRatios(wLens, bigShock, haitiCat, doShow, minMag)
	
#def mexicaliRatios(wLens=[25, 50, 75, 100, 250, 500], bigShock=5.0, cat='mexicali.cat', doShow=False, minMag=2.5, dosets=[1,0]):
def mexicaliRatios(wLens=[75, 100, 250, 500], bigShock=4.5, cat='mexicali.cat', doShow=False, minMag=4.0, dosets=[1,0]):
	# full set of mexicali rb ratios
	# plot GR distributions:
	rbMexical=getStandardMexicaliRB()
	rbMexical.GRshock(False, fname='images/mexicaliShock-GRdist.png')
	rbMexical.GRfullcat(False, fname='images/mexicaliFull-GRdist.png')
	rbMexical=None
	
	# calc ratios for wLens:
	print "calc ratios for wLens"
	for wLen in wLens:
		print "rb ratios for wLen=%d" % wLen
		getMexicaliRB('Mexicali', None, wLen, bigShock, cat, doShow, dtm.datetime.now(), minMag, dosets)
	#
	return 0

def convertSunshinesOtherCat():
	fin=open('shinesStuff/chichiareadata_2.dat')
	fout=open('shinesStuff/chighi2.dat', 'w')
	
	lcat=[]	# list version of catalog.
	
	for rw in fin:
		if rw[0]=='#': continue
		rws=rw.split('\t')
		yr=int(float(rws[1]))
		mnth=int(float(rws[2]))
		dy=int(float(rws[3]))
		hrs=int(float(rws[4]))
		mn=int(float(rws[5]))
		fsc=float(float(rws[6]))
		isc=int(fsc)
		msc=int((10**6)*fsc%1)
		lon=float(rws[7])
		lat=float(rws[8])
		mag=float(rws[9])
		#
		#thisdtm=dtm.datetime(yr, mnth, dy, hrs, mn, isc, msc)
		#
		#fout.write('%d/%d/%d\t%d:%d:%d.%d\t%f\t%f\t%f\n' %(yr, mnth, dy, hrs, mn, isc, msc, lat, lon, mag))
		fout.write('%d/%d/%d\t%d:%d:%f\t%f\t%f\t%f\n' %(yr, mnth, dy, hrs, mn, fsc, lat, lon, mag))
		evdtm=dtm.datetime(yr, mnth, dy, hrs, 0, 0, 0)
		evdtm=evdtm+dtm.timedelta(minutes=mn)
		evdtm=evdtm+dtm.timedelta(seconds=isc)
		evdtm=evdtm+dtm.timedelta(microseconds=msc)	# this dumb way of putting together a datetime object should circumvent the "secs between 0,59" problem.
		#
		lcat+=[[evdtm, lat, lon, mag]]
		
	#
	fin.close()
	fout.close()
	
	# and sort by date:
	lcat.sort(key=lambda rw: rw[0])	# by default, sort() should sort on first element, but use full syntax to be sure...
	
	return lcat

def getCatFromShineFile(fname='shinesStuff/chichiareadata_2.dat'):
	# get an eqcatalog() object from sunshines catalog. this is not the standard catalog format that i usually use.
	lcat=[]	# list version of catalog.
	fin=open(fname)
	for rw in fin:
		if rw[0]=='#': continue
		rws=rw.split('\t')
		yr=int(float(rws[1]))
		mnth=int(float(rws[2]))
		dy=int(float(rws[3]))
		hrs=int(float(rws[4]))
		mn=int(float(rws[5]))
		fsc=float(float(rws[6]))
		isc=int(fsc)
		msc=int((10**6)*fsc%1)
		lon=float(rws[7])
		lat=float(rws[8])
		mag=float(rws[9])
		#
		#thisdtm=dtm.datetime(yr, mnth, dy, hrs, mn, isc, msc)
		#
		#fout.write('%d/%d/%d\t%d:%d:%d.%d\t%f\t%f\t%f\n' %(yr, mnth, dy, hrs, mn, isc, msc, lat, lon, mag))
		#fout.write('%d/%d/%d\t%d:%d:%f\t%f\t%f\t%f\n' %(yr, mnth, dy, hrs, mn, fsc, lat, lon, mag))
		evdtm=dtm.datetime(yr, mnth, dy, hrs, 0, 0, 0)
		evdtm=evdtm+dtm.timedelta(minutes=mn)
		evdtm=evdtm+dtm.timedelta(seconds=isc)
		evdtm=evdtm+dtm.timedelta(microseconds=msc)	# this dumb way of putting together a datetime object should circumvent the "secs between 0,59" problem.
		#
		lcat+=[[evdtm, lat, lon, mag]]
		
	#
	lcat.sort(key=lambda rw: rw[0])	# by default, sort() should sort on first element, but use full syntax to be sure...
	fin.close()
	return yp.eqcatalog(lcat)


def convertSunshinesCatalog():
	f=open('cats/cwbeqp.txt')
	#f=open('cats/testcat.txt')
	fout=open('cats/taiwan.cat', 'w')
	#
	fout.write("#sunshine's taiwan catalog, converted to RB format\n")
	fout.write("#output:\n")
	# output format:
	fout.write("#YYY/MM/DD \\t HH:mm:SS.ss \\t lat \\t lon \\t mag \\t depth (optional)\n")
	for rw in f:
		#rw=thisrw
		if rw[0]=='#': continue
		#print "rw1: %s" % rw
		while "  " in rw: rw=rw.replace("  ", " ")
		while "\n" in rw: rw=rw.replace("\n", "")
		while "\r" in rw: rw=rw.replace("\r", "")
		rw=rw.replace(" ", "\t")
		#print "rw2: %s" % rw
		rws=rw.split('\t')
		#return rws
		lon=rws[0]
		lat=rws[1]
		strDt='%s/%s/%s' % (rws[2], rws[3], rws[4])
		strTm='%s:%s:%s' % (rws[5], rws[6], rws[7])
		mag=rws[9]
		depth=rws[8]
		#
		fout.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (strDt, strTm, lat, lon, mag, depth))
		#
	#
	f.close()
	fout.close()

def sunshine2SQL(incat='cats/cwbeqp.txt'):
	#f=open('cats/cwbeqp.txt')
	#f=open('cats/testcat.txt')
	
	#import MySQLdb
	sqlHost = 'localhost'
	sqlUser = 'myoder'
	sqlPassword = 'yoda'
	sqlPort = 3306
	sqlDB = 'QuakeData'
	catID=21	# no reason...
	myConn = yp.MySQLdb.connect(host=sqlHost, user=sqlUser, passwd=sqlPassword, port=sqlPort, db=sqlDB)
	# myConn.query(strDel)
	
	f=open(incat)
	#fout=open('cats/taiwan.cat', 'w')
	#
	for rw in f:
		#rw=thisrw
		if rw[0]=='#': continue
		#print "rw1: %s" % rw
		while "  " in rw: rw=rw.replace("  ", " ")
		while "\n" in rw: rw=rw.replace("\n", "")
		while "\r" in rw: rw=rw.replace("\r", "")
		rw=rw.replace(" ", "\t")
		#print "rw2: %s" % rw
		rws=rw.split('\t')
		#return rws
		lon=rws[0]
		lat=rws[1]
		strDt='%s/%s/%s' % (rws[2], rws[3], rws[4])
		strTm='%s:%s:%s' % (rws[5], rws[6], rws[7])
		mag=rws[9]
		depth=rws[8]
		#
		#fout.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (strDt, strTm, lat, lon, mag, depth))
		strIns="insert into Earthquakes (catalogID, eventDateTime, lat, lon, mag, depth) values (%d, \'%s %s\', %s, %s, %s, %s)" % (catID, strDt.replace('/', '-'), strTm, lat, lon, mag, depth)
		myConn.query(strIns)
		#return strIns
		#
	#
	f.close()
	myConn.close()
	#fout.close()

def getMexicaliCatPlot(doShow=False, doSave=True, objRB=None, saveName='images/mexicaliCatalogPlot.png'):
	if objRB==None: objRB=getStandardMexicaliRB()
	return getRBCatPlot(doShow, doSave, objRB, saveName)

def getMexicaliRB(eventName='Mexicali', objRB=None, wLen=250, bigShock=5.0, cat='mexicali.cat', doShow=True, maxDt=None, minMag=2.5, dosets=[1,0], startDate=None):
	obj=getStandardMexicaliRB()
	return getRBratios(eventName, obj, wLen, bigShock, cat, doShow, maxDt, minMag, dosets, startDate)

################################
def getTaiwanRB(eventName='Taiwan', objRB=None, wLen=250, bigShock=5.0, cat='mexicali.cat', doShow=True, maxDt=None, minMag=2.5, dosets=[1,0], startDate=None):
	obj=getStandardTaiwanRB()
	return getRBratios(eventName, obj, wLen, bigShock, cat, doShow, maxDt, minMag, dosets, startDate)

def taiwanRBscript2(winlen=1000, mc=2.0, stepsize=1, cat1=None, objrb=None):
	if cat1==None: cat1=getCatFromShineFile()	# get a catalog from sunshine's special catalog
	if objrb==None: objrb=getStandardTaiwanRB()	# we don't really use the catalog that this loads, but this way we know we load a clean RB object.
	#
	# plotIntervalRatios(self, minmag=3.0, windowLen=10, cat0=None, hitThreshold=1.0, bigmag=5.0, fignum=0, ratios=None, deltaipos=1):
	X1=cat1.plotIntervals([100, 250, 500, 1000, 2500], 2.0)
	X2=objrb.plotIntervalRatios(mc, winlen, cat1.cat, 1.0, 5.5, 2, None, stepsize)
	#
	return [cat1, objrb]
	
def taiwanRBscript(winlen=1000, cat1=None, objrb=None):
	if objrb==None: objrb=getStandardTaiwanRB()
	if cat1==None: 
		# original taiwan catalog:
		cat1=yp.eqcatalog()
		print "set taiwan catalog."
		cat1.setSpecialCatSQL('taiwan')
		#c1.addLatLonSubcat(subcatname='xysubcat-120-122-23-25', c1.cat, lats=[23,25], lons=[120, 122], llcols=[1,2]):
		# or use:
		#tcat1=getCatFromShineFile()	# which uses a special catalog provided by shine herself.
#	print "add xy subcat."
	
	#cat1.subcats+=[['xysubcat-120-122-23-25', cat1.getCatFromSQL(dtm.datetime(1995, 1,1), dtm.datetime(2005,1,1), [23,25], [120,121.5], 2.5, 'Earthquakes', 21)]]
	#cat1.addxytmSubcat('xytm-120-122-23-25-m25', cat1.cat, [dtm.datetime(1995,1,1), dtm.datetime(2004,12,31)], [23, 25], [120, 11.5], 2.5, [1,2,3])
#	subcat=cat1.getxytmSubcat(cat1.cat, [dtm.datetime(1990,1,1), dtm.datetime(2004,12,31)], [23, 25], [120, 121.5], 2.5, [1,2,3])
#	cat1.subcats+=[['xytmsubcat', subcat]]
	# addEllipCat(self, subcatname='newcat', fullcat=None, theta=0, clat=35.9, clon=-120.5, ra=1.0, rb=1.0):
#	cat1.addEllipCat('ellipcat1', subcat, 60, 23.85, 120.82, 1.5, .75)
	#return subcat
	#subcat=cat1.getCatFromSQL(dtm.datetime(1995, 1,1), dtm.datetime(2005,1,1), [23,25], [120,121.5], 2.5, 'Earthquakes', 21)
	#
	print "plot interval ratios"
	# plotIntervalRatios(self, minmag=3.0, windowLen=10, cat0=None, hitThreshold=1.0, bigmag=5.0, fignum=0, ratios=None, deltaipos=1):
	#A=objrb.plotIntervalRatios(2.5, winlen, subcat, 1.0, 5.5, 0)
	#A=objrb.plotIntervalRatios(2.5, winlen, objrb.fullCat, 1.0, 5.5, 0)
	A=objrb.plotIntervalRatios(2.5, winlen, cat1.cat, 1.0, 5.5, 0)
	#
	return [cat1, objrb]

###
# generally applicable tools...
def getRBCatPlot(doShow=False, doSave=True, objRB=None, saveName='images/chileCatalogPlot.png'):
	# objRB is a Chile specific (derived from haiti specific) intervalRecordBreaker() object.
	# note: this functionality has been integrated into the intervalRecordBreaker() object as obh.xyPlotCatalogs(doShow, doSave, saveName, ellipseCenter=[x,y])
	if objRB==None: objRB=getStandardchileRB()
	mainEvent=objRB.getMainEvent()
	
	bangDate=dtm.datetime(2010, 2, 27, 6, 34, 14)
	
	fcat=[]
	scat=[]
	aftershocks=[]
	for rw in objRB.shockCat:
		scat+=[rw[0:4]]
	for rw in objRB.fullCat:
		if rw not in scat: fcat+=[rw]
		if rw[0]>=bangDate: aftershocks+=[rw]
	
	#print "aftershocks: %d" % len(aftershocks)
	#return [scat, fcat]
	
	plt.figure(0)	
	plt.clf()
	#
	ax=plt.gca()
	
	tLat=None # 35.9
	tLon=None # -120.5
	tTheta=None # 40.0		#47?	note: tTheta is the angle CCW of the x' (transformed) axis from the x axis.
	tA=None # .4		# ellipse axes
	tB=None # .15
	
	#el = Ellipse((-72.533,18.457), 1.25, .25, 15, facecolor='r', alpha=0.5)
	#el = Ellipse((-72.719,-35.846), 2.0*objRB.tA, 2.0*objRB.tB, objRB.tTheta, facecolor='b', alpha=0.4)
	el = Ellipse((objRB.tLon, objRB.tLat), 2.0*objRB.tA, 2.0*objRB.tB, -objRB.tTheta, facecolor='b', alpha=0.3)
	ax.add_artist(el)
	#
	#plt.plot(map(operator.itemgetter(2), objRB.fullCat), map(operator.itemgetter(1), objRB.fullCat), '.')
	plt.plot(map(operator.itemgetter(2), fcat), map(operator.itemgetter(1), fcat), 'b.')
	plt.plot(map(operator.itemgetter(2), scat), map(operator.itemgetter(1), scat), 'g.')
	#plt.plot(map(operator.itemgetter(2), objRB.shockCat), map(operator.itemgetter(1), objRB.shockCat), '.')
	plt.plot(map(operator.itemgetter(2), aftershocks), map(operator.itemgetter(1), aftershocks), 'r.')
	
	#plt.plot(map(operator.itemgetter(2), fcat), map(operator.itemgetter(1), fcat), '+', label='Full Catalog')
	#plt.plot(map(operator.itemgetter(2), scat), map(operator.itemgetter(1), scat), '.', label='Aftershock zone')
	plt.plot([objRB.tLon], [objRB.tLat], 'r*', ms=15, label='M%s epicenter' % str(mainEvent[3]))
	plt.legend(loc='upper left', numpoints=1)
	if doSave: plt.savefig(saveName)
	if doShow: plt.show()
	
	return objRB

def getRBratios(eventName='Chile', objRB=None, wLen=25, bigShock=5.2, cat='chile.cat', doShow=True, maxDt=None, minMag=5.0, dosets=[1,1], startDate=None):
	# plotting note: #lats=map(operator.itemgetter(1), rbi.shockCat)
	# catFname='parkcat.cat', theta=tTheta, clat=tLat, clon=tLon, ra=tA, rb=tB
	# haiti epicenter bits (so far, very approximate):
	doBoxy=dosets[1]
	doShock=dosets[0]
	
	avlen=1
	hitThreshold=1.5
	legLoc='lower left'
	
	if cat==None: cat='chile.cat'
	reload(rbi)
	if objRB==None: objRB=getStandardchileRB()
	
	#bigShock=5.7	# magnitude of a big aftershock.
	bigShockShocks=[]
	bigFullShocks=[]
	#
	#
	hlat=objRB.tLat
	#hlon=-72.53	# actual epicenter
	hlon=objRB.tLon
	chileTheta=objRB.tTheta
	ra=objRB.tA
	rb=objRB.tB
	print "prams: %f, %f, %f, %f, %f" % (hlat, hlon, chileTheta, ra, rb)
	# hlat=-35.846, hlon=-72.719, chileTheta=-55, ra=4.5, rb=3.75
	objRB.setAftershockCatalog(cat, chileTheta, hlat, hlon, ra, rb, startDate, dtm.datetime.now(), 0)
	#
	objRB.GRshock(False, fname='images/%sShock-GRdist.png' % (eventName))
	objRB.GRfullcat(False, fname='images/%sFull-GRdist.png' % eventName)
	#
	# get lists of bigshocks:
	rnum=0
	for rw in objRB.shockCat:
		if rw[3]<minMag: continue
		if rw[3]>=bigShock: bigShockShocks+=[[rnum] + rw]
		rnum+=1
	rnum=0
	for rw in objRB.fullCat:
		if rw[3]<minMag: continue
		if rw[3]>=bigShock: bigFullShocks+=[[rnum] + rw]
		rnum+=1
	#
	#objRBSquare=rbi.intervalRecordBreaker(cat, chileTheta, hlat, hlon, ra, rb)	# a square catalog. we'll probably just impose a square catalog...
	# find mainshock(s):
	shockMainshockEvNum=0
	fullMainshockEvNum=0
	maxMag1=objRB.shockCat[0][3]
	maxMag2=objRB.fullCat[0][3]
	rwnum1=0
	rwnum2=0
	for rw in objRB.shockCat:
		if rw[3]<minMag: continue
		if rw[3]>maxMag1:
			maxMag1=rw[3]
			shockMainshockEvNum=rwnum1
			eventFloatDt=yp.datetimeToFloat(rw[0])
		rwnum1+=1
	#print "event date: %f" % (eventFloatDt)
	for rw in objRB.fullCat:
		if rw[3]<minMag: continue
		if rw[3]>maxMag2:
			maxMag2=rw[3]
			fullMainshockEvNum=rwnum2
		rwnum2+=1
	print "mainshock data: mag=%f/%f, evNum1=%d, evNum2=%d" % (maxMag1, maxMag2, shockMainshockEvNum, fullMainshockEvNum)
	#
	# now, get ratios:
	#wLen=10
	
	if doShock:
		print "get shock-zone ratios:"
		hratios=objRB.getIntervalRatios(minMag, wLen, objRB.shockCat)
		if maxDt!=None:
			if maxDt>hratios[-1][1]: hratios+=[[hratios[-1][0]+1, maxDt, hratios[-1][2]]]
		fdts=[]
		for rw in hratios:
			 fdts+=[rw[1].toordinal() + float(rw[1].hour)/24 + rw[1].minute/(24*60) + rw[1].second/(24*3600) + rw[1].microsecond/(24*3600000000)]
	
	if doBoxy:
		print "get boxy-ratios:"
		hratiosBoxy=objRB.getIntervalRatios(minMag,wLen,objRB.fullCat)
		if maxDt!=None:
			if maxDt>hratiosBoxy[-1][1]: hratiosBoxy+=[[hratiosBoxy[-1][0]+1, maxDt, hratiosBoxy[-1][2]]]
		fdtsBoxy=[]
		for rw in hratiosBoxy:
			 fdtsBoxy+=[rw[1].toordinal() + float(rw[1].hour)/24 + rw[1].minute/(24*60) + rw[1].second/(24*3600) + rw[1].microsecond/(24*3600000000)]
	
	print "ratios calculated. do plots."
	
	theseFigs=[]
	
	if doShock:
		hitThreshold=1.4
		f=plt.figure(0)
		theseFigs+=[f]
		plt.clf()
		#plt.semilogy(fdts, map(operator.itemgetter(2), hratios), 'k.')
		#plt.fill_between(fdts, map(operator.itemgetter(2), hratios), y2=1, color='b')
	
		#plt.fill([fdts[0]]+ fdts + [fdts[-1], fdts[0]], [1] + yp.getValsAbove(map(operator.itemgetter(2), hratios), 1) + [1, 1], 'b')
		#plt.fill([fdts[0]]+ fdts + [fdts[-1], fdts[0]], [1] + yp.getValsBelow(map(operator.itemgetter(2), hratios), 1) + [1, 1], 'r')
	
		plt.axvline(x=eventFloatDt, color='c', lw=3, label='mainshock' )
		#plt.fill_between(fdts, scipy.ones(len(fdts),int), map(operator.itemgetter(2), hratios), color='b', where=scipy.array([val>=1 for val in map(operator.itemgetter(2), hratios)]))
		#plt.fill_between(fdts, scipy.ones(len(fdts),int), map(operator.itemgetter(2), hratios), color='r', where=scipy.array([val<=1 for val in map(operator.itemgetter(2), hratios)]))
		reload(yp)
		
		ploty=yp.averageOver(map(operator.itemgetter(2), hratios), avlen)
		plt.fill_between(fdts, hitThreshold*scipy.ones(len(fdts),int), ploty[0], color='b', where=scipy.array([val>=hitThreshold for val in ploty[0]]))
		plt.fill_between(fdts, hitThreshold*scipy.ones(len(fdts),int), ploty[0], color='r', where=scipy.array([val<=hitThreshold for val in ploty[0]]))
		#plt.fill_between(fdts, scipy.ones(len(fdts),int), ploty[0], color='b', where=scipy.array([val>=1 for val in ploty[0]]))
		#plt.fill_between(fdts, scipy.ones(len(fdts),int), ploty[0], color='r', where=scipy.array([val<=1 for val in ploty[0]]))
		#
		#xvals=fdts[yp.greaterof(0, len(fdts)-wLen):]
		#plt.fill_between(xvals, scipy.ones(len(xvals),int), ploty[0][yp.greaterof(0, len(fdts)-wLen):], color='b', where=scipy.array([val>=1 for val in ploty[0][yp.greaterof(0, len(fdts)-wLen):]]))
		#plt.fill_between(xvals, scipy.ones(len(xvals),int), ploty[0][yp.greaterof(0, len(fdts)-wLen):], color='r', where=scipy.array([val<=1 for val in ploty[0][yp.greaterof(0, len(fdts)-wLen):]]))
	
		# note: we don't set the y-log scale in these "fill()" commands. we can do that with axis.set_yscale('log') i think.
		# we achieve this by doing semilogy() plots below.
		plt.title("%s rupture area, time-time, wLen=%d" % (eventName, wLen))
		plt.xlabel('time')
		plt.ylabel('$r=N_{rb-long} / N_{rb-short}$')
		plt.axvline(x=eventFloatDt)
		nbigshocks=0
		for rw in bigShockShocks:
			if nbigshocks==0:
				plt.axvline(x=yp.datetimeToFloat(rw[1]), color='g', label='m > %f' % bigShock)
				nbigshocks+=1
			else:
				plt.axvline(x=yp.datetimeToFloat(rw[1]), color='g')
				plt.plot([yp.datetimeToFloat(rw[1])], rw[4], '*')
			
		plt.semilogy([eventFloatDt], [1], 'r^', ms=10)
		#plt.axvline(x=eventFloatDt, color='r', lw=3, label='mainshock' )
		#	
		#plt.semilogy([y.datetimeToFloat(objRB.shockCat[25][0])], [1], 'o')
		#plt.axvline(x=yp.datetimeToFloat(objRB.shockCat[25][0]))
		plt.axhline(y=1, color='k')
		plt.legend(loc=legLoc, numpoints=2)
		ax=plt.gca()
		fg=plt.gcf()
		ax.xaxis.set_major_formatter(dates.DateFormatter('%Y-%b-%d'))
		ax.set_ylim([.1,10])
		fg.autofmt_xdate()
		plt.savefig('images/%sRuptureTimeTime-Wlen%d-mc%d.png' % (eventName, wLen, int(10*minMag)))
		#

		f=plt.figure(1)
		theseFigs+=[f]
		plt.clf()
		
		#plt.semilogy(map(operator.itemgetter(0), hratios), map(operator.itemgetter(2), hratios), 'k.')
		X=map(operator.itemgetter(0), hratios)
		#plt.fill([X[0]] + X+[X[-1], X[0]], [1] + yp.getValsAbove(map(operator.itemgetter(2), hratios), 1) + [1,1], 'b')
		#plt.fill([X[0]] + X+[X[-1], X[0]], [1] + yp.getValsBelow(map(operator.itemgetter(2), hratios), 1) + [1,1], 'r')
		#plt.bar(X, map(operator.itemgetter(2), hratios), width=1, bottom=1, color='b', align='edge')
		#plt.bar([X[0]] + X+[X[-1], X[0]], [1] + yp.getValsBelow(map(operator.itemgetter(2), hratios), 1) + [1,1], width=1, color='r', orientation='vertical', align='edge')
		#		
		#plt.fill_between(X, scipy.ones(len(X),int), map(operator.itemgetter(2), hratios), color='b', where=scipy.array([val>=1 for val in map(operator.itemgetter(2), hratios)]))
		#plt.fill_between(X, scipy.ones(len(X),int), map(operator.itemgetter(2), hratios), color='r', where=scipy.array([val<1 for val in map(operator.itemgetter(2), hratios)]))
		#
		ploty=yp.averageOver(map(operator.itemgetter(2), hratios), avlen)
		plt.axvline(x=shockMainshockEvNum, color='c', lw=3, label='mainshock')
		plt.axhline(y=hitThreshold, color='k')
		plt.fill_between(X, hitThreshold*scipy.ones(len(X),int), ploty[0], color='b', where=scipy.array([val>=hitThreshold for val in ploty[0]]))
		plt.fill_between(X, hitThreshold*scipy.ones(len(X),int), ploty[0], color='r', where=scipy.array([val<hitThreshold for val in ploty[0]]))
		
		#plt.semilogy([25], [1], 'o')
		#plt.axvline(x=25)
		nbigshocks=0
		for rw in bigShockShocks:
			if nbigshocks==0:
				plt.axvline(x=rw[0], color='g', label='m > %f' % bigShock)
				nbigshocks+=1
			else:
				plt.axvline(x=rw[0], color='g')
				plt.plot([rw[0]], rw[4], '*')
				
		plt.semilogy([shockMainshockEvNum], [1], 'r^', ms=10)
		#plt.axvline(x=shockMainshockEvNum, color='r', lw=3, label='mainshock')
		
		plt.title("%s rupture area, natural-time, wLen=%d" % (eventName, wLen))
		plt.xlabel("Number of Events, n")
		plt.ylabel('$r=N_{rb-long} / N_{rb-short}$')
		# don't over-crowd labels...
		xtks0=map(operator.itemgetter(0), hratios)
		xlbls0=map(operator.itemgetter(1), hratios)
		nskip=len(xtks0)/9	#15 labels seem to fit pretty well.
		if nskip<1:nskip=1
		xlbls=[]
		for i in xrange(len(xlbls0)):
			lbl=''
			if i%nskip==0: lbl=str(xlbls0[i]) + '(%d)' % i
			xlbls+=[lbl]
		plt.xticks(xtks0, xlbls)
		ax=plt.gca()
		ax.set_ylim([.1,10])
		fg=plt.gcf()
		#ax.xaxis.set_major_formatter(dates.DateFormatter('%Y-%b-%d'))
		fg.autofmt_xdate()
		plt.legend(loc=legLoc, numpoints=2)
		plt.savefig('images/rb%sRuptureNaturalTime-Wlen%d-mc%d.png' % (eventName, wLen, int(10*minMag)))		


	if doBoxy:
		f=plt.figure(2)
		theseFigs+=[f]
		plt.clf()
		#plt.semilogy(fdtsBoxy, map(operator.itemgetter(2), hratiosBoxy), '.-')
		boxyShocks=[[],[]]
		for i in xrange(len(fdtsBoxy)):
			if fdtsBoxy[i]>=eventFloatDt:
				boxyShocks[0]+=[fdtsBoxy[i]]
				boxyShocks[1]+=[hratiosBoxy[i][2]]
			
		#plt.fill([fdtsBoxy[0]] + fdtsBoxy + [fdts[-1], fdtsBoxy[0]], [1]+ yp.getValsAbove(map(operator.itemgetter(2), hratiosBoxy), 1) + [1,1], 'b')
		#plt.fill([fdtsBoxy[0]] + fdtsBoxy + [fdts[-1], fdtsBoxy[0]], [1]+ yp.getValsBelow(map(operator.itemgetter(2), hratiosBoxy), 1) + [1,1], 'r')
	
		#plt.fill([boxyShocks[0][0]] + boxyShocks[0] + [boxyShocks[0][-1], boxyShocks[0][0]], [1]+ yp.getValsAbove(boxyShocks[1], 1) + [1,1], 'b')
		#plt.fill([boxyShocks[0][0]] + boxyShocks[0] + [boxyShocks[0][-1], boxyShocks[0][0]], [1]+ yp.getValsBelow(boxyShocks[1], 1) + [1,1], 'r')
	
		plt.fill_between(boxyShocks[0], scipy.ones(len(boxyShocks[0]),int), boxyShocks[1], color='b', where=scipy.array([val>=1 for val in boxyShocks[1]]))
		plt.fill_between(boxyShocks[0], scipy.ones(len(boxyShocks[0]),int), boxyShocks[1], color='r', where=scipy.array([val<=1 for val in boxyShocks[1]]))

	
		#plt.bar([boxyShocks[0][0]] + boxyShocks[0] + [boxyShocks[0][-1], boxyShocks[0][0]], [1]+ yp.getValsAbove(boxyShocks[1], 1) + [1,1], width=1, color='b', orientation='vertical')
		#plt.bar([boxyShocks[0][0]] + boxyShocks[0] + [boxyShocks[0][-1], boxyShocks[0][0]], [1]+ yp.getValsBelow(boxyShocks[1], 1) + [1,1], width=1, color='r', orientation='vertical')
		#
		#print "fdts: %s" % str([fdts[0]] + fdtsBoxy + [fdts[-1], fdts[0]])
		#print "vals: %s" % str([1]+ yp.getValsBelow(map(operator.itemgetter(2), hratiosBoxy), 1) + [1,1])
		#plt.semilogy(fdtsBoxy, yp.getValsBelow(map(operator.itemgetter(2), hratiosBoxy), 1), '.-')
		# mainshock:
	
		#plt.axvline(x=eventFloatDt)
		plt.axvline(x=eventFloatDt, color='c', lw=3, label='mainshock')
		nbigshocks=0
		for rw in bigFullShocks:
			# print "shock events: %s (%f), %f" % (str(rw[1]), yp.datetimeToFloat(rw[1]), eventFloatDt)
			if nbigshocks==0:
				if yp.datetimeToFloat(rw[1])>eventFloatDt:
					plt.axvline(x=yp.datetimeToFloat(rw[1]), color='g', label='m > %f' % bigShock)
					nbigshocks+=1
			else:
				if yp.datetimeToFloat(rw[1])>eventFloatDt: plt.axvline(x=yp.datetimeToFloat(rw[1]), color='g')
		
		plt.semilogy([eventFloatDt], [1], 'r^', ms=10)
		#plt.axvline(x=eventFloatDt, color='r', lw=3, label='mainshock' )
		#nbigshokcs=0	
		#for rw in bigShockShocks:
		#	if nbigshocks==0:
		#		plt.axvline(x=yp.datetimeToFloat(rw[1]), color='g', label='m > %f' % bigShock)
		#	else:
		#		plt.axvline(x=yp.datetimeToFloat(rw[1]), color='g')
		#	nbigshocks+=1
		
		plt.axhline(y=1, color='k')
		plt.title("Box around %s, time-time Aftershocks, wLen=%d" % (eventName, wLen))
		plt.xlabel('time')
		plt.ylabel('$r=N_{rb-long} / N_{rb-short}$')
		plt.legend(loc='lower right', numpoints=2)
		ax=plt.gca()
		fg=plt.gcf()
		ax.xaxis.set_major_formatter(dates.DateFormatter('%Y-%b-%d'))
		ax.set_ylim([.1,10])
		fg.autofmt_xdate()
		#plt.semilogy([eventFloatDt], [1], 'r^')
		plt.savefig('images/rb%sBoxyTimeTimeAftershocks-Wlen%d.png' % (eventName, wLen) )
		#
		# boxy, time-time
		f=plt.figure(3)
		theseFigs+=[f]
		plt.clf()
		#plt.semilogy(fdtsBoxy, map(operator.itemgetter(2), hratiosBoxy), 'k.')
		#plt.fill([fdtsBoxy[0]] + fdtsBoxy + [fdts[-1], fdtsBoxy[0]], [1]+ yp.getValsAbove(map(operator.itemgetter(2), hratiosBoxy), 1) + [1,1], 'b')
		#plt.fill([fdtsBoxy[0]] + fdtsBoxy + [fdts[-1], fdtsBoxy[0]], [1]+ yp.getValsBelow(map(operator.itemgetter(2), hratiosBoxy), 1) + [1,1], 'r')
	
		plt.fill_between(fdtsBoxy, scipy.ones(len(fdtsBoxy),int), map(operator.itemgetter(2), hratiosBoxy), color='b', where=scipy.array([val>=1 for val in map(operator.itemgetter(2), hratiosBoxy)]))
		plt.fill_between(fdtsBoxy, scipy.ones(len(fdtsBoxy),int), map(operator.itemgetter(2), hratiosBoxy), color='r', where=scipy.array([val<=1 for val in map(operator.itemgetter(2), hratiosBoxy)]))

		#print "fdts: %s" % str([fdts[0]] + fdtsBoxy + [fdts[-1], fdts[0]])
		#print "vals: %s" % str([1]+ yp.getValsBelow(map(operator.itemgetter(2), hratiosBoxy), 1) + [1,1])
		#plt.semilogy(fdtsBoxy, yp.getValsBelow(map(operator.itemgetter(2), hratiosBoxy), 1), '.-')
	
		#plt.axvline(x=eventFloatDt)
		plt.axvline(x=eventFloatDt, color='c', lw=3, label='mainshock' )
		nbigshocks=0
		for rw in bigFullShocks:
			if nbigshocks==0:
				plt.axvline(x=yp.datetimeToFloat(rw[1]), color='g', label='m > %f' % bigShock)
				nbigshocks+=1
			else:
				plt.axvline(x=yp.datetimeToFloat(rw[1]), color='g')
	
		plt.semilogy([eventFloatDt], [1], 'r^', ms=10)
		#plt.axvline(x=eventFloatDt, color='r', lw=3, label='mainshock' )
	
		plt.axhline(y=1, color='k')
		plt.title("Box around %s, time-time, wLen=%d" % (eventName, wLen))
		plt.xlabel('time')
		plt.ylabel('$r=N_{rb-long} / N_{rb-short}$')
		plt.legend(loc='upper left', numpoints=2)
		ax=plt.gca()
		fg=plt.gcf()
		ax.xaxis.set_major_formatter(dates.DateFormatter('%Y-%b-%d'))
		ax.set_ylim([.1,10])
		fg.autofmt_xdate()
		#plt.semilogy([eventFloatDt], [1], 'r^')
		plt.savefig('images/rb%sBoxyTimeTime-Wlen%d-mc%d.png' % (eventName, wLen, int(10*minMag)))
	
		# boxy-natural time
		f=plt.figure(4)
		theseFigs+=[f]
		plt.clf()
		plt.axhline(y=1, color='k')
		#
		plt.axvline(x=fullMainshockEvNum, color='c', lw=3, label='mainshock' )
		#
		#plt.semilogy(map(operator.itemgetter(0), hratiosBoxy), map(operator.itemgetter(2), hratiosBoxy), 'k.')
		X=map(operator.itemgetter(0), hratiosBoxy)
		#plt.fill([X[0]] + X + [X[-1], X[0]], [1] + yp.getValsAbove(map(operator.itemgetter(2), hratiosBoxy), 1) + [1,1], 'b')
		#plt.fill([X[0]] + X + [X[-1], X[0]], [1] + yp.getValsBelow(map(operator.itemgetter(2), hratiosBoxy), 1) + [1,1], 'r')
		plt.fill_between(X, scipy.ones(len(X),int), map(operator.itemgetter(2), hratiosBoxy), color='b', where=scipy.array([val>=1 for val in map(operator.itemgetter(2), hratiosBoxy)]))
		plt.fill_between(X, scipy.ones(len(X),int), map(operator.itemgetter(2), hratiosBoxy), color='r', where=scipy.array([val<=1 for val in map(operator.itemgetter(2), hratiosBoxy)]))
		#
		#plt.axvline(x=fullMainshockEvNum)
		nbigshocks=0
		for rw in bigFullShocks:
			if nbigshocks==0:
				plt.axvline(x=rw[0], color='g', label='m > %f' % bigShock)
				nbigshocks+=1
			else:
				plt.axvline(x=rw[0], color='g')
		#
		plt.semilogy([fullMainshockEvNum], [1], 'r^', ms=10)
		#plt.axvline(x=fullMainshockEvNum, color='r', lw=3, label='mainshock' )
		#	
		plt.title("Box around %s, natural-time, wLen=%d" % (eventName, wLen))
		plt.xlabel("Number of Events, n")
		plt.ylabel('$r=N_{rb-long} / N_{rb-short}$')
		plt.legend(loc='upper left', numpoints=2)
		# don't over-crowd labels...
		xtks0=map(operator.itemgetter(0), hratiosBoxy)
		xlbls0=map(operator.itemgetter(1), hratiosBoxy)
		nskip=len(xtks0)/15	#15 labels seem to fit pretty well.
		if nskip<1:nskip=1
		xlbls=[]
		for i in xrange(len(xlbls0)):
			lbl=''
			if i%nskip==0: lbl=str(xlbls0[i]) + '(%d)' % i
			xlbls+=[lbl]
		plt.xticks(xtks0, xlbls)
		ax=plt.gca()
		fg=plt.gcf()
		ax.set_ylim([.1,10])
		#ax.xaxis.set_major_formatter(dates.DateFormatter('%Y-%b-%d'))
		fg.autofmt_xdate()
		plt.axis('auto')
		plt.savefig('images/rb%sBoxyNaturalTime-Wlen%d-mc%d.png' % (eventName, wLen, int(10*minMag)))
	
	
	
	if doShow: plt.show()
	'''
	for F in theseFigs:
		F.clear()
		plt.clf()
		#F.close()
		F=None
	theseFigs=None
	'''
	
	# note: look at the shockCat record-breaking intervals plot; we pretty well forecast the 5.75 aftershock (event 25)
	
	return [objRB, hratios]



