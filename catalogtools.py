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
