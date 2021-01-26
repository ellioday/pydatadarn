#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 01:11:27 2020

@author: elliott
"""

import spacepy.coordinates as coord
from spacepy.time import Ticktock

import datetime as dt
import numpy as np

import aacgmv2
import math

def mag_to_geo(alt, mlat, mlon):
	
	"""
	SPACEPY - converts geomagnetic coordinates to geographic
	
	Parameters
	----------
	
	alt: float
		altitude of coordinates (in km)
		
	mlat: float
		magnetic latitude of geomagnetic coordinates
		
	mlon: float
		magnetic latitude of geomagnetic coordinates
	"""
	Re = 6371
	cvals = coord.Coords([Re+alt, mlat, mlon], "MAG", "sph",["Re", "deg", "deg"])
	#set time epoch for coordinates
	cvals.ticks=Ticktock(["2013-10-02T00:00:00"], "ISO")
	gcoords = cvals.convert("GEO", "sph")
	
	return gcoords
	
def geo_to_mag(alt, lat, lon):
	
	"""
	SPACEPY - converts geographical coordinates to geomagnetic
	
	Parameters
	----------
	
	alt: float
		altitude of coordinates (in km)
		
	lat: float
		latitude of geographical coordinates
		
	lon: float
		longitude of geographical coordinates
	
	"""
	
	Re=6371
	cvals = coord.Coords([Re+alt, lat, lon], "GEO", "sph",["Re", "deg", "deg"])
	cvals.ticks=Ticktock(["2013-10-02T00:00:00"], "ISO")
	mcoords = cvals.convert("MAG", "sph")
	mlat = mcoords.lati[0]
	mlon = mcoords.long[0]
	return mlat, mlon

def polar_to_cart(colat, lon):
	
	"""
	converts polar coordinates to cartesian
	
	Parameters
	----------
	
	colat: float or array of floats
		colatitude of coordinates
		
	lon: float or array of floats
		longitude of coordinates (in degrees)
	"""
	lon_rad = np.deg2rad(lon)
	
	#check if array
	if not isinstance(colat, np.ndarray):
		x = colat*np.sin(lon_rad)
		y = -colat*np.cos(lon_rad)
	else:
		x, y = np.empty(len(colat)), np.empty(len(lon_rad))
		for i in range(len(x)):
			x[i] = colat[i]*np.sin(lon_rad[i])
			y[i] = -colat[i]*np.cos(lon_rad[i])
			
	return x, y

def aacgm_to_mlt(lon, time):
	
	"""
	convert from aacgm into magnetic local time
	
	Parameters
	----------
	
	lon: float or array of floats
		longitude or array of longitudes (in degrees) to convert to mlt 
		
	time: datetime object
		date and time for conversion (UT) in format 
		datetime.datetime(YY, MM, DD, hh, mm, ss)
		
	Returns
	-------
	
	mlt_lon: array of floats
		returns magnetic longitude as either an array-like float or array of
		floats (depending on input)
	"""
	
	mlt_lon = np.array(aacgmv2.convert_mlt(lon, time, m2a=False))
	
	return mlt_lon

def geo_to_aacgm(lat, lon, time, alt=0):
	
	"""
	convert from geographical latitude and longitude into altitude adjusted
	corrected geomagnetic coordinates (aacgm)
	
	Parameters
	----------
	
	lat: float or array of floats
		latitude of coordinates
		
	lon: float or array of floats
		longitude of coordinates (in degrees)
		
	time: datetime object
		datetime object of format datetime.datetime(YY, MM, DD, hh, mm, ss)
		
	Returns
	-------
	
	alat: float or array of floats
		latitude in aacgm
		
	alon: float or array of floats
		longitude in aacgm		
	"""
	
	alat, alon, r = aacgmv2.convert_latlon_arr(lat, lon, alt, time, method_code="G2A")
	
	return alat, alon

def rms(x):
	
	"""
	Returns the root mean square of input x values
	
	Parameters
	----------
	
	x: float array
		values of which to calculate rms
	"""
	
	n = len(x)
	rms = math.sqrt((1/n)*np.sum(x**2))
	return rms
	
#define a sine function
def sin(x, A, B, phi):
	
	"""
	Parameters
	----------
	A = median velocity
	B = amplitude (maximum |velocity|)
	phi = phase shift
	x = angle (in radians)
	"""
	return A+B*np.sin(x+phi)	

#restrict angle of arrays between -90 and 90 degrees (using sine curve)
def sin9090(thetas):
	
	"""
	Uses y = sin(theta) to change thetas >90 and <-90 to their equivalent y 
	values so that all thetas -90 <= theta <= 90
	e.g. 130deg -> 50deg, -110deg -> -70deg
	
	Parameters
	----------
	
	thetas: float array
		array of angles to restrict range between -90 and 90
	"""
	
	new_thetas = np.empty(len(thetas))
	
	for i in range(len(thetas)):
		theta = thetas[i]
		if theta < -90:
			new_thetas[i] = 2*(-90-theta)+theta
		elif theta > 90:
			new_thetas[i] = 2*(90-theta)+theta
		else:
			new_thetas[i] = theta
	
	return new_thetas

def cosine_rule(A, B, C, polar=False):
	
	"""
	Calculates angle of triangle ABC
	
	Parameters
	----------
	A: float array
		array containing x and y position of A [Bx, Ay]
	B: float array
		array containing x and y position of B [Bx, By]
	C: float array
		array containing x and y position of C [Cx, Cy]	
	Polar: Bool
		set to true if using polar lat/lon coordinates, x=lat, y=lon (in degrees)
	"""
	
	if not polar:
		c = math.sqrt((A[0]-B[0])**2 + (A[1]-B[1])**2)
		b = math.sqrt((A[0]-C[0])**2 + (A[1]-C[1])**2)
		a = math.sqrt((B[0]-C[0])**2 + (B[1]-C[1])**2)
		
	elif polar:
		
		
		c = math.sqrt(A[0]**2 + B[0]**2 - 2*A[0]*B[0]*np.cos(np.deg2rad(B[1]-A[1])))
		b = math.sqrt(A[0]**2 + C[0]**2 - 2*A[0]*C[0]*np.cos(np.deg2rad(C[1]-A[1])))
		a = math.sqrt(B[0]**2 + C[0]**2 - 2*B[0]*C[0]*np.cos(np.deg2rad(C[1]-B[1])))
	
	cosB = (c**2 + a**2 - b**2)/(2*c*a)
	B = np.arccos(cosB)
	
	return np.rad2deg(B)

def lon_look(lon0, lon1):
	
	"""
	Determines if lon(gitude)1 is due east or west of lon(gitude)0
	
	Parameters
	----------
	lon0: float
		longitude of reference point (in defrees)
	lon1: float
		longitude of point to determine direction with respect to lon(gitude)0
	"""
	
	#if longitude is given as -180 -> 180 then convert to 0 -> 360
	if lon0 < 0:
		lon0 = 180 + (180-abs(lon0))
		#print("lon0 = {}".format(lon0))
	if lon1 < 0:
		lon1 = 180 + (180-abs(lon0))
		#print("lon1 = {}".format(lon1))
		
	lon_diff = lon1 - lon0
	#print("lon_diff = {}".format(lon_diff))
	
	#set up different conditions depending on lon0 being located in east or west
	if 180 <= lon0 <= 360: #if lon0 is west
		if -180 <= lon_diff <= 0:
			direction = "W"
		elif 0 < lon_diff <= 180:
			direction = "E"
		elif -360 <= lon_diff < -180:
			direction = "E"
		else:
			print("unexpected outcome")
			return
		
	elif 0 <= lon0 < 180: #if lon0 is east#
		if 0 <= lon_diff <= 180:
			direction = "E"
		elif -180 <= lon_diff < 0:
			direction = "W"
		elif 180 < lon_diff <= 360:
			direction = "W"
		else:
			print("unexpected outcome")
			return
	
	return direction
	
def time_to_dtime(date):
	
	"""
	Converts a string time into a datetime object
	
	Parameters:
	-----------
	time: string
		time (in format "YYYY/MM/DD HH:mm:ss") 
	"""
	
	YY = int(date[0:4])
	MM = int(date[5:7])
	DD = int(date[8:10])
	HH = int(date[11:13])
	mm = int(date[14:16])
	ss = int(date[17:19])
	
	YY = int(date[0:4])
	MM = int(date[5:7])
	DD = int(date[8:10])
	HH = int(date[11:13])
	mm = int(date[14:16])
	ss = int(date[17:19])
	
	dtime = dt.datetime(YY, MM, DD, HH, mm, ss)
	
	return dtime
