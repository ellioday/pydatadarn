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
	range of -90 to 90 degrees e.g. 130deg -> 50deg, -110deg -> -70deg
	
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

def get_radar_azi(radar_lat, radar_lon, vec_lat, vec_lon, vec_azi):
	
	"""
	Calculates the radar look (beam) azimuth of backscatter measurement
	
	Parameters:
	-----------
	radar_lat: float
		latitude of radar
	radar_lon: float
		longitude of radar
	vec_lat: float
		latitude of backscatter measurement
	vec_lon: float
		longitude of backscatter measurement
	vec_azi: float
		azimuth angle (with respect to North) of backscatter measurement
		(in radians)
	"""
	
	#calculate length between North and vector location
	vec_len = math.sqrt(vec_lat**2 + vec_lon**2)
	#calculate length between North and radar location
	radar_len = math.sqrt(radar_lat**2 + radar_lon**2)
	#calculate radar look (beam) azimuth
	radar_azi = (radar_len/np.sin(vec_azi))*vec_len
	
	return radar_azi 