import aacgmv2

import numpy as np
import datetime as dt

def polar_to_cart(colat, lon):
	
	"""
	converts polar (r, theta) coordinates to cartesian (x, y)
	
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
	wrapper for converting from aacgm into magnetic local time
	
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
	wrapper for converting from geographical latitude and longitude into altitude adjusted
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
