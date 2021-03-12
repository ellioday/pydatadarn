#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 01:11:27 2020

@author: elliott
"""

import math

import datetime as dt
import numpy as np

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
	Returns the y values of sine function with input parameters

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
	
	if not isinstance(date, np.ndarray):

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

	elif isinstance(date, np.ndarray):

		dtimes = np.array([])
		for i in range(len(date)):

			YY = int(date[i][0:4])
			MM = int(date[i][5:7])
			DD = int(date[i][8:10])
			HH = int(date[i][11:13])
			mm = int(date[i][14:16])
			ss = int(date[i][17:19])
			
			YY = int(date[i][0:4])
			MM = int(date[i][5:7])
			DD = int(date[i][8:10])
			HH = int(date[i][11:13])
			mm = int(date[i][14:16])
			ss = int(date[i][17:19])
			
			dtime = dt.datetime(YY, MM, DD, HH, mm, ss)
			dtimes = np.append(dtimes, dtime)

		return dtimes
	
def vector_change(mcolats, mlons, los_vs, kvecs):
	
	"""
	Calculates the change in r and theta (polar coords) for a given vector
	
	mcolats: flt array
		colatitudes of vector location (degrees)
		
	mlons: flt array
		longitudes of vector location (degrees)
		
	los_vs: flt array
		magnitudes of vectors (m/s)
		
	kvecs: flt array
		kvectors of vectors (degrees)
		
	"""
	
	#calculate scale length
	vec_len = 2*500*los_vs/6371e3
	
	#obtain longitude and kvector in radians
	lon_rad = np.deg2rad(mlons)
	vec_azm = np.deg2rad(kvecs)
	
	#find latitude at end of vector
	colat = np.deg2rad(mcolats)
	cos_colat = np.cos(vec_len)*np.cos(colat) + np.sin(vec_len)*np.sin(colat)*np.cos(vec_azm)
	vec_colat = np.arccos(cos_colat)
	
	#find longitude of end of vector
	cos_dlon = (np.cos(vec_len)-np.cos(vec_colat)*np.cos(colat))/(np.sin(vec_colat)*np.sin(colat))
	delta_lon = np.arccos(cos_dlon)
	if isinstance(vec_azm, np.ndarray):
		for i in range(len(vec_azm)):
			if vec_azm[i] < 0: delta_lon[i] = -delta_lon[i]
	else:
		if vec_azm <0: delta_lon = -delta_lon
	vec_lon = lon_rad+delta_lon
	
	#calculate change in colatitude and angle (both in degrees)
	dr = np.rad2deg(vec_colat) - mcolats
	dtheta = np.rad2deg(vec_lon - np.deg2rad(mlons))
	
	return dr, dtheta

def lon360_to_180(lon):
	
	"""
	Converts the given longitude from 0->360 to -180->180 degrees
	
	Parameters
	----------
	
	lon: float
		longitude to convert
	"""
	
	remainder = (lon+180) % 360
	
	return remainder - 180

def lon180_to_360(lon):
	"""
	Converts the given longitude from -180->180 to 0->360
	
	Parameters
	----------
	
	lon: float
		longitude to convert
	"""
	
	return lon % 360