#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 11:41:06 2021

@author: elliott
"""

import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import matplotlib as mpl
import scipy.optimize as opt # for optimizing least square fit

from pydatadarn.utils import tools
from pydatadarn.utils import coordinate_transformations as coords

def vector_plot(mcolats, mlons, kvecs, los_vs, time, station_coords=[], station_names=False, mlt=True):
	
	"""
	Creates a polar plot of line of sight vectors

	
	Parameters
	----------
	
	mcolats: float array
		magnetic colatitude(s) of vector(s) in degrees
		
	mlons: float array
		magnetic longitude(s) of vector(s) in degrees
		
	kvecs: float array
		kvector of vector(s) (in degrees)
		
	los_vs: float array
		line of sight velocities (in m/s)
		
	station_coords: float array
		array containing magnetic colatitude and longitude of stations to plot
		onto vector plot. If multiple stations then use multiple rows for each
		station e.g. np.array([no.stations, 2])
		
	station_names: string array
		array containing names of respective stations in station_coords
		
	time: dtime object, optional
		time of measurement, used to convert from magnetic longitude to 
		magnetic local time (in format "YYYY/MM/DD HH:mm:ss"),
		only needed if mlt = True

	mlt: bool
		sets whether to convert from magnetic longitude into local time
		
	"""
	
	###############################################
	# calculate the change in dr/dtheta for vectors
	##############################################
	
	#calculate scale length
	vec_len = 2*500*abs(los_vs/6371e3)
	
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
	for i in range(len(vec_azm)):
		if vec_azm[i] < 0: delta_lon[i] = -delta_lon[i]
	vec_lon = lon_rad+delta_lon
	
	#calculate change in colatitude and angle (both in degrees)
	dr = np.rad2deg(vec_colat) - mcolats
	dtheta = np.rad2deg(vec_lon - np.deg2rad(mlons))
	
	print(dr[0:10], dtheta[0:10])
	
	####################
	# Plot the vectors #
	####################
	
	#create figure
	fig, ax = plt.subplots(1, 1, figsize=(6, 6),subplot_kw=dict(
		projection="polar"))
	
	#set plot
	ax.set_theta_offset(1.5*np.pi)
	ax.set_xticks(np.linspace(0, 2*np.pi, 4, endpoint=False))

	if mlt:
		dtime = tools.time_to_dtime(time)
		#convert mlons into mlts if so true
		mlons = coords.aacgm_to_mlt(mlons, dtime)*15
		ax.set_xticklabels(["00:00", "06:00", "12:00", "18:00"])
		
	ax.set_rlabel_position(135)
	ax.set_ylim(0, 50)
	ax.set_title("{} UT".format(time))
	ax.set_thetamin(0)
	ax.set_thetamax(90)
	
	#Define normalised scale
	cNorm = mpl.colors.Normalize(vmin=-600, vmax=600)

	cm = mpl.cm.jet
	#Create new axis at right hand side
	ax1 = fig.add_axes([0.9, 0.1, 0.03, 0.8])
	#plot colourmap in created axis
	cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cm, norm=cNorm)
	fig.subplots_adjust(left=0.05, right=0.85)
	
	#plot stations
	for i in range(len(station_coords)):
		station = station_coords[i]
		station_name = station_names[i]
		station_mcolat = station[0]
		station_mlon = station[1]
		if mlt:
			station_mlon = coords.aacgm_to_mlt(station_mlon, time)
			if isinstance(station_mlon, np.ndarray):
				station_mlon = station_mlon[0]
			ax.scatter(np.deg2rad(station_mlon), station_mcolat, 5)
			ax.annotate(station_name, [np.deg2rad(station_mlon), station_mcolat])	
	
	#plot vectors
	ax.quiver(np.deg2rad(mlons), mcolats, np.deg2rad(dtheta), dr, 
		width=0.0015, color=cm(cNorm(los_vs)), angles="xy", 
		scale_units="xy", scale=1)
	
	plt.show()

	return
	
def los_fit(azimuths, los_vs, time, resolution=100, mcolat_range=False, station_names=False, plot=False):	
	
	"""
	uses least square fitting to determine the longitudinal flow velocity
	of given radars
	
	Parameters
	----------
	
	azimuths: float array
		array containig azimuth values (in degrees). Must be between -90 and 90
	los_vs: float array
		array of line of sight velocities
	time: string
		time to plot data for (in format "YYYY/MM/DD HH:mm/ss")
	resolution: int
		number of x points to use in the best fitting model
	mcolat_range: array of float
		array containing either one magnetic colatitude [mcolat] to 
		retrieve data for or a low and high magnetic colatitude 
		[mcolat_low, mcolat_high] to limit the data.
	station_names: string array
		array containing names of respective stations in station_coords 
	plot: Bool
		default = False (no limit)
		set whether to plot the los fit

	"""
		
	#if there is no data or not enough data (need at least equal to the 
	#number of fitting parameters (3))
	if len(azimuths) < 4:
		#print("not enough data for requested time")
		fit = np.array([np.nan])
		w = [np.nan, np.nan, np.nan]
		return
	
	#create x(theta) series
	x_series = np.linspace(-180, 180, resolution, endpoint=True)
	
	#calculate rms for initial/guessed amplitude
	rms = tools.rms(los_vs)
	amp = rms*(2**0.5)
	
	#make the initial fit
	A = np.median(los_vs)
	B = amp
	phi = 0
	trial = tools.sin(np.deg2rad(x_series), A, B, phi)
	#optimise our fit
	
	#params = median velocity, amplitude, phase shift
	params = [A, B, phi]
	w, _ = opt.curve_fit(tools.sin, np.deg2rad(azimuths), 
					  los_vs, params, maxfev=50000)
	#print("Initial parameters = {}".format(params))
	#print("Estimated parameters = {}\n".format(self.w))
	fit = tools.sin(np.deg2rad(x_series), *w)
	
	if plot:		
		
		#make plot of los_v over angle
		fig, ax = plt.subplots(1, 1, figsize=[6, 6])
		
		#get indexes for location of stations
		if station_names:
			
			#get each unique station
			unique_stations = np.unique(station_names)
			if len(unique_stations) == 1:
				ax.scatter(azimuths, los_vs, marker = "+", label = unique_stations[0])
			elif len(unique_stations) == 2:
				station0_indexes = np.where(station_names == unique_stations[0])
				station1_indexes = np.where(station_names == unique_stations[1])
				ax.scatter(azimuths[station0_indexes], los_vs[station0_indexes], 
				   marker="+", label=unique_stations[0])
				ax.scatter(azimuths[station1_indexes], los_vs[station1_indexes], 
				   marker="+", label=unique_stations[1])
			
			ax.legend()
			
		else:
			ax.scatter(azimuths, los_vs, marker = "+")
		
		ax.plot(x_series, fit, "--", color="r")
		ax.set_xlim(-90, 90)
		ax.set_xlabel("azimuth")
		ax.set_ylabel("line of sight velocity (m/s)")
		
		#set title according to mcolat range
		if not (not mcolat_range):
			if len(mcolat_range) == 1:
				ax.set_title("mcolat range: {}: {} UT".format(mcolat_range[0], 
												  time))
			else:
				ax.set_title("mcolat range: {} - {}: {} UT".format(
					mcolat_range[0], mcolat_range[1], time))
		else:
			ax.set_title("mcolat range: {} - {}: {}".format(30.5, 40.5, time))		
			
		plt.show()
	
	return w
