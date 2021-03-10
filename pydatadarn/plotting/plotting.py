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
import cartopy.crs as ccrs

from pydatadarn.classes.station import Station
from pydatadarn.utils import tools as tools
from pydatadarn.utils import coordinate_transformations as coords

import fpipy
import cartopy
import aacgmv2

def vector_plot(mcolats, mlons, kvecs, los_vs, time, 
				station_names=[], FPI_names=[], FPI_kvecs=[], FPI_vels=[], boundary_mlats=[], boundary_mlons=[], 
				mlt=True, cart=False, mcolat_min=0, mcolat_max=50, theta_min=0, 
				theta_max=360, cbar_min=-600, cbar_max=600):
	
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
		
	station_names (optional): string array
		array containing names of respective stations in station_coords
		
	time: str, optional
		time of measurement, used to convert from magnetic longitude to 
		magnetic local time (in format "YYYY/MM/DD HH:mm:ss"),
		only needed if mlt = True
		
	station_names (optional): str array
		array of names of superDARN stations to plot		
		
	FPI_names (optional): str array
		array of names of FPI stations to plot	
		
	FPI_kvecs (optional): str array
		array of 2d Neutral vector kvectors (in same order as FPI_names)
		
	FPI_vels (optional): str array
		array of 2d neutral vector velocities (in same order as FPI_names)
		
	boundary_mlats (optional): float array
		array of boundary_mlats (degrees)
		
	boundary_mlons (optional): float array
		array of boundary_mlons (degrees)

	mlt (optional): bool
		sets whether to convert from magnetic longitude into local time 
		(default True)
		
	cart (optional): bool
		sets whether to plot over a cartographic map (defualt false)
	
	mcolat_min: float
		sets the minimum magnetic colatitude to be shown by the plot
	
	mcolat_max: float
		sets the maxmimum magnetic colatitude to be shown by the plot
		
	theta_min: float
		sets the minimum angle to be shown by the plot .
		(0, 90, 180, 270, 360 degrees = South, East, North, West, South 
		   if viewing plot as a compass)
		
	theta_max: float
		sets the maximum angle to be shown by the plot
		(0, 90, 180, 270, 360 degrees = South, East, North, West, South 
		   if viewing plot as a compass)
		
	cbar_min: float
		sets the minimum value for the colourbar (default -600 m/s)
		
	cbar_max: float
		sets the maximum value for the colourbar (default 600 m/s)
		
	"""
	
	####################
	# Plot the vectors #
	####################
	
	dtime = tools.time_to_dtime(time)
	
	if cart == False:
		#create figure
		fig, ax = plt.subplots(1, 1, figsize=(6, 6),subplot_kw=dict(
			projection="polar"))
		
		#set plot
		ax.set_theta_offset(1.5*np.pi)
		ax.set_xticks(np.linspace(0, 2*np.pi, 4, endpoint=False))
	
		if mlt == True:
			#convert mlons into mlts if so true
			mlons = coords.aacgm_to_mlt(mlons, dtime)*15
			ax.set_xticklabels(["00:00", "06:00", "12:00", "18:00"])
		
		ax.set_rlabel_position(135)
		ax.set_ylim(mcolat_min, mcolat_max)
		ax.set_title("{} UT".format(time))
		ax.set_thetamin(theta_min)
		ax.set_thetamax(theta_max)
	
	elif cart == True:
		
		fig, ax = plt.subplots(1, 1, figsize=(6, 6),subplot_kw=dict(projection=ccrs.Orthographic(-90, 90)))
		#ax.add_feature(cartopy.feature.LAKES)
		#ax.add_feature(cartopy.feature.COASTLINE)
		ax.set_global()
		ax.coastlines(color="gray")
		ax.gridlines()
	
	#Define normalised scale
	cNorm = mpl.colors.Normalize(vmin=cbar_min, vmax=cbar_max)

	cm = mpl.cm.jet
	#Create new axis at right hand side
	ax1 = fig.add_axes([0.9, 0.1, 0.03, 0.8])
	#plot colourmap in created axis
	cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cm, norm=cNorm)
	fig.subplots_adjust(left=0.05, right=0.85)
	
	#plot SuperDarn stations
	for i in range(len(station_names)):
		#get station data
		station_name = station_names[i]
		station_hdw = Station(station_name)
		#get station mcolat and mlon
		station_mlat, station_mlon =  station_hdw.get_coords(dtime, aacgm=True)
		station_mcolat = 90-station_mlat
		
		if cart == False:
			if mlt == True:
				station_mlon = coords.aacgm_to_mlt(station_mlon, dtime)*15
				if isinstance(station_mlon, np.ndarray):
					station_mlon = station_mlon[0]
			ax.scatter(np.deg2rad(station_mlon), station_mcolat, s=5)
			ax.annotate(station_name, [np.deg2rad(station_mlon), station_mcolat])	
			
		elif cart == True:
			station_lat, station_lon = station_hdw.get_coords(dtime, aacgm=False)
			station_colat = 90-station_lat
			ax.scatter(station_lon, station_lat, s=5, transform=ccrs.PlateCarree())
			ax.text(station_lon, station_lat, station_name, transform=ccrs.PlateCarree())
		
	FPI_mcolats = np.empty(len(FPI_names))
	FPI_mlons = np.empty(len(FPI_names))	
		
	#plot FPI stations
	for i in range(len(FPI_names)):
		#get FPI data
		FPI_name = FPI_names[i]
		FPI_hdw = fpipy.FPIStation(FPI_name)
		#get station mcolat and mlon
		FPI_mlat, FPI_mlon = FPI_hdw.get_coords(dtime, aacgm=True)
		FPI_mcolat = 90-FPI_mlat
		
		if cart == False:
			if mlt == True:
				FPI_mlon = coords.aacgm_to_mlt(FPI_mlon, dtime)*15
				if isinstance(FPI_mlon, np.ndarray):
					FPI_mlon = FPI_mlon[0]
			FPI_mcolats[i] = FPI_mcolat
			FPI_mlons[i] = FPI_mlon		
			ax.scatter(np.deg2rad(FPI_mlon), FPI_mcolat, marker="^", s=5)
			ax.annotate(FPI_name, [np.deg2rad(FPI_mlon), FPI_mcolat])
			
		if cart == True:
			FPI_lat, FPI_lon = FPI_hdw.get_coords(dtime, aacgm=False)
			ax.scatter(FPI_lon, FPI_lat, s=5, transform=ccrs.PlateCarree())
			ax.text(FPI_lon, FPI_lat, FPI_name, transform=ccrs.PlateCarree())
			
	#####################################
	### plot Heppner-Maynard Boundary ###
	#####################################
	
	if len(boundary_mlats) != len(boundary_mlons):
		print("boundary_mlats and boundary_mlons must be the same length")
		plt.close()
		return
		
	if cart == False:
		if mlt == True:
			boundary_mlons = coords.aacgm_to_mlt(boundary_mlons, dtime)*15
		ax.plot(np.deg2rad(boundary_mlons), 90-boundary_mlats, color="k", linestyle="--")
		
	elif cart == True:
		#convert hmb from aacgm to geographic
		boundary_lats, boundary_lons, alt = aacgmv2.convert_latlon_arr(boundary_mlats, boundary_mlons, 0, dtime)
		ax.plot(boundary_lons, boundary_lats, color="black", linestyle="--", transform=ccrs.PlateCarree())
		ax.scatter(0, 60, transform=ccrs.PlateCarree())
		ax.scatter(90, 60, transform=ccrs.PlateCarree())
		print("boundary_lons", boundary_lons, "\n")
		print("boundary_lats", boundary_lats, "\n")
	
	###############################################
	# calculate the change in dr/dtheta for vectors
	##############################################
	
	print("superDARN changes...")
	dr, dtheta = tools.vector_change(mcolats, mlons, los_vs, kvecs)
	
	print("FPI_mcolats", FPI_mcolats)
	print("FPI_mlons", FPI_mlons)
	print("FPI_vels", FPI_vels)
	print("FPI_kvecs", FPI_kvecs)
	
	if len(FPI_kvecs) > 0:
		FPI_dr, FPI_dtheta = tools.vector_change(FPI_mcolats, FPI_mlons, FPI_vels, FPI_kvecs)	
	
		print("FPI_dr", FPI_dr)
		print("FPI_dtheta", FPI_dtheta)
	
	#################
	# Make the Plot #
	#################
	
	if len(FPI_kvecs) > 0:
		#we need to scale between both the FPI and velocity data
		print("plotting FPI_vector")
		#plot vectors
		ax.quiver(np.deg2rad(FPI_mlons), FPI_mcolats, np.deg2rad(FPI_dtheta), FPI_dr, 
			width=0.0015, color="k", angles="xy", 
			scale_units="xy", scale=1)	
	
	#plot vectors
	print("plotting superDarn vectors")
	if cart == False:
		ax.scatter(np.deg2rad(mlons), mcolats, color=cm(cNorm(los_vs)), s=5)
		ax.quiver(np.deg2rad(mlons), mcolats, np.deg2rad(dtheta), dr, 
			width=0.0015, color=cm(cNorm(los_vs)), angles="xy", 
			scale_units="xy", scale=1, headaxislength=0)
		
	elif cart == True:
		#convert aacgm mlons/mlats in geographic lons/lats
		lats, lons, alts = aacgmv2.convert_latlon_arr(90-mcolats, mlons, 0, dtime, method_code="A2G")
		colats = 90-lats
		#our dr and dtheta have been calculated with respect to colatitude so
		#calculate end of vectors with respect to latitude instead
		colats_end = colats + dr
		lats_end = 90-colats_end
		lats_dr = lats_end-lats
		
		ax.scatter(lons, lats, color=cm(cNorm(los_vs)), s=5, transform=ccrs.PlateCarree())
		print("lons", lons)
		ax.quiver(lons, lats, dtheta, lats_dr, 
			width=0.0015, color=cm(cNorm(los_vs)), angles="xy", 
			scale_units="xy", headaxislength=0, transform=ccrs.PlateCarree())
	
	plt.show()

	print("\n")

	return dr, dtheta
	
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
		if station_names != False:
			
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