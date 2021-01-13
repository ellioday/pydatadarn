#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 08:46:24 2020

@author: elliott
"""
import matplotlib.pyplot as plt
from utils import geo_to_mag, mag_to_geo, mag_to_geo_array
from polar_map import Station
from data_classes import gridmap as GridData
from mpl_toolkits.basemap import Basemap
import numpy as np

bks_record = GridData("bks", "2013/10/02", "2013/10/02")
bks_data = bks_record.read_bz2(bks_record.file_list[0])

sample = bks_data[307]

"""
set up BKS station data
"""

BKS = Station(37.1, -77.95)
BKS.set_beam_number(24)
BKS.set_beam_width(3.24)
BKS.set_boresight(-40)
BKS.set_range_gates(110)

"""
Set up basemap coordinate system
"""

#obtain geographical coordinates of magnetic dipole ([0 mlat, 0 mlon])
geographic_dipole = [0, 0]
magnetic_dipole = [0, 0]
magnetic_dipole_coords = mag_to_geo(0, 90, 0)
magnetic_dipole[0]=magnetic_dipole_coords.lati[0]
magnetic_dipole[1]=magnetic_dipole_coords.long[0]

#set how far to go bound (in mlat)
magnetic_bounding_lat = 50
magnetic_bounding_coords = mag_to_geo(0, magnetic_bounding_lat, 0)
magnetic_bounding=magnetic_bounding_coords.lati[0]

#generate the basemap
m = Basemap(projection="ortho", lon_0=magnetic_dipole[1],
			lat_0=magnetic_dipole[0], resolution="c")
m.drawmapboundary()

"""
generate the 1 degree mlat, mlong grid
"""

#generate colatitude bin points
mcolatitude_bins = np.arange(0, int(magnetic_bounding_lat), 1)
#convert to geographic
colatitude_bins = mag_to_geo_array(0, mcolatitude_bins, 
 											np.zeros(len(mcolatitude_bins)))
colatitude_bins = colatitude_bins[0]
#convert colatitudes to radians
colatitude_bins_r = (colatitude_bins/180)*np.pi
#calculate number of longitude bins at each colatitude to acheive 1 degree accuracy
n_long_bins = np.empty(len(colatitude_bins), dtype = int)
for i in range(len(colatitude_bins)):
 	n_long_bins[i] = int(round(360*np.sin(colatitude_bins_r[i])))
 	
#convert colatitudes into latitudes
latitude_bins = 90 - colatitude_bins

for i in range(len(colatitude_bins)):
 	n_bins = n_long_bins[i]
 	for j in range(n_bins):
		  m.scatter(latitude_bins[i], j/360, latlon=True)

"""
obtain the data
"""

#obtain a 2-min sample of data
sample = bks_data[300]
mlats = sample["vector.mlat"]
mlons = sample["vector.mlon"]
kvecs =sample["vector.kvect"]

lats = np.empty([len(mlats)])
lons = np.empty([len(mlons)])

for i in range(len(lats)):
	coords = mag_to_geo(200, mlats[i], mlons[i])
	lats[i] = coords.lati[0]
	lons[i] = coords.long[0]

#add the bks radar station
m.scatter(BKS.long, BKS.lat, 5, latlon=True)
m.scatter(lons, lats, 1, latlon=True)
m.drawcoastlines()

plt.show()