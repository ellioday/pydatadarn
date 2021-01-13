#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 01:04:55 2020

@author: elliott
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 08:28:12 2020

@author: elliott
"""


import numpy as np
import data_classes
import matplotlib.pyplot as plt
from scipy import stats # for finding mode
import scipy.optimize as opt # for optimizing least square fit

import tools
import data_classes

from mpl_toolkits.axes_grid1 import make_axes_locatable

vectors = data_classes.GridData()

#add bks data
vectors.add_station("bks", 37.1, -77.95, "W")
#vectors.add_station("ade", 51.89, -176.63, "E")
#vectors.add_station("adw", 51.89, -176.63, "W")
#vectors.add_station("fhe", 38.859, -99.389, "E")
#vectors.add_station("fhw", 38.859, -99.389, "W")
#vectors.add_station("cve", 43.271, -120.358, "E")
#vectors.add_station("cvw", 43.271, -120.358, "W")
vectors.add_station("wal", 37.93, -75.47, "E")
vectors.add_data("2013/10/02 00:00:00", "2013/10/02 09:00:00", mod=False)

time_i = "2013/10/02 07:58:00"

#vectors.vector_plot(time_i)
vectors.los_fit(["bks", "wal"], time_i, plot=True)

times = np.array([])

#create set of parameters to extract data for following: [low, high, step]
hours = np.arange(0, 9, 1)
minutes = np.arange(0, 60, 2)
mcolats = np.arange(min(vectors.mcolats), max(vectors.mcolats)+1, 1)

#calculate how many time measurements we have done
no_measurements = len(hours)*len(minutes)
#calculate number of latitude measurements done
no_latitudes = len(mcolats)
velocities_fit = np.empty([no_latitudes, no_measurements])
velocities_fit.fill(np.nan)
velocities_raw = np.empty([no_latitudes, no_measurements])
velocities_raw.fill(np.nan)
time_count = 0

count = 0
for hour in hours:
	for minute in minutes:
		time = "2013/10/02 {:02d}:{:02d}:00".format(hour, minute)
		times = np.append(times, time)
		
		#Find indexes of matching raw velocites for this time
		velocities_raw_time_sample_index = np.where(vectors.times == time)[0]
		if len(velocities_raw_time_sample_index) == 0:
			time_count += 1
			continue
		#get velocites and mcolats for this time
		velocities_raw_time_sampled = vectors.los_vs[velocities_raw_time_sample_index]
		mcolats_raw_time_sampled = vectors.mcolats[velocities_raw_time_sample_index]

		#get which mcolats had data for this time sample
		mcolats_time_sampled = np.arange(min(mcolats_raw_time_sampled), max(mcolats_raw_time_sampled)+1, 1)
		for mcolat in mcolats_time_sampled:
			
			mcolat_index = np.where(mcolats == mcolat)[0][0]
			
			#get raw velocities for time and mcolat sample
			velocities_raw_time_mcolat_sampled_index = np.where(mcolats_raw_time_sampled == mcolat)[0]
			velocities_raw_time_mcolat_sampled = velocities_raw_time_sampled[velocities_raw_time_mcolat_sampled_index]
			#if no data then set as nan
			if len(velocities_raw_time_mcolat_sampled) == 0:
				velocities_raw_time_mcolat_sampled = [np.nan]
			#get and save largest velocity for time and mcolat sample
			velocity_raw_time_mcolat_sampled_max = np.nanmedian(velocities_raw_time_mcolat_sampled)
			velocities_raw[mcolat_index, time_count] = velocity_raw_time_mcolat_sampled_max
			
			#get velocites from best fit
			vectors.los_fit(["bks", "wal"], time, mcolat_range=mcolat)
			#obtain velocity magnitude from fit
			vel_mag = vectors.w[1]
			#save to array
			velocities_fit[mcolat_index, time_count] = vel_mag
			
			print("{} - {} - [{}, {}] - {} - {}".format(time, mcolat, mcolat_index, time_count, velocity_raw_time_mcolat_sampled_max, vel_mag))

		time_count += 1

#plot latitude against velocities with time

def make_colorbar_with_padding(ax):
	"""
	Create colorbar axis that fits the size of a plot
	"""
	
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="5%", pad=0.1)
	return cax

data_indexes = np.array([], dtype ="int")
#get only latitudes that had velocity measurements
for i in range(len(mcolats)):
	if False in np.isnan(velocities_fit[i]):
		print("data found at {}".format(i))
		data_indexes = np.append(data_indexes, i)
		
mcolat_indexes = mcolats[data_indexes]
velocities_indexed = velocities_fit[data_indexes]


#plot latitude vs time velocities
fig = plt.figure(figsize=[6, 1])
#get latitude vs velocity array

ax = fig.add_subplot(311)
plt.imshow(abs(velocities_raw), cmap="jet")
ax.set_aspect("auto")
#make the colourbar
cax = make_colorbar_with_padding(ax) #add colorbar within
cb0 = plt.colorbar(cax=cax, cmap = "jet")
cax.set_ylabel("Velocity (m/s)", rotation = -90)
#set y axis
ax.set_ylabel("Magnetic Latitude")
ax.set_yticks((mcolats-mcolats[0])[::5]) # must be set for ax.set_yticklabels to work properly
ax.set_yticklabels((90-mcolats)[::5])
#set x axis
ax.axes.get_xaxis().set_ticks([])

ax0 = fig.add_subplot(312)
plt.imshow(abs(velocities_fit), cmap="jet")
ax0.set_aspect("auto")
#make the colourbar
cax0 = make_colorbar_with_padding(ax0) #add colorbar within
cb0 = plt.colorbar(cax=cax0, cmap = "jet")
cax0.set_ylabel("Velocity (m/s)", rotation = -90)
#set y axis
ax0.set_ylabel("Magnetic Latitude")
ax.set_yticks((mcolats-mcolats[0])[::5]) # must be set for ax.set_yticklabels to work properly
ax.set_yticklabels((90-mcolats)[::5])
#set x axis
ax0.axes.get_xaxis().set_ticks([])

#get median velocity for each time
vel_median = np.empty(time_count)
for i in range(len(vel_median)):
	vel_median[i] = np.nanmedian(abs(velocities_indexed[:,i]))

ax1 = fig.add_subplot(313)
ax1.plot(np.arange(0, time_count, 1), vel_median)
ax1.set_xlim(0, time_count)
#padd empty space so x axis aligns with ax2
cax1 = make_colorbar_with_padding(ax1)
cax1.spines["top"].set_visible(False)
cax1.spines["right"].set_visible(False)
cax1.spines["bottom"].set_visible(False)
cax1.spines["left"].set_visible(False)
cax1.axes.get_xaxis().set_visible(False)
cax1.axes.get_yaxis().set_visible(False)
#set y axis
ax1.set_ylabel("Velocity (m/s)")
#set x axis
time_labels = []
hour_labels = np.append(hours, hours[len(hours)-1]+1)
for hour in hour_labels:
	time_labels.append("{:02d}".format(hour))
ax1.set_xticks(hour_labels*len(minutes))
ax1.set_xticklabels(time_labels)
ax1.set_xlabel("Time (UT)")
ax1.xaxis.set_tick_params(labeltop="on")

fig.subplots_adjust(hspace=0)

plt.show()
