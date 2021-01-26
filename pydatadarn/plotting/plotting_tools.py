#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 03:51:55 2020

@author: elliott
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

class BeamRangeVelocityMap:
	def __init__(self, beam_numbers, range_gates):
		
		self.beam_numbers = beam_numbers
		self.range_gates = range_gates
		
		#create array from 0 to final beam number
		self.beam_range = np.arange(0, beam_numbers[::-1][0]+1)
		
		#create beam_numbers x range_gates array
		self.br_array = np.empty([len(self.beam_range), len(self.range_gates)])
		self.br_array[:] = np.nan
		
		return
		
	def set_time(self, time):
		
		"""
		sets the time of the map
		
		input must be an array of format [YY, MM, DD, HH, mm, ss, us]
		
		YY = year
		MM = month
		DD = day
		HH = hour
		mm = minute
		ss = second
		us = microsecond
		
		"""
		
		YY = time[0]
		MM = time[1]
		DD = time[2]
		HH = time[3]
		mm = time[4]
		ss = time[5]
		us = time[6]
		
		self.time = "{:02d}:{:02d}:{:02d} UT {:02d}/{:02d}/{}".format(HH, mm, ss, DD, MM, YY)
		
		return
	
	def set_velocity(self, velocity, beam_number, range_gate):
		
		self.br_array[beam_number, range_gate] = velocity
		
		return
	
	def plot(self, v_min=-200, v_max=200):
		
		blue = [0/255, 0/255, 254/255]
		light_blue = [36/255, 72/255, 182/255]
		light_green = [71/255, 145/255, 108/255]
		green = [108/255, 217/255, 36/255]
		light_orange = [255/255, 176/255, 1/255]
		orange = [255/255, 117/255, 0/255]
		dark_orange = [255/255, 59/255, 0/255]
		red = [254/255, 0/255, 0/255]
		grey = [128/255, 128/255, 128/255]
		white = [255/255, 255/255, 255/255]
		
		colors = [blue, light_blue, light_green, green, light_orange, orange,
			dark_orange, red]
		colors=colors[::-1]
		
		supderdarn_cmap = ListedColormap(colors)
		
		plt.figure(figsize = [8, 8])
		plt.imshow(self.br_array, cmap = supderdarn_cmap, vmin = v_min, vmax=v_max)
		plt.colorbar()
		plt.axis("scaled")
		#plt.title(self.time)
		
		return
	
	def cmap(self, v_min = -200, v_max = 200):
		
		blue = [0, 0, 254]
		light_blue = [36, 72, 182]
		light_green = [71, 145, 108]
		green = [108, 217, 36]
		light_orange = [255, 176, 1]
		orange = [255, 117, 0]
		dark_orange = [255, 59, 0]
		red = [254, 0, 0]
		grey = [128, 128, 128]
		white = [255, 255, 255]
		
		step = (abs(v_min)+abs(v_max))/8
		steps = np.arange(v_min, v_max, step)
		
		self.cmapping = np.empty([len(self.beam_range), len(self.range_gates), 3], dtype = np.uint8)
		self.cmapping[:] = np.nan
		
		for i in range(len(self.beam_range)):
			for j in range(len(self.range_gates)):
				
				los_v = self.br_array[i][j]
				if abs(los_v) < 20:
					self.cmapping[i][j] = grey
				elif los_v <= steps[0]:
					self.cmapping[i][j] = red
				elif los_v <= steps[1]:
					self.cmapping[i][j] = dark_orange
				elif los_v <= steps[2]:
					self.cmapping[i][j] = orange
				elif los_v <= steps[3]:
					self.cmapping[i][j] = light_orange
				elif los_v <= steps[4]:
					self.cmapping[i][j] = green
				elif los_v <= steps[5]:
					self.cmapping[i][j] = light_green
				elif los_v <= steps[6]:
					self.cmapping[i][j] = light_blue
				elif los_v > steps[7]:
					self.cmapping[i][j] = blue
				else:
					self.cmapping[i][j] = white
					
		plt.figure()
		plt.imshow(self.cmapping)
		plt.title(self.time)
				
		return
		
	def cmap_experimental(self, v_min = -200, v_max = 200):
		
		blue = [0/255, 0/255, 254/255]
		light_blue = [36/255, 72/255, 182/255]
		light_green = [71/255, 145/255, 108/255]
		green = [108/255, 217/255, 36/255]
		light_orange = [255/255, 176/255, 1/255]
		orange = [255/255, 117/255, 0/255]
		dark_orange = [255/255, 59/255, 0/255]
		red = [254/255, 0/255, 0/255]
		grey = [128/255, 128/255, 128/255]
		white = [255/255, 255/255, 255/255]
		
		colors = [blue, light_blue, light_green, green, light_orange, orange,
			dark_orange, red]
		
		color_bar= np.empty([8, 1, 3])
		
		step = (abs(v_min)+abs(v_max))/8
		steps = np.arange(v_min, v_max, step)
		
		self.cmapping = np.empty([len(self.beam_range), len(self.range_gates), 3])
		self.cmapping[:] = np.nan
		
		for i in range(len(self.beam_range)):
			for j in range(len(self.range_gates)):
				
				los_v = self.br_array[i][j]
				if abs(los_v) < 20:
					self.cmapping[i][j] = grey
				elif los_v <= steps[0]:
					self.cmapping[i][j] = red
				elif los_v <= steps[1]:
					self.cmapping[i][j] = dark_orange
				elif los_v <= steps[2]:
					self.cmapping[i][j] = orange
				elif los_v <= steps[3]:
					self.cmapping[i][j] = light_orange
				elif los_v <= steps[4]:
					self.cmapping[i][j] = green
				elif los_v <= steps[5]:
					self.cmapping[i][j] = light_green
				elif los_v <= steps[6]:
					self.cmapping[i][j] = light_blue
				elif los_v > steps[7]:
					self.cmapping[i][j] = blue
				else:
					self.cmapping[i][j] = white
					
		plt.figure(figsize = [8, 8])
		plt.imshow(self.cmapping)
		plt.colorbar(plt.cm.ScalarMappable(cmap=ListedColormap(color_bar[::-1])))
		#plt.colorbar(cmap=ListedColormap(color_bar[::-1]))
		#plt.title(self.time)
				
#		fig = plt.figure(constrained_layout = True)			
# 		spec = gridspec.Gridspec(ncols=2, nrows=1, figure=fig)
# 		ax1 = fig.add_subplot(spec[0, 0])
# 		ax2 = fig.add_subplot(spec[0, 1])
# 		
# 		ax1.imshow(self.cmapping)
# 		ax1.set_title(self.time)
# 		ax2.imshow(color_bar)
# 		ax2.yaxis.tick_right()
# 		ax2.tick_params(axis="x", which="both", bottom = False, top = False, labelbottom = False)
# 		#ax2.set_yticks(np.arange(v_min+step, v_max-step, step))
# 		#plt.title(self.time)
# 				
		
		return
	
def beam_range_velocities(data):
	
	"""
	Takes input (loaded) data file and returns 2-dimensional array of the
	beam_number x range_gate LOS velocities
	"""
	
	#find out how many beams are in the data
	beams_determined = False
	beam_count = 1
	beam_list = [data[0]["bmnum"]] #get first beam number
	while beams_determined == False:
		beam = data[beam_count]["bmnum"]
		beam_list.append(beam)
		beam_count += 1
		#all beam numbers are found when there is a repeat occurence
		if beam == beam_list[0]:
			beams_determined = True
			
	#obtain number of range gates in data
	range_gates = data[0]["nrang"]
	#order beam_ranges and create a list of range_gates
	range_gates_list = np.arange(1, range_gates+1)
	
	#create array of BeamRangeVelocityMaps
	los_maps = []
	beam_cycles = 0
	
	#keep a record of minimum and maximum velocities
	v_min = 0
	v_max = 0
	
	#for all samples
	for i in range(len(data)):
		
		#obtain single sample
		sample = data[i]
		
		#until the beam repeats keep the same map
		#when the beam repeats it indicates a new time map is needed
		beam_number = int(sample["bmnum"])
		
		#if we have fully cycled through the beams then reset
		if beam_number == beam_list[0]:
			beam_cycles += 1
			#create the map
			map_instance = BeamRangeVelocityMap(beam_list, range_gates_list)
			
			if beam_cycles > 1:
				los_maps.append(map_instance)
			
			beam_record = []
		
		beam_record.append(beam_number)
		
		#get LOS velocity, time, beam_numbers and range_gates of the sample
		YY = int(sample["time.yr"])
		MM = int(sample["time.mo"])
		DD = int(sample["time.dy"])
		HH = int(sample["time.hr"])
		mm = int(sample["time.mt"])
		ss = int(sample["time.sc"])
		us = int(sample["time.us"])
		time = [YY, MM, DD, HH, mm, ss, us]
		
		range_gates = sample["slist"]
		los_v = sample["v"]
		
		#for all range gates in sample
		for i in range(len(range_gates)):
			
			if los_v[i] < v_min:
				v_min = los_v[i]
			elif los_v[i] > v_max:
				v_max = los_v[i]
				
			map_instance.set_velocity(float(los_v[i]), beam_number, int(range_gates[i]))
			
		map_instance.set_time(time)
		
	return los_maps, [v_min, v_max]