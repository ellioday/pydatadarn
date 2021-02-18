"""
Created on Thu Nov  5 04:17:39 2020

@author: elliott
"""

# file format = YYYYMMDD.HHmm.ss.bks.fitacf.bz2

import os
import pydarn
import bz2

import scipy.optimize as opt # for optimizing least square fit

import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib as mpl

from pydatadarn.utils import tools
from pydatadarn.utils import coordinate_transformations as coords
from pydatadarn.classes.station import Station

#get path to superdarn data
luna_path_file = open("/home/elliott/Documents/python_analysis/luna_path.txt", "r")
luna_path = luna_path_file.read()
if luna_path[len(luna_path)-1:len(luna_path)] == "\n":
	luna_path = luna_path[0:len(luna_path)-1]
#paths for fitacf and grdmap files
fitacf_path = luna_path + "fitacf/"
grdmap_path = luna_path + "users/daye1/Superdarn/Data/grid/"
	
class LoadGridmap:
	
	"""
	A class used to load fitacf data
	"""
	
	def __init__(self, station, start_date, end_date):
		
		"""
		Returns a list of all file names/paths to fitacf data between the 
		requested dates
		
		Parameters
		----------
		
		station: string
			3-letter station name to load data
			
		start_date: string
			"YYYY/MM/DD/hh.mm.ss"
			
		end_date: string
			"YYYY/MM/DD/hh.mm.ss"
		"""
		
		#extract time parameters into [start, end]
		YY = [int(start_date[0:4]), int(end_date[0:4])]
		MM = [int(start_date[5:7]), int(end_date[5:7])]
		DD = [int(start_date[8:10]), int(end_date[8:10])]
	
		self.path = "{}{}/{}/{}/".format(grdmap_path, station, YY[0], MM[0])	
		print(self.path)
		
		files = sorted(os.listdir(self.path))
		print(files)
		num_files = len(files)
		
		#index where files are within time bounds
		bound_start = False
		bound_end = False
		bounds = [0, 0]
		i = 0
		
		while bound_end == False and i <= num_files:
			
			if i == num_files:
				bounds[1]=i
				break
			#files are named in format: YYYYMMDD.HHmm.ss.bks.fitacf.bz2
		#	print("Checking file number {}...".format(i))
			
			print(files[i])
			file_MM = int(files[i][4:6])
			file_DD = int(files[i][6:8])
			
		#	print("file date = {:02d}/{:02d}".format(file_MM, file_DD))
			
			if not bound_start:
				#check if day is within bound
				if file_MM >= MM[0] and file_DD >= DD[0]:
						bound_start = True
						bounds[0] = i
						print("start file found")
			
			if bound_start:
				#check if month then day and then time is within end bound
				if file_MM > MM[1]:
					bound_end = True
					bounds[1] = i
				elif file_DD > DD[1]:
					bound_end = True
					bounds[1] = i
						
			if bound_end:
				#print("Data retrieved between {} and {}".format(start_date, end_date))
				break
			
			i += 1
			
		else: 
			if not bound_start:
				#print("No files were found between the entered dates")
				return
		
		self.file_list = files[bounds[0]:bounds[1]]
		
		return
	
	def read(self, fname):
		
		"""
		Takes input fitacf file (fname) and returns the data from the file
		
		Parameters
		----------
		
		fname: string
			name of file
		"""
		
		fname = self.path+fname
		
		file_data = pydarn.SDarnRead(fname)
		grid_data = file_data.read_grid()
		
		return grid_data
	
	def read_bz2(self, fname):
		
		"""
		Takes input fitacf.bz2 file (fname) and returns the data from the file
		
		Parameters
		----------
		
		fname: string
			name of file
		"""
		
		fname = self.path+fname
		
		with bz2.open(fname) as fp:
			bz2_stream = fp.read()
			
		file_data = pydarn.SDarnRead(bz2_stream, True)
		grid_data = file_data.read_grid()
		
		return grid_data
	

class GridData():
	
	"""
	A class used to load .grdmap data into easy to access arrays
	
	Parameters
	----------
	
	None
	
	"""
	def __init__(self):

		self.station_metadata = dict() # dictionary to store information for each station
		self.stations = np.array([]) # store which station data comes from
		self.mlats = np.array([])
		self.mlons = np.array([])
		self.kvecs = np.array([])
		self.los_vs = np.array([])
		self.los_e = np.array([]) # los velocity standard deviation
		self.times = np.array([]) # store time as a string
		self.dtimes = np.array([]) # store time as a datetime object
		
		return
	
	def add_station(self, sname, lat, lon, look, boresight = False):
		
		"""
		Adds station data to object
		"""
		
		self.station_metadata[sname] = Station(lat, lon, look)
		if not (not boresight):
			self.station_metadata[sname].set_boresight = boresight
		
		return
	
	def add_data(self, start_date, end_date, mod=True, get_rad_azms=True):
		
		"""
		Adds los data between dates for stations currently in the object
		
		Parameters
		----------
		
		start_date: string
			inclusive start date to get data from 
			(in format "YYYY/MM/DD HH:mm:ss")
			
		end_date: string
			inclusive end date to get data from
			(in format "YYYY/MM/DD HH:mm:ss")
			
		mod: bool
			if True, will apply the modulus/absolute to the los velocities so
			that line of sight velocities are all positive. If false, will set
			velocities towards the radar as negative and velocities away 
			positive (default = True)
		"""
		self.mod = mod
		
		start_YY = int(start_date[0:4])
		start_MM = int(start_date[5:7])
		start_DD = int(start_date[8:10])
		start_HH = int(start_date[11:13])
		start_mm = int(start_date[14:16])
		start_ss = int(start_date[17:19])
		
		end_YY = int(end_date[0:4])
		end_MM = int(end_date[5:7])
		end_DD = int(end_date[8:10])
		end_HH = int(end_date[11:13])
		end_mm = int(end_date[14:16])
		end_ss = int(end_date[17:19])
		
		start_dtime = dt.datetime(start_YY, start_MM, start_DD, start_HH,
							start_mm, start_ss)
		end_dtime = dt.datetime(end_YY, end_MM, end_DD, end_HH, end_mm, end_ss)
		
		for station in self.station_metadata:
			#access files
			self.data_files = LoadGridmap(str(station), start_date, end_date)
			#access data in files
			for i in range(len(self.data_files.file_list)):
				
				data_file = self.data_files.file_list[i]
				if data_file[len(data_file)-4:] == ".bz2":
					self.grid_data = self.data_files.read_bz2(data_file)
				else:
					self.grid_data = self.data_files.read(data_file)
				
				for j in range(len(self.grid_data)):
					#obtain individual sample
					self.sample = self.grid_data[j]
				
					#get data
					try:
	
						mlats = self.sample["vector.mlat"]
						mlons = self.sample["vector.mlon"]
						kvecs = self.sample["vector.kvect"]
						los_vs = self.sample["vector.vel.median"]
						los_e = self.sample["vector.vel.sd"]
						
						YY = int(self.sample["start.year"])
						MM = int(self.sample["start.month"])
						DD = int(self.sample["start.day"])
						hh = int(self.sample["start.hour"])
						mm = int(self.sample["start.minute"])
						ss = int(self.sample["start.second"])
						#access time
						start_day = "{:02d}/{:02d}/{:02d}".format(YY, MM, DD)
						start_time = "{:02d}:{:02d}:{:02d}".format(hh, mm, ss)
						full_time = "{} {}".format(start_day, start_time)
						dtime = dt.datetime(YY, MM, DD, hh, mm, ss)
						#access the look of the station
						look = self.station_metadata[station].look
						self.sample_access = self.sample
						#print("try succesfull")						
						#only save the data if it lies between the requested times
						
						#make sure kvecs are consistent between look directions
						#(aka -ve kvecs are east look, +ve kvecs are west look)
						#use los_v sign to denote direction of flow, negative 
						#velocities = towards radar station.
						if not mod:
							if look == "E":
								#for east look kvec always needs to be negative
								for i in range(len(kvecs)):
									if kvecs[i] < 0:
										los_vs[i] = -los_vs[i]
							elif look == "W":
								#for west look kvec always needs to be positive
								for i in range(len(kvecs)):
									if kvecs[i] > 0:
										los_vs[i] = -los_vs[i]		
						
						print("start_dtime", start_dtime)
						print("dtime", dtime)
						print("end_dtime", end_dtime)
						
						if start_dtime <= dtime <= end_dtime:
							print("time within bounds")
							self.mlats = np.append(self.mlats, mlats)
							self.mlons = np.append(self.mlons, mlons)
							self.kvecs = np.append(self.kvecs, kvecs)
							self.los_vs = np.append(self.los_vs, los_vs)
							self.los_e = np.append(self.los_e, los_e)
							for i in range(len(mlats)):
								self.times = np.append(self.times, full_time)
								self.dtimes = np.append(self.dtimes, dtime)
								self.stations = np.append(self.stations, station)
						else:
							print("time not within bounds")
						
					except:
						continue
			
			print("station {} data added.".format(station))
		
		#calculate colatitudes
		self.mcolats = 90-self.mlats	
						
		#calculate vector and radar azimuths from mag north, radar and vector points
		self.vec_azms = np.empty(len(self.kvecs))
		self.rad_azms = np.empty(len(self.kvecs))
		for i in range(len(self.vec_azms)):

			#get radar 
			radar_name = self.stations[i]
			radar = self.station_metadata[radar_name]
			#get mcolat and mlon of radar for correct time
			radar.set_aacgm(self.dtimes[i])
			
			#get vector position and azimuth
			vec_mcolat = self.mcolats[i]
			vec_mlon = self.mlons[i]
			
			#get vector azms
			vec_azm = tools.cosine_rule([0, 0], [vec_mcolat, vec_mlon], [radar.mcolat, radar.mlon], polar=True)
			#get radar azms
			radar_azm = tools.cosine_rule([0, 0], [radar.mcolat, radar.mlon], [vec_mcolat, vec_mlon], polar=True)
			
			if tools.lon_look(radar.mlon, vec_mlon) == "E":
				radar_azm = -radar_azm
			
			self.vec_azms[i] = vec_azm
			self.rad_azms[i] = radar_azm
			
		#calculate vector azimuths from kvector
		#(-ve vec_azms are east look, +ve vec_azms are west look)
		#then use los_v sign to denote direction of flow, -ve los_v = towards
		#radar
		self.vec_azms = np.array(self.kvecs)
		for i in range(len(self.vec_azms)):
 			look_direction = self.station_metadata[self.stations[i]].look
 			vec_azm = self.vec_azms[i]
 			if look_direction == "E":
 					if vec_azm > 0:
						  self.vec_azms[i] -= 180
 			elif look_direction == "W":	
 					if vec_azm < 0:
						 self.vec_azms[i] += 180				
			
		#restrict azimuths between -90 and 90 degrees
		self.vec_azms = tools.sin9090(self.vec_azms)
		
		return
	
	def get_station_data(self, station_name):
		
		"""
		Returns a dictionary of all data captured from one station
		
		Parameters
		----------
		
		station_name: string
			name of the station to retrieve data for
		"""
		
		#find which indexes correspond to the requested station
		indexes = np.where(self.stations == station_name)
		
		data_dict = dict()
		data_dict["mlats"] = self.mlats[indexes]
		data_dict["mcolats"] = self.mcolats[indexes]
		data_dict["mlons"] = self.mlons[indexes]
		data_dict["kvecs"] = self.kvecs[indexes]
		data_dict["los_vs"] = self.los_vs[indexes]
		data_dict["los_e"] = self.los_e[indexes]
		data_dict["times"] = self.times[indexes]
		data_dict["dtimes"] = self.dtimes[indexes]
		data_dict["vec_azms"] = self.vec_azms[indexes]
		data_dict["rad_azms"] = self.rad_azms[indexes]
		
		return data_dict