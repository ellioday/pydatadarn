"""
Created on Thu Nov  5 04:17:39 2020

@author: elliott
"""

# file format = YYYYMMDD.HHmm.ss.bks.fitacf.bz2

import os
import pydarn
import bz2
import elliotools

import numpy as np
import datetime as dt

#get path to superdarn data
superdarn_data_path = open("/home/elliott/Documents/python_analysis/superdarn_data_path.txt", "r")
superdarn_path = superdarn_data_path.read()
#paths for fitacf and grdmap files
fitacf_path = superdarn_path + "fitacf/"
grdmap_path = superdarn_path + "grdmap/"


class LoadFitacf:
	
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
			
		start_date: string, required
			"YYYY/MM/DD/hh.mm.ss"
			
		end_date: string, required
			"YYYY/MM/DD/hh.mm.ss"
		"""
		
		#extract time parameters into [start, end]
		YY = [int(start_date[0:4]), int(end_date[0:4])]
	
		#convert time to dtime
		start_dtime = elliotools.time_to_dtime(start_date)
		end_dtime = elliotools.time_to_dtime(end_date)
	
		self.path = "{}{}/{}/".format(fitacf_path, station, YY[0])	
		
		files = sorted(os.listdir(self.path))
		print(files)
		num_files = len(files)
		
		#index where files are within time bounds
		bound_start = False
		bound_end = False
		bounds = [0, 0]
		i = 0
		
		while bound_end == False or i < num_files:
			
			#files are named in format: YYYYMMDD.HHmm.ss.bks.fitacf.bz2
			print("Checking validity of file number {}...".format(i))
			
			#some files can be a funny format e.g.:
			#.YYYYMMDD.HHmm.ss.a.bks.fitacf.bz2.sCXyq
			#check the first value in the file name to check that it is of a
			#standard form, if it is not standard, then skip this file
			try:
				int(files[i][0])
			except ValueError:
				i+=1
				continue
			
			file_YY = int(files[i][0:4])
			file_MM = int(files[i][4:6])
			file_DD = int(files[i][6:8])
			file_HH = int(files[i][9:11])
			file_mm = int(files[i][11:13])
			file_ss = int(files[i][14:16])
			file_time = "{:04d}/{:02d}/{:02d} {:02d}:{:02d}:{:02d}.".format(
				file_YY, file_MM, file_DD, file_HH, file_mm, file_ss)
			file_dtime = elliotools.time_to_dtime(file_time)
			
			#print("file date = {:02d}/{:02d}/{:02d}.{:02d}.{:02d}".format(file_MM, file_DD, file_HH, file_mm, file_ss))
			
			if not bound_start:
				#check if day and then time is within bound
				if file_dtime >= start_dtime:
					bound_start = True
					bounds[0] = i
					print("start file found")
			
			if bound_start:
				#check if month then day and then time is within end bound
				if file_dtime > end_dtime:
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
		grid_data = file_data.read_fitacf()
		
		return grid_data
	
	def read_field(self, field, pydarn_data):
		
		"""
		returns the requested field from the input pydarn_data e.g. los_v
		
		Parameters
		----------
		
		field: string
			name of field to access
			
		pydard_data:
			dataset where field is contained e.g. after running running read_bz2
		"""
		
		size = len(pydarn_data)
		
		field_type = type(pydarn_data[0][field])
		print("field_type = {}".format(field_type))
		
		#if field type is a scalar, store each value in an array
		if field_type != np.ndarray:
			data = np.empty(size, type = field_type)
			field_type = type(pydarn_data[0][field][0])
			for i in range(size):
				data[i] = pydarn_data[i][field]
		
		#if the field type is a vector
		elif field == np.ndarray:
			#obtain type inside vector
			field_type = type(pydarn_data[0][field][0])
			#create ndarray of [data_size x vector_size]
			data = np.empty([size, pydarn_data[0][field].shape[0]], dtype = field_type)
			for i in range(size):
				data = np.empty([size, pydarn_data[0][field].shape[0]], dtype=field_type)
	
			
		return data
	
	def fread(self, files):
		
		"""
		Takes an input array of file names (files) and returns all data from each
		file into a dictionary
		
		Parameters
		----------
		
		files: array of strings
			array containing names of files to load
		"""
		
		data = dict()
		
		for i in range(len(files)):
			data[i] = self.read_bz2(self.path+files[i])
		
		self.data = data
			
		return self.data

class FitData():
	
	"""
	A class used to load .fitacf data into easy to access arrays
	
	Paramteers
	----------
	
	None
	"""
	
	def __init__(self):
		
		self.station_metadata = dict()
		self.stations = np.array([])
		self.bmazm = np.array([])
		self.bmnum = np.array([])
		self.nrang = np.array([])
		self.range_gates = np.array([])	
		self.los_vs = np.array([])	
		self.times = np.array([])
		self.dtimes = np.array([])
		
		return
		
	def add_station(self, sname, lat, lon, look):
	
		"""
		Adds station data to object
		"""
		
		self.station_metadata[sname] = Station(lat, lon, look)
	
		return
	
	def add_data(self, start_date, end_date):
		
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
		"""
		
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
			print("station = {}\nstart date = {}\nend date = {}\n".format(station, start_date, end_date))
			self.data_files = LoadFitacf(station, start_date, end_date)
			#access data in files
			for i in range(len(self.data_files.file_list)):
				self.fitacf_data = self.data_files.read_bz2(self.data_files.file_list[i])
				
				for j in range(len(self.fitacf_data)):
					#obtain individual sample
					self.sample = self.fitacf_data[j]
				
					#get data
					try:
	
						bmazm = self.sample["bmazm"]
						bmnum = self.sample["bmnum"]
						range_gates = self.sample["slist"]
						los_vs = self.sample["los_vs"]
						
						YY = int(self.sample["time.yr"])
						MM = int(self.sample["time.mo"])
						DD = int(self.sample["time.dy"])
						hh = int(self.sample["time.hr"])
						mm = int(self.sample["time.mt"])
						ss = int(self.sample["time.sc"])
						
						start_day = "{:02d}/{:02d}/{:02d}".format(YY, MM, DD)
						start_time = "{:02d}:{:02d}:{:02d}".format(hh, mm, ss)
						full_time = "{} {}".format(start_day, start_time)
						dtime = dt.datetime(YY, MM, DD, hh, mm, ss)
						
						if start_dtime <= dtime <= end_dtime:
						
							self.bmazm = np.append(self.bmazm, bmazm)
							self.bmnum = np.append(self.bmnum, bmnum)
							self.range_gates = np.append(self.range_gates, range_gates)
							self.los_vs = np.append(self.los_vs, los_vs)
							for i in range(len(bmazm)):
								self.times = np.append(self.times, full_time)
								self.dtimes = np.append(self.dtimes, dtime)
								self.stations = np.append(self.stations, station)
						
					except:
						continue
			
		return
