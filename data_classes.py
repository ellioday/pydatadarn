"""
Created on Thu Nov  5 04:17:39 2020

@author: elliott
"""

# file format = YYYYMMDD.HHmm.ss.bks.fitacf.bz2

import os
import pydarn
import bz2
import tools

import scipy.optimize as opt # for optimizing least square fit

import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib as mpl


#get path to superdarn data
data_path = open("/home/elliott/Documents/python_analysis/data_path.txt", "r")
superdarn_path = data_path.read()
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
		MM = [int(start_date[5:7]), int(end_date[5:7])]
		DD = [int(start_date[8:10]), int(end_date[8:10])]
		HH = [int(start_date[11:13]), int(end_date[11:13])]
		mm = [int(start_date[14:16]), int(end_date[14:16])]
		ss = [int(start_date[17:19]), int(end_date[17:19])]
	
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
			
			file_MM = int(files[i][4:6])
			file_DD = int(files[i][6:8])
			file_HH = int(files[i][9:11])
			file_mm = int(files[i][11:13])
			file_ss = int(files[i][14:16])
			
			#print("file date = {:02d}/{:02d}/{:02d}.{:02d}.{:02d}".format(file_MM, file_DD, file_HH, file_mm, file_ss))
			
			if not bound_start:
				#check if day and then time is within bound
				if file_MM >= MM[0] and file_DD >= DD[0]:
					if file_HH >= HH[0] and file_mm >= mm[0] and file_ss >= ss[0]:
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
				elif file_DD == DD[1]:
					if file_HH >= HH[1]:
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
	
		self.path = "{}{}/{}/".format(grdmap_path, station, YY[0])	
		
		files = sorted(os.listdir(self.path))
		num_files = len(files)
		
		#index where files are within time bounds
		bound_start = False
		bound_end = False
		bounds = [0, 0]
		i = 0
		
		while bound_end == False or i < num_files:
			
			#files are named in format: YYYYMMDD.HHmm.ss.bks.fitacf.bz2
		#	print("Checking file number {}...".format(i))
			
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
	
class Station():
	
	
	"""
	creates a station containing data with coordinates in:
		geographical (lat/lon)
		aacgm (mlat/mlon)
	"""
	
	def __init__(self, lat, lon, look):
		
		"""
		Parameters
		----------
		
		lat: float
			geographical latitude
			
		lon: float
			geographical longitude
		"""
		
		self.lat = lat
		self.lon = lon
		self.look = look
		
		return
		
	def set_aacgm(self, time):
		
		"""
		gets aacgmv2 coordinates from geographical
		
		Parameters
		----------
		
		time: datetime object
			datetime object of format datetime.datetime(YY, MM, DD, hh, mm, ss)
		"""
		self.mlat, self.mlon = tools.geo_to_aacgm(self.lat, self.lon, time)
		if isinstance(self.mlat, np.ndarray):
			self.mlat = self.mlat[0]
		if isinstance(self.mlon, np.ndarray):
			self.mlon = self.mlon[0]
		self.mcolat = 90 - self.mlat
		
		return
	
	def get_mlt(self, time):
		
		"""
		gets magnetic local time from magnetic longitude
		
		Parameters
		----------
		time: datetime object
			datetime object of format datetime.datetime(YY, MM, DD, hh, mm, ss)
		"""
		
		mlt = tools.aacgm_to_mlt(self.mlon, time)
		
		if isinstance(mlt, np.ndarray):
			mlt = mlt[0]
		
		return mlt
		
	def set_boresight(self, boresight, time=dt.datetime(2010, 1, 1, 1, 1, 1)):
		
		"""
		Takes geographical boresight and calculates magnetic north boresight
		
		Parameters
		----------
		
		boresight: float
			angle of boresight (in radians)
		time: dtime object
			time of interest
		"""
		
		self.boresight_geo = boresight
		self.boresight_lon = tools.geo_to_aacgm(0, 90, self.boresight)
		
		self.boresight_mag = tools.geo_to_mag_azm(self.boresight_geo, self.lat, self.lon)
		
		return
	
	def set_beam_number(self, beam_number):
		self.beam_number = beam_number
		return
	def set_beam_width(self, beam_width):
		self.beam_width = beam_width
		return
	def set_range_gates(self, range_gates):
		self.range_gates = range_gates
	
		#create polar field of view meshgrid.
		beam_from = self.boresight-(0.5*self.beam_number*self.beam_width)
		beam_to = self.boresight+(0.5*self.beam_number*self.beam_width)
		self.beam_range = np.arange(beam_from, beam_to, self.beam_width)
		self.gates_range = np.arange(0, self.range_gates+1, 1)
		
		return	

class FitacfData():
	
	"""
	A class used to load .fitacf dat into easy to access arrays
	
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
				self.grid_data = self.data_files.read_bz2(self.data_files.file_list[i])
				
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
						
						if start_dtime <= dtime <= end_dtime:
							#print("time within bounds")
							self.mlats = np.append(self.mlats, mlats)
							self.mlons = np.append(self.mlons, mlons)
							self.kvecs = np.append(self.kvecs, kvecs)
							self.los_vs = np.append(self.los_vs, los_vs)
							self.los_e = np.append(self.los_e, los_e)
							for i in range(len(mlats)):
								self.times = np.append(self.times, full_time)
								self.dtimes = np.append(self.dtimes, dtime)
								self.stations = np.append(self.stations, station)
						
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
			
			#get radar azms
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
	
		
		"""
		uses least square fitting to determine the longitudinal flow velocity
		of given radars
		
		Parameters
		----------
		
		stations: array of string
			name of station(s) in a list. Must contain the name of 1 or 2 
			station(s)
		time: string
			time to plot data for (in format "YYYY/MM/DD HH:mm/ss")
		resolution: int
			number of x points to use in the line of sight model
		mcolat_range: array of float
			array containing either one magnetic colatitude [mcolat] to 
			retrieve data for or a low and high magnetic colatitude 
			[mcolat_low, mcolat_high] to limit the data. 
			default = False (no limit)
		plot: Bool
			set whether to plot the los fit
		use_radar_azi: Bool
			set whether to use radar azimuth (default False). If set to False
			then vector azimuth will be used for the fit
			*NOTE* vector azimuth will provide most reliable results for a
			single latitude only.
		"""
		
		#get data for specified stations
		if len(stations) == 1:
			self.station_data = self.get_station_data(stations)
			#add array to tell which station the data belongs to
			self.station_data["station_name"] = np.full(len(self.station_data["mcolats"]), stations[0])
		#if multiple stations are provided then concatenate into one dictionary
		elif len(stations) == 2:
			self.station_data = self.get_station_data(stations[0])
			self.station_data["station_name"] = np.full(len(self.station_data["mcolats"]), stations[0])
			temp_data = self.get_station_data(stations[1])
			temp_data["station_name"] = np.full(len(temp_data["mcolats"]), stations[1])
			print("station_data_1_length = {} station_data_2_length = {}".format(len(self.station_data["mcolats"]), len(temp_data["mcolats"])))
			for key in self.station_data:
				print("key = {}\n{} {}".format(key, len(self.station_data[key]), len(temp_data[key])))
				self.station_data[key] = np.append(self.station_data[key], temp_data[key])
				print("full length = {}".format(len(self.station_data["mcolats"])))	
			print(len(self.station_data["mcolats"]))	
				
		else:
			raise Exception("number of stations needs to be 1 or 2")
			
		#restrict data for requested time
		time_indices = np.where(self.station_data["times"] == time)[0]
		for key in self.station_data:
			self.station_data[key] = self.station_data[key][time_indices]
			
		#restrict data between requested mcolat_ranges
		if not (not mcolat_range):
			if isinstance(mcolat_range, float) or isinstance(mcolat_range, int):
				mcolat_indices = np.where((self.station_data["mcolats"] == mcolat_range))
			elif len(mcolat_range) == 2:
				#get indexes that are withing mcolat_range
				mcolat_indices = np.where((self.station_data["mcolats"] >= mcolat_range[0]) & \
							(self.station_data["mcolats"] <= mcolat_range[1]))
			else:
				raise Exception("mcolat range must be either 1 or 2")
			for key in self.station_data:
				self.station_data[key] = self.station_data[key][mcolat_indices]
		
		#if there is no data or not enough data (need at least equal to the 
		#number of fitting parameters (3))
		if len(stations) == 1:
			if len(self.station_data["los_vs"]) < 4:
				#print("not enough data for requested time")
				self.fit = np.array([np.nan])
				self.w = [np.nan, np.nan, np.nan]
				return
		
		#if two stations are provided then make sure there is enough data to
		#do a fit of both of them (individually)
		#for quality purposes, we will make sure there are twice as many points
		#as parameters (6)
		elif len(stations) == 2:
			station0_index = np.where(self.station_data["station_name"] == stations[0])
			station1_index = np.where(self.station_data["station_name"] == stations[1])
			if len(self.station_data["station_name"][station0_index]) < 7 or len(self.station_data["station_name"][station1_index]) < 7:
				#print("not enough data for requested time")
				self.fit = np.array([np.nan])
				self.w = [np.nan, np.nan, np.nan]
				return
		
		#create x(theta) series
		x_series = np.linspace(-180, 180, resolution, endpoint=True)
		
		if use_radar_azi:
			azms = self.station_data["rad_azms"]
		else:
			azms = self.station_data["vec_azms"]
		
		#calculate rms for initial/guessed amplitude
		rms = tools.rms(self.station_data["los_vs"])
		amp = rms*(2**0.5)
		
		#make the initial fit
		A = np.median(self.station_data["los_vs"])
		B = amp
		phi = 0
		trial = tools.sin(np.deg2rad(x_series), A, B, phi)
		#optimise our fit
		
		#params = median velocity, amplitude, phase shift
		params = [A, B, phi]
		self.w, _ = opt.curve_fit(tools.sin, np.deg2rad(azms), 
						  self.station_data["los_vs"], params, maxfev=50000)
		#print("Initial parameters = {}".format(params))
		#print("Estimated parameters = {}\n".format(self.w))
		self.fit = tools.sin(np.deg2rad(x_series), *self.w)
		
		if plot:
		
			#get indexes for location of stations
			station0_indexes = np.where(self.station_data["station_name"] == stations[0])
			if len(stations) == 2:
				station1_indexes = np.where(self.station_data["station_name"] == stations[1])		
			
			#make plot of los_v over angle
			fig, ax = plt.subplots(1, 1, figsize=[6, 6])		
			ax.scatter(azms[station0_indexes], 
				 self.station_data["los_vs"][station0_indexes], marker="+",
				 label = self.station_data["station_name"][station0_indexes][0])
			#ax.plot(x_series, trial, "--", color="k", label="initial fit")
			ax.plot(x_series, self.fit, "--", color="r")
			ax.set_xlim(-90, 90)
			ax.set_xlabel("azimuth")
			ax.set_ylabel("line of sight velocity (m/s)")
			
			#set title according to mcolat range
			if not (not mcolat_range):
				if isinstance(mcolat_range, float) or isinstance(mcolat_range, int):
					ax.set_title("mcolat range: {}: {} UT".format(mcolat_range, time))
				else:
					ax.set_title("mcolat range: {} - {}: {} UT".format(mcolat_range[0], 
												 mcolat_range[1], time))
			else:
				ax.set_title("mcolat range: {} - {}: {}".format(min(self.station_data["mcolats"]), 
					max(self.station_data["mcolats"]), time))		
			
			#if there was a second station plot that
			if len(stations) == 2:
				ax.scatter(azms[station1_indexes], 
					 self.station_data["los_vs"][station1_indexes], marker="+",
					 label = self.station_data["station_name"][station1_indexes][0])
				
			ax.legend()
			plt.show()
		
		return