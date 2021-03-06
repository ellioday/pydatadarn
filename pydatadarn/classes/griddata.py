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
luna_path_file = open("/home/elliott/Documents/python_analysis/luna_path.txt", "r")
luna_path = luna_path_file.read()
if luna_path[len(luna_path)-1:len(luna_path)] == "\n":
	luna_path = luna_path[0:len(luna_path)-1]
	
#paths for fitacf and grid/map files
fitacf_path = luna_path + "fitacf/"
grid_path = luna_path + "users/daye1/Superdarn/Data/grid/"
map_path = luna_path + "users/daye1/Superdarn/Data/map/"
	
class LoadGridmap:
	
	"""
	A class used to load fitacf data
	"""
	
	def __init__(self, station, start_date, end_date, is_map=False):
		
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
			
		is_map: bool (optional)
			sets whether files are .grd (True) or .map (False) files 
			(default = False)
		"""
		
		#extract time parameters into [start, end]
		YY = [int(start_date[0:4]), int(end_date[0:4])]
		MM = [int(start_date[5:7]), int(end_date[5:7])]
		DD = [int(start_date[8:10]), int(end_date[8:10])]
	
		if is_map == False:
			self.path = "{}{}/{}/{}/".format(grid_path, station, YY[0], MM[0])	
			files = sorted(os.listdir(self.path))
		else:
			self.path = "{}{}/{}/".format(map_path, YY[0], MM[0])
			all_files = sorted(os.listdir(self.path))
			files = [i for i in all_files if "north" in i]
			
		print(self.path)	

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
			file_MM = int(files[i][4:6])
			file_DD = int(files[i][6:8])

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
				break
			i += 1
			
		else: 
			if not bound_start:
				return
		
		self.file_list = files[bounds[0]:bounds[1]]
		
		return
	
	def read(self, fname, is_map=False):
		
		"""
		Takes input fitacf file (fname) and returns the data from the file
		
		Parameters
		----------
		
		fname: string
			name of file
		"""
		
		fname = self.path+fname
		
		file_data = pydarn.SDarnRead(fname)
		if is_map == False:
			grid_data = file_data.read_grid()
		else:
			grid_data = file_data.read_map()
		
		return grid_data
	
	def read_bz2(self, fname, is_map=False):
		
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
		if is_map == False:
			grid_data = file_data.read_grid()
		else:
			grid_data = file_data.read_map()
		
		return grid_data
	
class GridData():
	
	"""
	A class used to load .grdmap data into easy to access arrays
	
	Parameters
	----------
	
	is_map: bool (optional)
		sets whether files are .grd (True) or .map (False) files 
		(default = False)
	
	"""
	
	def __init__(self, is_map=False):

		self.is_map = is_map
		self.station_metadata = dict() # dictionary to store information for each station
		self.stations = np.array([]) # store which station data comes from
		self.mlats = np.array([])
		self.mlons = np.array([])
		self.kvecs = np.array([])
		self.los_vs = np.array([])
		self.look = np.array([]) # for determining if data are due east or west of station
		self.los_e = np.array([]) # los velocity standard deviation
		self.times = np.array([]) # store time as a string
		self.dtimes = np.array([]) # store time as a datetime object
		
		#if this is map data then add extra arrays
		if is_map == True:
			self.order = dict()
			self.N = dict()
			self.N1 = dict()
			self.N2 = dict()
			self.N3 = dict()
			self.imfBx = dict()
			self.imfBy = dict()
			self.imfBz = dict()
			self.imfTilt = dict()
			self.boundary_mlats = dict()
			self.boundary_mlons = dict()
			self.mod_kvecs = dict()
			self.mod_mlats = dict()
			self.mod_mcolats = dict()
			self.mod_mlons = dict()
			self.mod_los_vs = dict()
			self.mod_times = dict()
			self.mod_dtimes = dict()
			self.latmin = dict()
			self.lon_shft = dict()
			self.lat_shft = dict()
			
		return
	
	def add_data(self, sname, start_date, end_date, mod=True, get_rad_azms=True):
		
		"""
		Adds los data between dates for stations currently in the object
		
		Parameters
		----------
		
		sname: string
			3-letter station code for radar data to add e.g. "ade". to add data
			from all available radars use "all"
		
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
		
		self.data_files = LoadGridmap(str(sname), start_date, end_date, self.is_map)
		
		#access data in files
		for i in range(len(self.data_files.file_list)):
			
			data_file = self.data_files.file_list[i]
			
			if data_file[len(data_file)-4:] == ".bz2":
				self.grid_data = self.data_files.read_bz2(data_file, self.is_map)
			else:
				self.grid_data = self.data_files.read(data_file, self.is_map)
			
			#obtain 2-min sample
			for j in range(len(self.grid_data)):
				self.sample = self.grid_data[j]
			
				#get data
				try:

					#los_v is +ve if pointing Westwards
					#los_v is -ve if pointing Eastwards					

					mlats = self.sample["vector.mlat"]
					mlons = self.sample["vector.mlon"]
					kvecs = self.sample["vector.kvect"]
					los_vs = self.sample["vector.vel.median"]
					los_e = self.sample["vector.vel.sd"]
					
					if self.is_map == True:
						
						order = self.sample["fit.order"]
						N = self.sample["N"]
						N1 = self.sample["N+1"]
						N2 = self.sample["N+2"]
						N3 = self.sample["N+3"]
						imfBx = self.sample["IMF.Bx"]
						imfBy = self.sample["IMF.By"]
						imfBz = self.sample["IMF.Bz"]
						imfTilt = self.sample["IMF.tilt"]
						boundary_mlats = self.sample["boundary.mlat"]
						boundary_mlons = self.sample["boundary.mlon"]
						mod_mlats = self.sample["model.mlat"]
						mod_mlons = self.sample["model.mlon"]
						mod_kvecs = self.sample["model.kvect"]
						mod_los_vs = self.sample["model.vel.median"]
						latmin = self.sample["latmin"]
						lon_shft = self.sample["lon.shft"]
						lat_shft = self.sample["lat.shft"]

					#access time
					YY = int(self.sample["start.year"])
					MM = int(self.sample["start.month"])
					DD = int(self.sample["start.day"])
					hh = int(self.sample["start.hour"])
					mm = int(self.sample["start.minute"])
					ss = int(self.sample["start.second"])
					start_day = "{:02d}/{:02d}/{:02d}".format(YY, MM, DD)
					start_time = "{:02d}:{:02d}:{:02d}".format(hh, mm, ss)
					full_time = "{} {}".format(start_day, start_time)
					dtime = dt.datetime(YY, MM, DD, hh, mm, ss)

					#add skip and stop condition
					if dtime < start_dtime:
						continue
					if dtime > end_dtime:
						break

					#if using grid of all stations then we do not know which
					#station the data belongs to
					if sname != "all":
						
						look = np.array([])
						self.station_metadata[sname] = Station(sname)
	
						#find which station data to use (time)
						station = self.station_metadata[sname]
						#calculate if the data point (mlon) is due east or
						#west from the station
						for i in range(len(mlats)):
							#get station magnetic coordinates
							station_mlat, station_mlon = station.get_aacgm(dtime)
							#get look direction
							look_direction = elliotools.lon_look(station_mlon, mlons[i])	
							look = np.append(look, look_direction)	
					
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
							self.stations = np.append(self.stations, sname)
						if sname != "all":
							self.look = np.append(self.look, look)
						
						if self.is_map == True:
							self.order[full_time] = order
							self.N[full_time] = N
							self.N1[full_time] = N1
							self.N2[full_time] = N2
							self.N3[full_time] = N3
							self.imfBx[full_time] = imfBx
							self.imfBy[full_time] = imfBy
							self.imfBz[full_time] = imfBz
							self.imfTilt[full_time] = imfTilt
							self.boundary_mlats[full_time] = boundary_mlats
							self.boundary_mlons[full_time] = boundary_mlons
							self.mod_mlats[full_time] = mod_mlats
							self.mod_mcolats[full_time] = 90 - mod_mlats
							self.mod_mlons[full_time] = mod_mlons
							self.mod_kvecs[full_time] = mod_kvecs
							self.mod_los_vs[full_time] = mod_los_vs
							self.latmin[full_time] = latmin
							self.lon_shft[full_time] = lon_shft
							self.lat_shft[full_time] = lat_shft
							self.mod_dtimes[full_time] = dtime
						
					else:
						print("time not within bounds")
					
				except:
					continue
		
		#calculate magnetic colatitudes
		self.mcolats = 90-self.mlats
		
		return
		
	def get_azms(self):
		
		"""
		Calculates both vector and radar azimuths for velocity data. ONLY WORKS
		IF DATA HAS BEEN GIVEN AS INDIVIDUAL STATIONS
		
		also converts los_vs so that instead of -ve = eastwards and +ve = westwards
		now, -ve = towards radar and +ve = away,
		"""
						
		#calculate vector and radar azimuths from mag north, radar and vector points
		self.vec_azms = np.empty(len(self.kvecs))
		self.rad_azms = np.empty(len(self.kvecs))
		self.los_vs_rad_point = np.empty(len(self.kvecs))
		for i in range(len(self.vec_azms)):

			#get radar 
			radar_name = self.stations[i]
			radar = self.station_metadata[radar_name]
			#get mcolat and mlon of radar for correct time
			radar_mlat, radar_mlon = radar.get_aacgm(self.dtimes[i])
			radar_mcolat = 90-radar_mlat
			
			#get vector position and azimuth
			vec_mcolat = self.mcolats[i]
			vec_mlon = self.mlons[i]
			los_vs_rad_point = self.los_vs[i]
			
			#get vector azms
			vec_azm = elliotools.cosine_rule([0, 0], [vec_mcolat, vec_mlon], [radar_mcolat, radar_mlon], polar=True)
			#get radar azms
			radar_azm = elliotools.cosine_rule([0, 0], [radar_mcolat, radar_mlon], [vec_mcolat, vec_mlon], polar=True)
			#get look direction from radar to los measurement
			radar_look = elliotools.lon_look(radar_mlon, vec_mlon)
			
			if radar_look == "E":
				radar_azm = -radar_azm
				los_vs_rad_point = -los_vs_rad_point
			
			self.vec_azms[i] = vec_azm
			self.rad_azms[i] = radar_azm
			self.los_vs_rad_point[i] = los_vs_rad_point
			
		#calculate vector azimuths from kvector
		#(-ve vec_azms are east look, +ve vec_azms are west look)
		#then use los_v sign to denote direction of flow, -ve los_v = towards
		#radar
		self.vec_azms = np.array(self.kvecs)
		for i in range(len(self.vec_azms)):
 			look_direction = self.look[i]
 			vec_azm = self.vec_azms[i]
 			if look_direction == "E":
 					if vec_azm > 0:
						  self.vec_azms[i] -= 180
 			elif look_direction == "W":	
 					if vec_azm < 0:
						 self.vec_azms[i] += 180				
			
		#restrict azimuths between -90 and 90 degrees
		self.vec_azms = elliotools.sin9090(self.vec_azms)
		
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
		try:
			data_dict["vec_azms"] = self.vec_azms[indexes]
			data_dict["rad_azms"] = self.rad_azms[indexes]
			data_dict["los_vs_rad_point"] = self.los_vs_rad_point[indexes]
		except:
			0
		data_dict["snames"] = self.stations[indexes]
		
		return data_dict

	def get_time_data(self, time):
		
		"""
		Returns a dictionary of all data captured at a particular time
		
		Parameters
		----------
		
		time: string
			time to retrieve data for (YYYY/MM/DD hh:mm:ss)
		"""
		
		indexes = np.where(self.times == time)
		
		data_dict = dict()
		data_dict["mlats"] = self.mlats[indexes]
		data_dict["mcolats"] = self.mcolats[indexes]
		data_dict["mlons"] = self.mlons[indexes]
		data_dict["kvecs"] = self.kvecs[indexes]
		data_dict["los_vs"] = self.los_vs[indexes]
		data_dict["los_e"] = self.los_e[indexes]
		data_dict["times"] = self.times[indexes]
		data_dict["dtimes"] = self.dtimes[indexes]
		data_dict["stations"] = self.stations[indexes]
		try:
			data_dict["vec_azms"] = self.vec_azms[indexes]
			data_dict["rad_azms"] = self.rad_azms[indexes]
			data_dict["los_vs_rad_point"] = self.los_vs_rad_point[indexes]
		except:
			0
		data_dict["snames"] = self.stations[indexes]
		
		if self.is_map == True:
		
			data_dict["order"] = self.order[time]
			data_dict["N"] = self.N[time]
			data_dict["N1"] = self.N1[time]
			data_dict["N2"] = self.N2[time]
			data_dict["N3"] = self.N3[time]
			data_dict["imfBx"] = self.imfBx[time]
			data_dict["imfBy"] = self.imfBy[time]
			data_dict["imfBz"] = self.imfBz[time]
			data_dict["imfTilt"] = self.imfTilt[time]
			data_dict["boundary_mlats"] = self.boundary_mlats[time]
			data_dict["boundary_mlons"] = self.boundary_mlons[time]
			data_dict["mod_kvecs"] = self.mod_kvecs[time]
			data_dict["mod_mlats"] = self.mod_mlats[time]
			data_dict["mod_mcolats"] = self.mod_mcolats[time]
			data_dict["mod_mlons"] = self.mod_mlons[time]
			data_dict["mod_los_vs"] = self.mod_los_vs[time]
			data_dict["mod_times"] = time
			data_dict["mod_dtimes"] = self.mod_dtimes[time]
			data_dict["latmin"] = self.latmin[time]
			data_dict["lon_shft"] = self.lon_shft[time]
			data_dict["lat_shft"] = self.lat_shft[time]			 
		
		return data_dict
	
	def spherical_fit(self, time):
		
		return
	