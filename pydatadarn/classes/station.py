import numpy as np
from pydatadarn.utils import tools
from pydatadarn.utils import coordinate_transformations as coords
import os
import datetime

#get path to superdarn rst
rst = os.environ["RSTPATH"]
hdw_path = rst + "/tables/superdarn/hdw/"

#generate arrays for the data
#for more information see hdw tables in supderdarn rst

class Station():
	
	"""
	A class used to access and store hardware data for superdarn radar stations. 
	Data is collected from the superdarn rst (required) hdw tables.
	"""
	
	def __init__(self, rad):
		
		"""
		Parameters
		----------
		
		rad: string
			3 letter station code for requested radar data (e.g. "ade")
		"""
		
		#find hdw file
		fname = "hdw.dat.{}".format(rad)
		file = hdw_path + fname
		
		#generate arrays for the data
		#for more information see hdw tables in supderdarn rst
		self.stat_id = np.array([], dtype = "int")
		self.years = np.array([], dtype = "int")
		self.secs = np.array([], dtype = "int")
		self.glats = np.array([], dtype = "float")
		self.glons = np.array([], dtype = "float")
		self.alts = np.array([], dtype = "float")
		self.boresight = np.array([], dtype = "float")
		self.beam_sep = np.array([], dtype = "float")
		self.vel_sign = np.array([], dtype = "float")
		self.an_rx_atten_db = np.array([], dtype = "float")
		self.tdiff = np.array([], dtype = "float")
		self.phase_sign = np.array([], dtype = "float")
		self.interf_offset_x = np.array([], dtype = "float")
		self.interf_offset_y = np.array([], dtype = "float")
		self.interf_offset_z = np.array([], dtype = "float")
		self.an_rx_rise_time = np.array([], dtype = "float")
		self.an_atten_stages = np.array([], dtype = "int")
		self.max_range_gates = np.array([], dtype = "int")
		self.max_beams = np.array([], dtype = "int")
		self.dtimes = np.array([])

		#open the file		
		with open(file) as fp:
			lines = fp.readlines()
			for line in lines:
				if line[0:5] == "# EOF":
					break
				elif line[0] == "#":
					continue
				else:
					#add data to arrays
					data = line.split()
					self.stat_id = np.append(self.stat_id, int(data[0]))
					
					#if final entry for year is a ridiculously large number
					#we can run into compatibility issues so set it to current
					#date
					year = int(data[1])
					sec = int(data[2])
					if year > datetime.datetime.now().year :
						year = datetime.datetime.now().year
						#also set the seconds to be the time now
						sec = datetime.datetime.now().timetuple().tm_yday*24*3600
					
					self.years = np.append(self.years, year)
					self.secs = np.append(self.secs, sec)
					self.glats = np.append(self.glats, float(data[3]))
					self.glons = np.append(self.glons, float(data[4]))
					self.alts = np.append(self.alts, float(data[5]))
					self.boresight = np.append(self.boresight, float(data[6]))
					self.beam_sep = np.append(self.beam_sep, float(data[7]))
					self.vel_sign = np.append(self.vel_sign, float(data[8]))
					self.an_rx_atten_db = np.append(self.an_rx_atten_db, float(data[9]))
					self.tdiff = np.append(self.tdiff, float(data[10]))
					self.phase_sign = np.append(self.phase_sign, float(data[11]))
					self.interf_offset_x = np.append(self.interf_offset_x, float(data[12]))
					self.interf_offset_y = np.append(self.interf_offset_y, float(data[13]))
					self.interf_offset_z = np.append(self.interf_offset_z, float(data[14]))
					self.an_rx_rise_time = np.append(self.an_rx_atten_db, float(data[15]))
					self.an_atten_stages = np.append(self.an_atten_stages, int(data[16]))
					self.max_range_gates = np.append(self.max_range_gates, int(data[17]))
					self.max_beams = np.append(self.max_beams, int(data[18]))
					#get datetimes
					dtime = datetime.datetime(year, 1, 1, 0, 0, 0) + datetime.timedelta(seconds=int(data[2]))
					self.dtimes = np.append(self.dtimes, dtime)
		
	def get_coords(self, dtime, aacgm=True):			
	
		"""
		Calculates and returns aacmgv2 coordinates for this radar from the
		specified time
		
		Parameters
		----------
		
		time: datetime object
			datetime object of format datetime.datetime(YYYY, MM, DD, hh, mm, ss)
			
		aacgm: bool
			if True, will return coordinates as aacgm. If false will return coordinates
			as geographic (default=True)
		"""
	
		#get first time in station metadata that is greater
		#than time of data requested. The correct data will
		#then be the metadata before that
		station_time_index = min(np.where(self.dtimes > dtime)[0])-1
		glat = self.glats[station_time_index]
		glon = self.glons[station_time_index]	
		alt = self.alts[station_time_index]

		if aacgm == True:
			#get mlat and mlon
			mlat, mlon = coords.geo_to_aacgm(glat, glon, dtime, alt)
			if isinstance(mlat, np.ndarray):
				mlat = mlat[0]
			if isinstance(mlon, np.ndarray):
				mlon = mlon[0]
			return mlat, mlon
		
		else:
			return glat, glon
	