import numpy as np

from pydatadarn.utils import tools
from pydatadarn.utils import coordinate_transformations as coords

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
		self.mlat, self.mlon = coords.geo_to_aacgm(self.lat, self.lon, time)
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
		
		mlt = coords.aacgm_to_mlt(self.mlon, time)
		
		if isinstance(mlt, np.ndarray):
			mlt = mlt[0]
		
		return mlt
		
	def set_boresight(self, boresight):
		
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
