#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 02:02:41 2021

@author: elliott
"""

import numpy as np
import tools
import data_classes
import plotting

#create simple ion flow of v = 400km/s along a single line of latitude
v = 400
mcolat = 37.5

#create a sample radar with measurements taken at beam azimuths
beam_azms = np.array([-10, -20, -30, -40])

#calculate what velocities should be, v = vmax at theta = 90deg, v = vmin at theta = 0deg
beam_velocities = v*np.sin(np.deg2rad(beam_azms))

#now test the line of sight fit
beam_w = plotting.los_fit(beam_azms, beam_velocities, mcolat_range = [0], time="beam azimuths", plot=True)

#lets check vector azimuths
radar_mcolat = 55
radar_mlon = 20

#use law of sines to calculate vector azimuth
vec_azms = np.rad2deg(np.arcsin((radar_mcolat/mcolat)*np.sin(np.deg2rad(beam_azms))))

vec_w = plotting.los_fit(vec_azms, beam_velocities, mcolat_range = [0], time = "vector azimuths", plot=True)