#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 08:52:01 2020

@author: elliott
"""

import data_classes
import supermaps
import numpy as np

if __name__ == "__main__":

	bks = data_classes.LoadFitacf("bks", "2013/10/02/00.00.00", "2013/10/02/02.00.00")
	#obtain the data from all the retrieved files
	data = bks.fread(bks.file_list)
	data0 = data[0]
	
	beam_velocities, v_range = supermaps.beam_range_velocities(data0)
	for i in range(len(beam_velocities)):
		beam_velocities[i].cmap()