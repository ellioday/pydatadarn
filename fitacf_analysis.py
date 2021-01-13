#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 07:12:25 2020

@author: elliott
"""

import data_classes

fitacf_data = data_classes.LoadFitacf("ade", "2013/10/02/09.00.00", "2013/10/02/10.00.00")
file_list = fitacf_data.file_list
fname = file_list[0]
data = fitacf_data.read_bz2(fname)

sample = data[307]