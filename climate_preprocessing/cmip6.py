#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 13:41:04 2022

@author: carmensteinmann
"""
from pathlib import Path
import xarray as xr

path_cmip =  Path('/Volumes/Files/WCR/2022/sequential_heat/data/cmip6/mrsos_Lmon_IPSL-CM6A-LR_ssp126_r1i1p1f1_gr_20500116-20501216_v20190903.nc') 

 

"""Executing functions"""
#Read harvest dates, temperature and soil moisture data
ds_cmip = xr.open_dataset(path_cmip)