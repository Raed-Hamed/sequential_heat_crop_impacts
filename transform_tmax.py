#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 15:50:07 2022

@author: carmensteinmann
"""

from pathlib import Path
import xarray as xr
import numpy as np


#specify temperature and soil data path and harvest date file
path_tmax =  Path('/Users/carmensteinmann/Documents/CLIMADA/own_projects/sequential_heat_crop_impacts/data/cpc/tmax_365')


"""Read data files"""

#read tmax files
files = [f.name for f in path_tmax.iterdir() if f.is_file() and not f.name.startswith('.')]
files.sort()

#get data dimensions from one tmax file
ds_test = xr.open_dataset(Path(path_tmax, files[0]))
tmax_test = np.roll(ds_test.tmax.values, 360, axis=2) 
lat_tmax = ds_test.lat.values
lon_tmax = np.roll(ds_test.lon.values, 360) 
lon_tmax[0:360] =  lon_tmax[0:360]-360

for idx_file, file in enumerate(files):
    ds = xr.open_dataset(Path(path_tmax, file))
    tmax = np.roll(ds.tmax.values, 360, axis=2)
    output_path = Path('/Users/carmensteinmann/Documents/CLIMADA/own_projects/sequential_heat_crop_impacts/data/output/tmax180')
    ds_output = xr.Dataset(data_vars=dict(tmax=(["time", "lat", "lon"], tmax)),
                           coords=dict(
                               lat=(["lat"], lat_tmax),
                               lon=(["lon"], lon_tmax),
                               time=(["time"], ds_test.time.values)
                               )
                           )

    ds_output.to_netcdf(Path(output_path, file))
    ds_output.close()

