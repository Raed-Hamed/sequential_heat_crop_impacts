#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 15:50:07 2022

@author: carmensteinmann
"""

from pathlib import Path
import xarray as xr
import numpy as np

"""Tmax"""
#specify temperature and soil data path and harvest date file
path_tmax =  Path('/Users/carmensteinmann/Documents/CLIMADA/own_projects/sequential_heat_crop_impacts/data/cpc/tmax_365')
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
                                time=(["time"], ds.time.values)
                                )
                            )

    ds_output.to_netcdf(Path(output_path, file))
    ds_output.close()

"""Tmin"""
path_tmin =  Path('/Users/carmensteinmann/Documents/CLIMADA/own_projects/sequential_heat_crop_impacts/data/cpc/tmin_365')
#read tmax files
files = [f.name for f in path_tmin.iterdir() if f.is_file() and not f.name.startswith('.')]
files.sort()
#get data dimensions from one tmax file
ds_test = xr.open_dataset(Path(path_tmin, files[0]))
tmin_test = np.roll(ds_test.tmin.values, 360, axis=2) 
lat_tmin = ds_test.lat.values
lon_tmin = np.roll(ds_test.lon.values, 360) 
lon_tmin[0:360] =  lon_tmin[0:360]-360

for idx_file, file in enumerate(files):
    ds = xr.open_dataset(Path(path_tmin, file))
    tmin = np.roll(ds.tmin.values, 360, axis=2)
    output_path = Path('/Users/carmensteinmann/Documents/CLIMADA/own_projects/sequential_heat_crop_impacts/data/output/tmin180')
    ds_output = xr.Dataset(data_vars=dict(tmin=(["time", "lat", "lon"], tmin)),
                            coords=dict(
                                lat=(["lat"], lat_tmin),
                                lon=(["lon"], lon_tmin),
                                time=(["time"], ds.time.values)
                                )
                            )

    ds_output.to_netcdf(Path(output_path, file))
    ds_output.close()


"""Plotting to check shifting coordinates"""
# import climada.util.plot as u_plot
# def plot_tmax(t, lat, lon):
#     lon1, lat1 = np.meshgrid(lon, lat)
#     coords = np.append(lat1.reshape(360*720,1), lon1.reshape(360*720,1), axis=1)

#     col_name = 't max'
#     l_title = 'T max global'
#     u_plot.geo_im_from_array(t[0,:,:].reshape(360*720,1), coords, col_name,
#                                     l_title)

# plot_tmax(tmax_test, lat_tmax, lon_tmax)
# plot_tmax(ds_test.tmax.values, ds_test.lat.values, ds_test.lon.values)

# ds_test2 = xr.open_dataset(Path(output_path, files[10]))
# plot_tmax(ds_test2.tmax.values, ds_test2.lat.values, ds_test2.lon.values)

