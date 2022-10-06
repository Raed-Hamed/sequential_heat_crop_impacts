#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 16:19:10 2022

@author: carmensteinmann
"""

from pathlib import Path
import xarray as xr
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.transforms as transforms


def read_t_data(path, variable, grid_cells):
    #read tmax files
    files = [f.name for f in path.iterdir() if f.is_file() and not f.name.startswith('.')]
    files.sort()
    #get data dimensions from one tmax file
    ds_test = xr.open_dataset(Path(path, files[0]))
    lat_tmax = ds_test.lat.values
    lon_tmax = ds_test.lon.values
    
    #load temperature data for all years
    nr_years_tmax= len(files)
    tmax = np.empty(np.append(ds_test[variable].values.shape, nr_years_tmax))
    time = []
    for idx_file, file in enumerate(files):
        ds = xr.open_dataset(Path(path, file))
        tmax[:, grid_cells[0],grid_cells[1], idx_file] = ds[variable].values[:,grid_cells[0], grid_cells[1]]
        time.append(str(pd.to_datetime(ds.time.values[0]).year))
        ds.close()
        
    return tmax, lat_tmax, lon_tmax, time


def kdd_gdd_per_gridcell(grid_cells, harvest_end_mean, tmax, tmin, thr_spring, thr_summer, spring_dates, summer_dates):
    """Compute kdd for summer and spring for all grid cells with a harvest date"""
    t_base_spring, t_optimum_spring, t_high_spring, gdd_max_spring = thr_spring
    t_base_summer, t_optimum_summer, t_high_summer, gdd_max_summer = thr_summer
    spring_start, spring_end = spring_dates
    summer_start, summer_end = summer_dates
    
    nr_years = tmax.shape[3]
    nr_lon = tmax.shape[2]
    nr_lat = tmax.shape[1]
    
    shape_kdd = (nr_years, nr_lat, nr_lon)
    
    kdd_spring = np.empty(shape_kdd)
    kdd_spring[:,:,:] = np.nan
    kdd_summer = np.empty(shape_kdd)
    kdd_summer[:,:,:] = np.nan
    
    gdd_spring = np.empty(shape_kdd)
    gdd_spring[:,:,:] = np.nan
    gdd_summer = np.empty(shape_kdd)
    gdd_summer[:,:,:] = np.nan
    
    # nr_years_moisture = 30
    # moisture_spring = np.zeros((nr_years_moisture, 360, 720))
    # moisture_summer = np.zeros((nr_years_moisture, 360, 720))
    
    for idx, _ in enumerate(grid_cells[0]):
        #grid cell
        lat = grid_cells[0][idx]
        lon = grid_cells[1][idx]
        harvest_grid_cell = harvest_end_mean[lat, lon]
    
        #spring days and temperatures for current grid cell
        day_spring_start = int(harvest_grid_cell-spring_start)
        day_spring_end = int(harvest_grid_cell-spring_end)
        tmax_spring = tmax[day_spring_start:day_spring_end, lat, lon, :]
        tmin_spring = tmin[day_spring_start:day_spring_end, lat, lon, :]
        tmean_spring = (tmax_spring + tmin_spring)/2 
        
        #summer days and temperatures for current grid cell
        day_summer_start = int(harvest_grid_cell-summer_start)
        day_summer_end = int(harvest_grid_cell-summer_end)
        tmax_summer = tmax[day_summer_start:day_summer_end, lat, lon, :]
        tmin_summer = tmin[day_summer_start:day_summer_end, lat, lon, :]
        tmean_summer = (tmax_summer + tmin_summer)/2 
        
        #spring and summer kdd and gdd
        kdd_spring[:, lat, lon] = compute_kdd(tmax_spring, t_high_spring)
        kdd_summer[:, lat, lon] = compute_kdd(tmax_summer, t_high_summer)
        gdd_spring[:, lat, lon] = compute_gdd(tmean_spring, t_base_spring, gdd_max_spring)
        gdd_summer[:, lat, lon] = compute_gdd(tmean_summer, t_base_summer, gdd_max_summer)
        
        # #soil moisture
        # moisture_spring_gridcell = moisture[day_spring_start:day_spring_end, lat, lon, :]
        # moisture_spring[:,lat, lon] = np.mean(moisture_spring_gridcell, axis=0)
        # moisture_summer_gridcell = moisture[day_summer_start:day_summer_end, lat, lon, :]
        # moisture_summer[:,lat, lon] = np.mean(moisture_summer_gridcell, axis=0)
        
    return kdd_spring, kdd_summer, gdd_spring, gdd_summer


def compute_kdd(tmax_season, t_high_season):
    kdd_daily = tmax_season - t_high_season
    kdd_daily[kdd_daily<=0] = 0
    #kdd_summer_gridcell[np.isnan(kdd_summer_gridcell)] = 0
    kdd_annual = np.sum(kdd_daily, axis=0)
    
    return kdd_annual

def compute_gdd(tmean_season, t_base_season, gdd_max_season):
    gdd_daily= tmean_season - t_base_season
    gdd_daily[gdd_daily<=0] = 0
    gdd_daily[gdd_daily>=gdd_max_season] = gdd_max_season
    gdd_annual = np.sum(gdd_daily, axis=0)
    
    return gdd_annual


def save_outputs(path_output, filename_output, sum_kdd_spring, sum_kdd_summer, gdd_spring, gdd_summer, 
                 time_tmax, lat_tmax, lon_tmax):                                                    
    """Save output"""
    ds_output = xr.Dataset(data_vars=dict(kdd_spring=(["time", "lat", "lon"], sum_kdd_spring),
                                   kdd_summer=(["time", "lat", "lon"], sum_kdd_summer),
                                   gdd_spring=(["time", "lat", "lon"], gdd_spring), 
                                   gdd_summer=(["time", "lat", "lon"], gdd_summer)),
                           coords=dict(
                               time=(["time"], time_tmax),
                               lat=(["lat"], lat_tmax),
                               lon=(["lon"], lon_tmax),
                               )
                           )
    
    ds_output.to_netcdf(Path(path_output, filename_output))
    ds_output.close()


def plot_gdd_kdd(year, lat, lon, tmax, tmin, harvest_end_mean, thr_spring, thr_summer, spring_dates, 
                 summer_dates, season, filename):
    """PLOTTING gdd and kdd"""
    year = 0
    lat = 82
    lon = 140
    
    t_base_spring, t_optimum_spring, t_high_spring, gdd_max_spring = thr_spring
    t_base_summer, t_optimum_summer, t_high_summer, gdd_max_summer = thr_summer
    spring_start, spring_end = spring_dates
    summer_start, summer_end = summer_dates
    
    harvest_grid_cell = harvest_end_mean[lat, lon]
    # day_spring_start = int(harvest_grid_cell-spring_start)
    # day_spring_end = int(harvest_grid_cell-spring_end)
    day_summer_start = int(harvest_grid_cell-summer_start)
    day_summer_end = int(harvest_grid_cell-summer_end)
    tmax_summer = tmax[day_summer_start:day_summer_end, lat, lon, year]
    tmean_summer = (tmax[day_summer_start:day_summer_end, lat, lon, year] + 
                    tmin[day_summer_start:day_summer_end, lat, lon, year])/2 
    x_axis = np.arange(day_summer_start, day_summer_end)
    y4 = np.minimum(np.ones(30)*t_optimum_summer, tmean_summer)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hlines(t_high_summer, day_summer_start, day_summer_end, color='k', linestyles = 'dashed')
    ax.hlines(t_optimum_summer, day_summer_start, day_summer_end, color='k', linestyles = 'dashed')
    ax.hlines(t_base_summer, day_summer_start, day_summer_end, color='k', linestyles = 'dashed')
    ax.plot(np.arange(day_summer_start, day_summer_end), tmax_summer, color='#e6550d', label='$T_{max, gridcell}$')
    ax.plot(np.arange(day_summer_start, day_summer_end), tmean_summer, color='#fdae6b', label='$T_{mean, gridcell}$')
    ax.fill_between(x_axis, tmax_summer, t_high_summer, where=tmax_summer>=t_high_summer, color='k', alpha=0.3, label='$KDD_{summer}$', interpolate=True)
    ax.fill_between(x_axis, y4, t_base_summer, color='k', alpha=0.1, label='$GDD_{summer}$', interpolate=True)
    plt.ylabel('Temperature [°C]')
    plt.xlabel('Day of the year')
    plt.xlim(day_summer_start, day_summer_end-1)
    
    trans = transforms.blended_transform_factory(
        ax.get_yticklabels()[0].get_transform(), ax.transData)
    ax.text(1, t_high_summer+0.8, '$T_{high}$' , color="k", transform=trans, 
            ha="right", va="center")
    ax.text(1, t_optimum_summer+0.8, '$T_{optimum}$' , color="k", transform=trans, 
            ha="right", va="center")
    ax.text(1, t_base_summer +0.8, '$T_{base}$' , color="k", transform=trans, 
            ha="right", va="center")
    plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
    plt.savefig('kdd_gdd_vis.pdf', format='pdf', bbox_inches='tight')



"""TESTING"""
# test_kdd = sum_kdd_summer[6,107,130]
# year = 6
# lat = 107
# lon = 130
# harvest_grid_cell = harvest_end_mean[lat, lon]
# day_summer_start = int(harvest_grid_cell-summer_start)
# day_summer_end = int(harvest_grid_cell-summer_end)
# kdd_summer_gridcell = tmax[day_summer_start:day_summer_end, lat, lon, :] - t_high_summer
# kdd_summer_gridcell[kdd_summer_gridcell<=0] = 0
# #kdd_summer_gridcell[np.isnan(kdd_summer_gridcell)] = 0
# size_kdd = kdd_summer_gridcell.shape
# total_kdd = np.sum(kdd_summer_gridcell, axis=0)

# #testing gdd
# np.max(gdd_spring[~np.isnan(gdd_spring)]) / gdd_max_spring
# np.max(gdd_summer[~np.isnan(gdd_summer)]) / gdd_max_summer
# max_gdd = np.max(gdd_summer[~np.isnan(gdd_summer)])
# median_gdd = np.median(gdd_summer[~np.isnan(gdd_summer)])
# #grid_cells = np.where(gdd_summer == max_gdd)
# grid_cells = np.where(gdd_summer >= median_gdd) and np.where(sum_kdd_summer>=20)
# year = 0
# lat = 82
# lon = 140
# # year = 0
# # lat = 47
# # lon = 66
# harvest_grid_cell = harvest_end_mean[lat, lon]
# day_spring_start = int(harvest_grid_cell-spring_start)
# day_spring_end = int(harvest_grid_cell-spring_end)
# day_summer_start = int(harvest_grid_cell-summer_start)
# day_summer_end = int(harvest_grid_cell-summer_end)
# tmax_gridcell = tmax[day_spring_start:day_spring_end, lat, lon, :]
# tmin_gridcell = tmin[day_spring_start:day_spring_end, lat, lon, :]
# gdd_spring_gridcell = (tmax[day_spring_start:day_spring_end, lat, lon, :] + 
#                 tmin[day_spring_start:day_spring_end, lat, lon, :])/2 - t_base_spring
# gdd_spring_gridcell = (tmax_gridcell + tmin_gridcell)/2 - t_base_spring
# gdd_spring_gridcell[gdd_spring_gridcell<=0] = 0
# gdd_spring_gridcell[gdd_spring_gridcell>=gdd_max_spring] = gdd_max_spring
# gdd_spring[:, lat, lon] = np.sum(gdd_spring_gridcell, axis=0)

# gdd_summer_gridcell = (tmax[day_summer_start:day_summer_end, lat, lon, :] + 
#                 tmin[day_summer_start:day_summer_end, lat, lon, :])/2 - t_base_summer
# gdd_summer_gridcell[gdd_summer_gridcell<=0] = 0
# gdd_summer_gridcell[gdd_summer_gridcell>=gdd_max_summer] = gdd_max_summer
# gdd_summer[:, lat, lon] = np.sum(gdd_summer_gridcell, axis=0)

# kdd_max_year = tmax_gridcell[:, int(year)]
# gdd_spring_year = gdd_spring_gridcell[:, int(year)]



"""PLOTTING ONLY KDD"""
# year_max, lat_max, lon_max = np.where(sum_kdd_summer == np.max(sum_kdd_summer[~np.isnan(sum_kdd_summer)]))
# harvest_grid_cell = harvest_end_mean[int(lat_max), int(lon_max)]
# day_summer_start = int(harvest_grid_cell-summer_start)
# day_summer_end = int(harvest_grid_cell-summer_end)
# kdd_max = tmax[day_summer_start:day_summer_end, lat_max, lon_max, int(year_max)]


# x_axis = np.arange(day_summer_start, day_summer_end)
# plt.hlines(t_high_summer, day_summer_start, day_summer_end, color='#fdae6b', label='$T_{high, summer}$')
# plt.plot(np.arange(day_summer_start, day_summer_end), kdd_max, color='#e6550d', label='$T_{max, gridcell}$')
# plt.fill_between(x_axis, kdd_max.flatten(), np.ones(30)*t_high_summer, color='#e6550d', alpha=0.1, label='$KDD_{summer}$')
# plt.ylabel('Temperature [°C]')
# plt.xlabel('Day of the year')
# plt.legend(loc='upper right')
# plt.xlim(day_summer_start, day_summer_end-1)
# plt.savefig('kdd_vis.png', format='png', bbox_inches='tight')


# import climada.util.plot as u_plot
# def plot_tmax(t, lat, lon):
#     lon1, lat1 = np.meshgrid(lon, lat)
#     coords = np.append(lat1.reshape(360*720,1), lon1.reshape(360*720,1), axis=1)

#     col_name = 't max'
#     l_title = 'T max global'
#     u_plot.geo_im_from_array(t[0,:,:].reshape(360*720,1), coords, col_name,
#                                     l_title)

# plot_tmax(sum_kdd_summer, lat_tmax, lon_tmax)
