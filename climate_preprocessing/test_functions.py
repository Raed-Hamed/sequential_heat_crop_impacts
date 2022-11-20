#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 13:39:02 2022

@author: carmensteinmann
"""

"""TESTING"""
# year = 6
# lat = 107
# lon = 130
# test_kdd = kdd_summer[year,lat,lon]

# harvest_grid_cell = harvest_end_mean[lat, lon]
# day_summer_start = int(harvest_grid_cell-summer_start)
# day_summer_end = int(harvest_grid_cell-summer_end)
# kdd_summer_gridcell = tmax[day_summer_start:day_summer_end, lat, lon, :] - t_high_summer
# kdd_summer_gridcell[kdd_summer_gridcell<=0] = 0
# kdd_summer_gridcell[np.isnan(kdd_summer_gridcell)] = 0
# size_kdd = kdd_summer_gridcell.shape
# total_kdd = np.sum(kdd_summer_gridcell, axis=0)

# #testing gdd
# year = 0
# lat = 82
# lon = 140
# # year = 0
# # lat = 47
# # lon = 66
# np.max(gdd_spring[~np.isnan(gdd_spring)]) / gdd_max_spring
# np.max(gdd_summer[~np.isnan(gdd_summer)]) / gdd_max_summer
# max_gdd = np.max(gdd_summer[~np.isnan(gdd_summer)])
# median_gdd = np.median(gdd_summer[~np.isnan(gdd_summer)])
# #grid_cells = np.where(gdd_summer == max_gdd)
# grid_cells = np.where(gdd_summer >= median_gdd) and np.where(kdd_summer>=20)
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
                
#testing soil moisture
# year = 0
# lat = 82
# lon = 140
# harvest_grid_cell = harvest_end_mean[lat, lon]
# day_spring_start = int(harvest_grid_cell-spring_start)
# day_spring_end = int(harvest_grid_cell-spring_end)
# moisture_spring_gridcell = moisture[day_spring_start:day_spring_end, lat, lon, :]
# year0 = np.mean(moisture_spring_gridcell[:,3])
# mean_sm = np.mean(moisture_spring_gridcell, axis=0).T