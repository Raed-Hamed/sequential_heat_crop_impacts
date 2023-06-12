#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 10:45:18 2023

@author: carmenst
"""

import os
import xarray as xr
import numpy as np
from glob import iglob
from pathlib import Path
import geopandas as gpd
import xesmf as xe
import pandas as pd

def create_mean(rootdir_glob, nr_years, old_time_period):

    file_list = [f for f in iglob(rootdir_glob, recursive=True) if os.path.isfile(f)]


    for file in file_list:
        new_filename = file.replace("input", "output")
        new_filename = new_filename.replace(old_time_period, "mean")
        
        ds_tmax= xr.open_dataset(file)
        tmax = ds_tmax.tasmax.values[-nr_years:, :, :]
        mean_tmax = np.mean(tmax, axis=0)
        

        ds_output = xr.Dataset(data_vars=dict(tmax_mean=(["lat", "lon"], mean_tmax),
                                       ),
                               coords=dict(                           
                                       lat=(["lat"], ds_tmax.lat.values),
                                       lon=(["lon"], ds_tmax.lon.values),
                                       )    
                               )
        ds_output.attrs = ds_tmax.attrs
        
        Path(os.path.dirname(new_filename)).mkdir(parents=True, exist_ok=True)
        
        ds_output.to_netcdf(new_filename)
        ds_output.close()


def compute_difference(path_future_output):
    file_list = [f for f in iglob(path_future_output, recursive=True) if os.path.isfile(f)]
    
    # nr_runs = len(file_list)
    
    # ds_mean_fut = xr.open_dataset(file_list[0])
    # size = ds_mean_fut.tmax_mean.values.shape
    
    # multi_model_diff = np.empty((nr_runs, size[0], size[1]))
    
    for idx_file, fut_file in enumerate(file_list):
        
        new_filename = fut_file.replace("future", "future-historic")
        new_filename1 = new_filename.replace("mean", "difference")
        
        if not os.path.isfile(new_filename1):
            Path(os.path.dirname(new_filename)).mkdir(parents=True, exist_ok=True)
        
            his_file0 = fut_file.replace("future", "historic")
            his_file1 = his_file0.replace("ssp245/", "")
            his_file2 = his_file1.replace("ssp370/", "")
            his_file3 = his_file2.replace("ssp119/", "")
            his_file4 = his_file3.replace("ssp126/", "")
            
    
            
            directory0, directory1, directory2, directory3, directory4, str1, str2, str3, _ , _, _  = his_file4.split('_')
            
            his_file = directory0+'_'+directory1+'_'+directory2+ '_'+directory3+'_'+directory4+'_'+str1+'_'+str2+'_'+str3+'_'+'hist_r1i1p1f1_mean.nc'
            
            ds_mean_fut = xr.open_dataset(fut_file)
            ds_mean_his = xr.open_dataset(his_file)
            
            mean_fut = ds_mean_fut.tmax_mean.values
            mean_his = ds_mean_his.tmax_mean.values
            
            difference = mean_fut-mean_his
            
            
            ds_output = xr.Dataset(data_vars=dict(diff=(["lat", "lon"], difference),
                                           ),
                                   coords=dict(                           
                                           lat=(["lat"], ds_mean_fut.lat.values),
                                           lon=(["lon"], ds_mean_fut.lon.values),
                                           )    
                                   )
            
            ds_output.attrs = ds_mean_his.attrs
            
    
            ds_output.to_netcdf(new_filename)
            ds_output.close()
        
        # multi_model_diff[idx_file, :, :] = 
        
    
    
    # multi_model_mean = np.mean(multi_model_diff, axis=0)
    
    # return multi_model_mean


# def compute_mean_per_county(path_input):
    
#     us_shp_file = '/Users/carmenst/Documents/Polybox/WCR/Conferences_summerschools/2022_10_Como/Project/Sequential_heat_crops/Data/shapefiles/input/gadm36_USA_2.shp'
#     crops_shapefile = gpd.read_file(us_shp_file)
 
    
#     dir_list = [f for f in iglob(path_input, recursive=True) if os.path.isdir(f)][2:]
    
#     for directory in dir_list:
#         file_list = [f for f in iglob(directory+'/*', recursive=True) if os.path.isfile(f)]
    
#         multi_models = np.zeros((crops_shapefile.shape[0], len(file_list)))
        
#         for file_idx, file in enumerate(file_list):
            
#             # open file containing diffence in mean per model
#             ds = xr.open_dataset(file)
    
#             _, _, _, _, _, _, _, _, _, _, _, _, _, _, reference, percentile, f_increase, ssp,filename = file.split('/')
#             _, _, _, model, _, _, nr_percen = filename.split('_')
            
#             mean_per_county = fct_mean_per_county(ds, model, 'diff', crops_shapefile)            
#             multi_models[:, file_idx] = mean_per_county
        
#         _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, ssp, season = directory.split('/')
#         crops_shapefile[ssp+'_'+season] = np.nanmean(multi_models, axis=1)
    
#     return crops_shapefile



def compute_diff_per_model(path_input):
    
    us_shp_file = '/Users/carmenst/Documents/Polybox/WCR/Conferences_summerschools/2022_10_Como/Project/Sequential_heat_crops/Data/shapefiles/input/gadm36_USA_2.shp'
    crops_shapefile = gpd.read_file(us_shp_file)
    
    dir_list = [f for f in iglob(path_input+'/**/*', recursive=True) if os.path.isdir(f)]
    
    for directory in dir_list:
        file_list = [f for f in iglob(directory+'/*', recursive=True) if os.path.isfile(f)]
    
        # model_dif = np.zeros((crops_shapefile.shape[0]))
        multi_models = np.zeros((crops_shapefile.shape[0], len(file_list)))
        
        _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, ssp, season = directory.split('/')
        
        for file_idx, file in enumerate(file_list):
            
            #get model name for regridding
            _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, filename = file.split('/')
            _, _, model, _, _, _ = filename.split('_')
            
            # open file containing diffence in mean per model
            ds = xr.open_dataset(file)

            mean_per_county = fct_mean_per_county(ds, model, 'diff', crops_shapefile)            
            multi_models[:, file_idx] = mean_per_county
            crops_shapefile[ssp[3:]+'_'+season[:2]+'_M'+str(file_idx)] = mean_per_county
    
        crops_shapefile[ssp[3:]+'_'+season[:2]+'_mean'] = np.nanmean(multi_models, axis=1)
            
    return crops_shapefile


def compute_frequency(path_input, nr_years, percentile, reference, ssps):
    
    #reference can be the given percentile for the hist or the fut time span
    
    dir_list = [f for f in iglob(path_input, recursive=True) if os.path.isdir(f)]
    
    for directory in dir_list:
        file_list = [f for f in iglob(directory+'/*', recursive=True) if os.path.isfile(f)]
    
        # multi_models = np.zeros((crops_shapefile.shape[0], len(file_list)))
        
        for file_idx, file in enumerate(file_list):
            
            # open file containing diffence in mean per model
            ds_his = xr.open_dataset(file)
            tmax_his = ds_his.tasmax.values[-nr_years:, :, :]
            percentile_his = np.percentile(tmax_his, percentile, axis=0)            
            
            diff_his_perc = tmax_his[:] - percentile_his
            diff_his_perc[diff_his_perc<=0] = 0
            diff_his_perc[diff_his_perc>0] = 1
            # greater_than_perc = np.sum(diff_his_perc, axis=0)
            
            new_file = file.replace("historic/input", "frequencies/reference_"+reference+"/percentile_"+str(percentile)+"/exceedances/hist")
            new_file1 = new_file.replace("1961-2014", str(percentile))
            
            ds_output = xr.Dataset(data_vars=dict(frequency=(["time", "lat", "lon"], diff_his_perc),
                                           ),
                                   coords=dict(                           
                                           time=(["time"], ds_his.time.values[-nr_years:]),
                                           lat=(["lat"], ds_his.lat.values),
                                           lon=(["lon"], ds_his.lon.values),
                                           )    
                                   )
            
            ds_output.attrs = ds_his.attrs
            
            Path(os.path.dirname(new_file1)).mkdir(parents=True, exist_ok=True)
            
            ds_output.to_netcdf(new_file1)
            ds_output.close()
                
            #get corresponding future file:
            for ssp in ssps:
                fut_file0 = file.replace("historic/input", "future/input/"+ssp)
                fut_file1 = fut_file0.replace("1961-2014", "2015-2100")
                fut_file2 = fut_file1.replace("hist", ssp)
                
                try:
                   if os.path.isfile(fut_file2):
                       ds_fut = xr.open_dataset(fut_file2)
                       tmax_fut = ds_fut.tasmax.values[-nr_years:, :, :]
                       if reference == 'historic':
                           diff_fut_perc = tmax_fut[:] - percentile_his
                           diff_fut_perc[diff_fut_perc<=0] = 0
                           diff_fut_perc[diff_fut_perc>0] = 1
                       elif reference == 'future':
                           percentile_fut = np.percentile(tmax_fut, percentile, axis=0) 
                           diff_fut_perc = tmax_fut[:] - percentile_fut
                           diff_fut_perc[diff_fut_perc<=0] = 0
                           diff_fut_perc[diff_fut_perc>0] = 1
                           
                           
                       new_file = file.replace("historic/input", "frequencies/reference_"+reference+"/percentile_"+str(percentile)+"/exceedances/"+ssp)
                       new_file1 = new_file.replace("1961-2014", str(percentile))
                       new_file2 = new_file1.replace("/hist", "/"+ssp)
                    
                       ds_output = xr.Dataset(data_vars=dict(frequency=(["time", "lat", "lon"], diff_fut_perc),
                                                    ),
                                            coords=dict(                           
                                                    time=(["time"], ds_fut.time.values[-nr_years:]),
                                                    lat=(["lat"], ds_fut.lat.values),
                                                    lon=(["lon"], ds_fut.lon.values),
                                                    )    
                                            )
                    
                       ds_output.attrs = ds_fut.attrs
                    
                       Path(os.path.dirname(new_file2)).mkdir(parents=True, exist_ok=True)
                    
                       ds_output.to_netcdf(new_file2)
                       ds_output.close()
                           
                       
                except Exception:
                   print('Does not exist')
                
                

def compute_frequency_change(path_input):
    
    file_list = [f for f in iglob(path_input+'/*/spring/*', recursive=True) if os.path.isfile(f)]
        
    for file in file_list:
        ds_spring = xr.open_dataset(file)
        f_spring = ds_spring.frequency.values
        
        file_su0 = file.replace("spring", "summer")
        file_su = file_su0.replace("sp_", "su_")
        ds_su = xr.open_dataset(file_su)
        f_su = ds_su.frequency.values
        
        coincide = np.where(f_spring+f_su ==2)
        coincide_matrix = np.zeros(f_spring.shape)
        coincide_matrix[coincide] = 1
        nr_coincides = np.sum(coincide_matrix, axis=0)
        
        
        new_file = file.replace("exceedances", "frequency")
        new_file1 = new_file.replace("spring/sp_", "coincide_")
        
        dict_var = dict(frequency=(["lat", "lon"], nr_coincides))
        create_ds_2D(ds_spring.lat.values, ds_spring.lon.values, dict_var, new_file1)
        
        # ds_output = xr.Dataset(data_vars=dict(frequency=(["lat", "lon"], nr_coincides),
        #                                 ),
        #                         coords=dict(                           
        #                                 lat=(["lat"], ds_spring.lat.values),
        #                                 lon=(["lon"], ds_spring.lon.values),
        #                                 )    
        #                         )
        
        # ds_output.attrs = ds_spring.attrs
        
        # Path(os.path.dirname(new_file1)).mkdir(parents=True, exist_ok=True)
        
        # ds_output.to_netcdf(new_file1)
        # ds_output.close()
        

def comparison_fut_his(path, ssp):
    
    files_his = [f for f in iglob(path+'/'+ssp+'/*', recursive=True) if os.path.isfile(f)]

    for file in files_his:
        ds_fut = xr.open_dataset(file)
        co_fut = ds_fut.frequency.values
        
        file_his = file.replace(ssp, "hist")
        ds_his = xr.open_dataset(file_his)
        co_his = ds_his.frequency.values
        
        f_increase = (co_fut/co_his)-1
        
        new_file = file.replace("frequency", "frequency increase")
        new_file1 = new_file.replace("coincide_", "f_change_")
        
        ds_output = xr.Dataset(data_vars=dict(f_change=(["lat", "lon"], f_increase),
                                        ),
                                coords=dict(                           
                                        lat=(["lat"], ds_his.lat.values),
                                        lon=(["lon"], ds_his.lon.values),
                                        )    
                                )
        
        # ds_output.attrs = ds_fut.attrs
        
        Path(os.path.dirname(new_file1)).mkdir(parents=True, exist_ok=True)
        
        ds_output.to_netcdf(new_file1)
        ds_output.close()
            

def compute_f_mean(path):
        
    us_shp_file = '/Users/carmenst/Documents/Polybox/WCR/Conferences_summerschools/2022_10_Como/Project/Sequential_heat_crops/Data/shapefiles/input/gadm36_USA_2.shp'
    crops_shapefile = gpd.read_file(us_shp_file)
    
    files_his = [f for f in iglob(path+'/*', recursive=True) if os.path.isfile(f)]
    multi_models = np.zeros((crops_shapefile.shape[0], len(files_his)))
    
    
    for idx, file in enumerate(files_his):
        ds = xr.open_dataset(file)
        # ds['f_change'].values[ds['f_change'].values>=100] = 100000
        ds['f_change'].values[ds['f_change'].values>=100] = np.nan
        
        _, _, _, _, _, _, _, _, _, _, _, _, _, _, reference, percentile, f_increase, ssp,filename = file.split('/')
        _, _, _, model, _, _, nr_percen = filename.split('_')
        
        mean_counties = fct_mean_per_county(ds, model, 'f_change', crops_shapefile)
        
        multi_models[:, idx] = mean_counties  

    new_file1 = file.replace(reference+'/'+percentile+'/'+f_increase+'/'+ssp,'/result')
    new_file2 = new_file1.replace(filename, 'f_change_relative_'+ssp+'_'+reference+'_'+percentile+'.shp')
    
    # idx_unprece = np.where(multi_models == 100000)
    # mask_unprece = np.zeros(multi_models.shape)
    # mask_unprece[idx_unprece] = 1
    # perc_unprece = (np.sum(mask_unprece, axis=1))/len(files_his)
    
    increase_agreement = np.zeros(multi_models.shape)
    increase_agreement[np.where(multi_models >0)] = 1
    
    # remove nan values
    # multi_models[idx_unprece] = np.nan
    
    # multi model mean (disregarding nan values) and percentage of models predicting an increase in co-occurence
    crops_shapefile['mean_incre'] = np.nanmean(multi_models, axis=1)
    crops_shapefile['%M_incre'] = np.sum(increase_agreement, axis=1)/len(files_his)
    
    # percentage of models agreeing in amplitude of frequency change
    oom_of_mean = np.zeros(multi_models.shape)
    std = np.nanstd(multi_models, axis=1)
    larger_than_min = np.where(multi_models >(crops_shapefile['mean_incre'].values-std)[:,None])
    smaller_than_max = np.where(multi_models <(crops_shapefile['mean_incre'].values+std)[:,None])
    oom_of_mean[larger_than_min and smaller_than_max] = 1
    percen_m_increase = (np.sum(oom_of_mean, axis=1))/len(files_his)
    crops_shapefile['%M_ampli'] = percen_m_increase
    
    # save percentage of models predicting unprecedented hot-hot extremes
    # crops_shapefile['%M_unpre'] = perc_unprece
    
    crops_shapefile.to_file(new_file2)
    

def fct_mean_per_county(ds, model, var_name, crops_shapefile):
    ds_out = regrid_cmip6(ds, model, var_name)
    gdf = create_gdf_from_ds(ds_out, var_name)            
    points_polys = gpd.sjoin(gdf, crops_shapefile, how="left")
    stats_pt = points_polys.groupby('index_right')[0].agg(['mean']) #nanmean
    result = pd.merge(crops_shapefile, stats_pt , left_index=True, right_index=True, how='outer')
    mean_per_county = np.asarray(result['mean'].values)
    
    return mean_per_county

def regrid_cmip6(ds, model, var_name):
    
    regridder_dir = '/Users/carmenst/Documents/Polybox/WCR/Conferences_summerschools/2022_10_Como/Project/Sequential_heat_crops/Data/Climate/regridder'
    regridder_file = regridder_dir+'/regridder_'+model+'.nc'
    
    ds_lat = ds.lat.values
    ds_lon = ds.lon.values-360
    
    ds_in = xr.Dataset(
        {
            "lat": (["lat"], ds_lat, {"units": "degrees_north"}),
            "lon": (["lon"], ds_lon, {"units": "degrees_east"}),
        }
    )
    
    ds_out = xr.Dataset(
        {
            "lat": (["lat"], np.arange(np.min(ds_lat), np.max(ds_lat), 0.1), {"units": "degrees_north"}),
            "lon": (["lon"], np.arange(np.min(ds_lon), np.max(ds_lon), 0.1), {"units": "degrees_east"}),
        }
    )
    
    
    if os.path.isfile(regridder_file):
        regridder = xe.Regridder(ds_in, ds_out, 'conservative', weights=regridder_file)
    else:
        regridder = xe.Regridder(ds_in, ds_out, "conservative", ignore_degenerate=True)
        regridder.to_netcdf(regridder_file)        
    
    regridded_ds = regridder(ds[var_name].values)    
    ds_out[var_name] = (('lat', 'lon'), regridded_ds)
    
    return ds_out

def create_gdf_from_ds(ds, var_name):
    lon1, lat1 = np.meshgrid(ds.lon.values, ds.lat.values)
    geometry = gpd.points_from_xy(lon1.flatten(), lat1.flatten())
    gdf = gpd.GeoDataFrame(ds[var_name].values.flatten(), crs="EPSG:4326", geometry=geometry)
    
    return gdf


def create_ds_2D(lat, lon, dict_var, file_path):
    ds_output = xr.Dataset(data_vars=dict_var, coords=dict(                           
                                    lat=(["lat"], lat),
                                    lon=(["lon"], lon),
                                    )    
                            )
    
    Path(os.path.dirname(file_path)).mkdir(parents=True, exist_ok=True)
    
    ds_output.to_netcdf(file_path)
    ds_output.close()
    

"""Compute mean per time period and grid cell for historic and future time period"""
# # future
# path_future = '/Users/carmenst/Documents/Polybox/WCR/Conferences_summerschools/2022_10_Como/Project/Sequential_heat_crops/Data/Climate/CMIP6/future/input/**/*' 
# nr_years = 40
# create_mean(path_future, nr_years, "2015-2100")


# # historic
# path_historic = '/Users/carmenst/Documents/CLIMADA/own_projects/sequential_heat_crop_impacts/CMIP6_DATA/historic/input/**/*' 
# # '/Users/carmenst/Documents/CLIMADA/own_projects/sequential_heat_crop_impacts/CMIP6_DATA/output'

# create_mean(path_historic, nr_years, "1961-2014")


# # """Compute multi model mean"""
# path_future_output = '/Users/carmenst/Documents/Polybox/WCR/Conferences_summerschools/2022_10_Como/Project/Sequential_heat_crops/Data/Climate/CMIP6/future/output/**/*' 
# compute_difference(path_future_output)


"""Compute value for each county"""
path_difference = '/Users/carmenst/Documents/Polybox/WCR/Conferences_summerschools/2022_10_Como/Project/Sequential_heat_crops/Data/Climate/CMIP6/future-historic/output/**/*' 
# # # crops_shapefile = compute_mean_per_county(path_difference)
# # # crops_shapefile.to_file('counties_multi_model_mean.shp')
crops_shapefile = compute_diff_per_model(path_difference)
crops_shapefile.to_file('counties_difference_per_model.shp')

#compute exceedances of thr
path_historic = '/Users/carmenst/Documents/Polybox/WCR/Conferences_summerschools/2022_10_Como/Project/Sequential_heat_crops/Data/Climate/CMIP6/historic/input/**/*' 
nr_years = 40
percentiles = [50, 75]
references = ['historic', 'future']
#ssps= ['ssp245', 'ssp370']
ssps= ['ssp119', 'ssp126']

us_shp_file = '/Users/carmenst/Documents/Polybox/WCR/Conferences_summerschools/2022_10_Como/Project/Sequential_heat_crops/Data/shapefiles/input/gadm36_USA_2.shp'
shp_US = gpd.read_file(us_shp_file)


for percentile in percentiles: 
    percentile_str = '/percentile_'+str(percentile)
    for reference in references:
        reference_str = '/reference_'+reference
        
        # frequency
        compute_frequency(path_historic, nr_years, percentile, reference, ssps)
        
        # frequency change
        path_scenario = '/Users/carmenst/Documents/Polybox/WCR/Conferences_summerschools/2022_10_Como/Project/Sequential_heat_crops/Data/Climate/CMIP6/frequencies' +(
            reference_str+ percentile_str +'/exceedances')
        compute_frequency_change(path_scenario)
        
        for ssp in ssps:
            # occurence of hot-hot events in the past and future
            path_cooccurrence = '/Users/carmenst/Documents/Polybox/WCR/Conferences_summerschools/2022_10_Como/Project/Sequential_heat_crops/Data/Climate/CMIP6/frequencies' + (
                reference_str+ percentile_str +'/frequency')
            comparison_fut_his(path_cooccurrence, ssp)
            
            # compare changes in co-occurence
            path = '/Users/carmenst/Documents/Polybox/WCR/Conferences_summerschools/2022_10_Como/Project/Sequential_heat_crops/Data/Climate/CMIP6/frequencies'+(
                reference_str+percentile_str +'/frequency increase/'+ssp)
            compute_f_mean(path)

# import numpy as np
# import geopandas as gpd
# import matplotlib.pyplot as plt
# import cartopy.crs as ccrs
# import copy
# # us_shp_file ='/Users/carmenst/Documents/Polybox/WCR/Conferences_summerschools/2022_10_Como/Project/Sequential_heat_crops/Data/shapefiles/output/2023_06_T_diff/counties_difference_per_model.shp'

# scenario = 'f_change_relative_ssp370_reference_historic_percentile_50.shp'
# us_shp_file ='/Users/carmenst/Documents/Polybox/WCR/Conferences_summerschools/2022_10_Como/Project/Sequential_heat_crops/Data/shapefiles/output/2023_06_f_change/'+scenario
# shp_US = gpd.read_file(us_shp_file)
 

# shp_US.plot(column='%M_incre', cmap ='seismic' ,legend=True, vmin=0, vmax=1)

# mean = shp_US['mean_incre'].values
# perc_increase = shp_US['%M_incre'].values
# idx_decrease  = np.where(mean <= 0)
# percentage_agreement =  copy.deepcopy(perc_increase)
# percentage_agreement[idx_decrease] = 1-perc_increase[idx_decrease]
# shp_US['%M_agree'] = percentage_agreement
# shp_US.plot(column='%M_agree', cmap ='seismic' ,legend=True, vmin=0, vmax=1)
# shp_US.to_file('/Users/carmenst/Desktop/'+scenario)

# missing_kwds = dict(color='grey', label='No Data')
# shp_US.plot(column = '245_su_M0', missing_kwds=missing_kwds)

# shp_US.plot(column = 'mean_incre')
