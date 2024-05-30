import os
import re
import numpy as np
import xarray as xr
from pathlib import Path
import geopandas as gpd
import xesmf as xe
import pandas as pd

SEASONS = ['spring', 'summer']

CROPS_DICT = {'wheat': {'spring_season': [3,4],
                        'spring_month': 4,
                        'threshold_spring': 15,
                        'summer_season': [5,6],
                        'summer_month': 6,
                        'threshold_summer': 20.5},
              'maize': {'spring_season': [5,6],
                        'spring_month': 6,
                        'threshold_spring': 28.8,
                        'summer_season': [7,8],
                        'summer_month': 8,
                        'threshold_summer': 30.5}
              }


def extract_bbox(file, lon_bounds, lat_bounds):
    
    dataset  = xr.open_dataset(file)
    
    # Adjust the longitude values from 0-360 to -180-180
    dataset['lon'] = xr.where(dataset['lon'] > 180, dataset['lon'] - 360, dataset['lon'])
    
    # Sort the longitude coordinates
    ds = dataset.sortby('lon')
    
    study_area = ds.sel(lon=slice(lon_bounds[0], lon_bounds[1]), 
                            lat=slice(lat_bounds[0], lat_bounds[1])) 
    
    return study_area



def monthly_t_metrics(path, lon_bounds, lat_bounds, months, thresholds, years):
    files = [f.name for f in path.iterdir() if f.is_file() and not f.name.startswith('.')]
    files.sort()

    datasets = []
    for file in files:
        print("Processing file", file)
        file_path = Path(path, file)
        ds_file = compute_t_metrics(file_path, lon_bounds, lat_bounds, months, thresholds)
        datasets.append(ds_file)
        
    results_ds = xr.concat(datasets, dim='time')
    
    output = results_ds.rename_dims({'time': 'year'})
    output['year'] = np.arange(years[0], years[1])
    
    return output

def compute_t_metrics(file, lon_bounds, lat_bounds, months, thresholds, variable='tmax'):
                                         
    study_area = extract_bbox(file, lon_bounds, lat_bounds)
    
    # Select data for the months from March to August
    months_selection = study_area.sel(time=study_area['time.month'].isin(months))
    
    # Compute mean per month
    monthly_mean = months_selection.groupby('time.month').mean(dim='time')
    output_ds = monthly_mean.rename_vars({variable: 'monthly_mean'})

    for threshold in thresholds: 
        """Sum above threshold"""
        # Compute difference between values and threshold, and keep positive differences
        values_above_threshold = (months_selection - threshold).clip(min=0).fillna(0)
        # Sum values above threshold per month
        values_above_threshold_per_month = values_above_threshold.groupby('time.month').sum(dim='time')
        output_ds['values_above_thr_'+str(threshold)] = values_above_threshold_per_month[variable]
        
        """Days above threshold"""
        # Create a boolean mask indicating values above the threshold
        above_threshold_mask = months_selection > threshold
        # Count the number of days per month that exceed the threshold
        days_above_threshold_per_month = above_threshold_mask.groupby('time.month').sum(dim='time')
        output_ds['days_above_thr_'+str(threshold)] = days_above_threshold_per_month[variable]
    
    
    return output_ds

def annual_t_metrics(t_monthly, crops, metrics):
    sum_metric = []
    name_metric = []
    for crop in crops:
        crop_var = CROPS_DICT[crop]
        for season in SEASONS:
            for metric in metrics:
                metric_values = t_monthly[metric+str(crop_var[f'threshold_{season}'])]
                sum_metric.append(metric_values.sel(month=slice(*crop_var[f'{season}_season'])).groupby('year').sum('month'))
                name_metric.append(crop+'_'+season+'_'+metric+'season')
                
                #per month
                sum_metric.append(metric_values.sel(month=crop_var[f'{season}_month']))
                name_metric.append(crop+'_'+season+'_'+metric+'month')
            
            #monthly mean across months
            metric_values = t_monthly['monthly_mean']
            sum_metric.append(metric_values.sel(month=slice(*crop_var[f'{season}_season'])).groupby('year').mean('month'))
            name_metric.append(crop+'_'+season+'_'+'monthly_mean_'+'season')
            
            #monthly mean per month
            metric_values = t_monthly['monthly_mean']
            sum_metric.append(metric_values.sel(month=crop_var[f'{season}_month']))
            name_metric.append(crop+'_'+season+'_'+'monthly_mean_'+'month')
            
    # Create an empty dataset
    final_dataset = xr.Dataset()
    # Add variables to the dataset
    for name, metric in zip(name_metric, sum_metric):
        final_dataset[name] = metric

    final = final_dataset.drop_vars('month')
    
    
    return final


def metrics_from_timperiod(file, lon_bounds, lat_bounds, crops, variable='SMroot'):

    study_area = extract_bbox(file, lon_bounds, lat_bounds)
    
    sum_metric = []
    name_metric = []
    for crop in crops:
        crop_var = CROPS_DICT[crop]
        for season in SEASONS:
            #monthly mean across months
            months = crop_var[f'{season}_season']
            months_selection = study_area.sel(time=study_area['time.month'].isin(months))
            sum_metric.append(months_selection.groupby('time.year').mean('time')[variable])
            name_metric.append(crop+'_'+season)
            #REST OF THE NAME: +'_'+'monthly_mean_'+'season'
            
            # #monthly mean per month
            # months = crop_var[f'{season}_month']
            # months_selection = study_area.sel(time=study_area['time.month'].isin(months))
            # sum_metric.append(months_selection.groupby('time.year').mean('time')[variable])
            # name_metric.append(crop+'_'+season+'_'+'monthly_mean_'+'month')
            
    # Create an empty dataset
    final_dataset = xr.Dataset()
    # Add variables to the dataset
    for name, metric in zip(name_metric, sum_metric):
        final_dataset[name] = metric
    
    return final_dataset


"""CMIP 6"""

def get_model_list(dir_cmip):
    """Get dictionary with all CMIP6 files per climate model"""    
    historic_dir = Path(os.path.join(dir_cmip, 'historical'))
    future_dir = Path(os.path.join(dir_cmip, 'future'))
    
    climate_model_files = {}
    
    # Step 1: Determine climate model used in historic files and store historic file path
    for root, _, files in os.walk(historic_dir):
        for filename in files:
            if filename == '.DS_Store':
                continue  # Skip the ds.store file
            parts = filename.split('_')
            if len(parts) >= 2:
                climate_model = parts[1]
                historic_file = os.path.join(root, filename)
                climate_model_files[climate_model] = {"historic_file": historic_file, "future_files": []}
    
    # Step 2: Find future files for each climate model
    for climate_model in climate_model_files.keys():
        for root, _, files in os.walk(future_dir):
            for filename in files:
                # Extract the climate model part from the filename
                match = re.match(r'tasmax_([^_]+)_ssp\d+.*\.nc', filename)
                if match and match.group(1) == climate_model:
                    future_file = os.path.join(root, filename)
                    climate_model_files[climate_model]["future_files"].append(future_file)

    return climate_model_files

def get_monthly_sub(t_his_path, lon_bounds, lat_bounds, crops, variable='tasmax', nr_years=40):
    
    t_monthly = metrics_from_timperiod(t_his_path, lon_bounds, lat_bounds, crops, variable)
    
    # drop additional coordinate axis if present for this model
    try:
        t_monthly = t_monthly.drop_vars('height')
    except ValueError:
        pass     
    try:
        t_monthly = t_monthly.drop_vars('depth')
    except ValueError:
        pass     
    
    # extract variable for nr-years (last x years in the dataset)
    final_year = t_monthly.year.max()
    start_year = final_year-nr_years+1
    t_sub = t_monthly.sel(year=slice(start_year, final_year))
    
    return t_sub


def compute_cooccurence(t_sub, t_quantile, crops):
    binary_data = xr.where(t_sub >= t_quantile, 1, 0)

    # Initialize an empty dictionary to store the results
    concurrence_results = {}
    
    # Loop over each crop to compute concurrences
    for crop in crops:
        spring_key = f"{crop}_spring"
        summer_key = f"{crop}_summer"
        
        # Logical AND for spring and summer
        both = (binary_data[spring_key] == 1) & (binary_data[summer_key] == 1)
        
        # Count the number of concurrences per lat/lon by summing along the 'year' dimension
        concurrence = both.sum(dim='year')
        
        # Store the result in the dictionary
        concurrence_results[f"{crop}_concurrence"] = concurrence
    
    # Create a new Dataset to store these results
    concurrence_dataset = xr.Dataset(concurrence_results, coords={"lat": binary_data.lat, "lon": binary_data.lon})
    
    
    return concurrence_dataset

def extract_ssp(file_path):
    match = re.search(r'ssp(\d{3})', file_path)
    return match.group(0) if match else None

def get_data_model(climate_model_files, model, quantile, lon_bounds_t, lat_bounds, crops, variable, coord_drop='height', nr_years=40):
    his_path = climate_model_files[model]['historic_file']
    var_historic = get_monthly_sub(his_path, lon_bounds_t, lat_bounds, crops, variable, coord_drop, nr_years)
    # compute the quantile
    his_quantile = var_historic.quantile(q=quantile, dim='year') 
    mean_his = var_historic.mean(dim='year')
    cooccurence_his = compute_cooccurence(var_historic, his_quantile, crops)
    
    list_data_delta = []
    list_data_frequency = []
    for file in climate_model_files[model]['future_files']:
        ssp = extract_ssp(file)
        var_future = get_monthly_sub(file, lon_bounds_t, lat_bounds, crops, variable, coord_drop, nr_years)
        
        mean_fut = var_future.mean(dim='year')
        diff_mean = mean_fut-mean_his
        # Create a dictionary to map old variable names to new variable names
        rename_dict = {var: f"delta_{variable}_{ssp}_{var}" for var in diff_mean.data_vars}
        # Rename the variables in the Dataset
        diff_mean_renamed = diff_mean.rename(rename_dict)
        list_data_delta.append(diff_mean_renamed)
        
        cooccurence_fut_ref_his = compute_cooccurence(var_future, his_quantile, crops)
        diff_freq_his = (cooccurence_fut_ref_his/cooccurence_his)-1
        # Create a dictionary to map old variable names to new variable names
        rename_dict = {var: f"his-reference_{variable}_{ssp}_{var}" for var in diff_freq_his.data_vars}
        # Rename the variables in the Dataset
        diff_freq_his_renamed = diff_freq_his.rename(rename_dict)
        
        fut_quantile = var_future.quantile(q=quantile, dim='year') 
        cooccurence_fut_ref_fut = compute_cooccurence(var_future, fut_quantile, crops)
        diff_freq_fut = (cooccurence_fut_ref_fut/cooccurence_his)-1
        # Create a dictionary to map old variable names to new variable names
        rename_dict = {var: f"fut-reference_{variable}_{ssp}_{var}" for var in diff_freq_fut.data_vars}
        # Rename the variables in the Dataset
        diff_freq_fut_renamed = diff_freq_fut.rename(rename_dict)
        
        merged = xr.merge([diff_freq_his_renamed, diff_freq_fut_renamed])
        
        list_data_frequency.append(merged)
    
    data_delta = xr.merge(list_data_delta)  
    data_frequency = xr.merge(list_data_frequency)  
    
    return data_delta, data_frequency



# """old functions to get mean per county"""

# def fct_mean_per_county(ds, model, var_name, crops_shapefile):
#     ds_out = regrid_cmip6(ds, model, var_name)
#     gdf = create_gdf_from_ds(ds_out, var_name)            
#     points_polys = gpd.sjoin(gdf, crops_shapefile, how="left")
#     stats_pt = points_polys.groupby('index_right')[0].agg(['mean']) #nanmean
#     result = pd.merge(crops_shapefile, stats_pt , left_index=True, right_index=True, how='outer')
#     mean_per_county = np.asarray(result['mean'].values)
    
#     return mean_per_county

# def regrid_cmip6(ds, model, var_name):
    
#     regridder_dir = '/Users/carmenst/Documents/Polybox/WCR/Conferences_summerschools/2022_10_Como/Project/Sequential_heat_crops/Data/Climate/regridder'
#     regridder_file = regridder_dir+'/regridder_'+model+'.nc'
    
#     ds_lat = ds.lat.values
#     ds_lon = ds.lon.values-360
    
#     ds_in = xr.Dataset(
#         {
#             "lat": (["lat"], ds_lat, {"units": "degrees_north"}),
#             "lon": (["lon"], ds_lon, {"units": "degrees_east"}),
#         }
#     )
    
#     ds_out = xr.Dataset(
#         {
#             "lat": (["lat"], np.arange(np.min(ds_lat), np.max(ds_lat), 0.1), {"units": "degrees_north"}),
#             "lon": (["lon"], np.arange(np.min(ds_lon), np.max(ds_lon), 0.1), {"units": "degrees_east"}),
#         }
#     )
    
    
#     if os.path.isfile(regridder_file):
#         regridder = xe.Regridder(ds_in, ds_out, 'conservative', weights=regridder_file)
#     else:
#         regridder = xe.Regridder(ds_in, ds_out, "conservative", ignore_degenerate=True)
#         regridder.to_netcdf(regridder_file)        
    
#     regridded_ds = regridder(ds[var_name].values)    
#     ds_out[var_name] = (('lat', 'lon'), regridded_ds)
    
#     return ds_out

# def create_gdf_from_ds(ds, var_name):
#     lon1, lat1 = np.meshgrid(ds.lon.values, ds.lat.values)
#     geometry = gpd.points_from_xy(lon1.flatten(), lat1.flatten())
#     gdf = gpd.GeoDataFrame(ds[var_name].values.flatten(), crs="EPSG:4326", geometry=geometry)
    
#     return gdf