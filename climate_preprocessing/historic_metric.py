import os
import re
import numpy as np
import xarray as xr
from pathlib import Path

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
    
    ds  = xr.open_dataset(file)
    
    # Select data within the specified bounds
    if lon_bounds[0] > lon_bounds[1]:
        # As the CPC data coordinates range from 0 to 360, this step might be 
        # necessary in case the study area crosses the meridian line
        study_area = ds.sel(lon=slice(lon_bounds[0], 360), 
                                     lat=slice(lat_bounds[0], lat_bounds[1])
                                     ).combine_first(ds.sel(
                                         lon=slice(0, lon_bounds[1]),
                                         lat=slice(lat_bounds[0], lat_bounds[1])))
    else:
        study_area = ds.sel(lon=slice(lon_bounds[0], lon_bounds[1]), 
                                lat=slice(lat_bounds[0], lat_bounds[1])) 
                                         
    
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


def sm_metrics(path_moisture, lon_bounds, lat_bounds, crops, variable='SMroot'):
    
    ds_moisture = xr.open_dataset(path_moisture)
    
    study_area = ds_moisture.sel(lon=slice(lon_bounds[0], lon_bounds[1]), 
                                 lat=slice(lat_bounds[0], lat_bounds[1])) 

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

# def get_cmip6_per_model(dir_cmip):
#     """Get dictionary with all CMIP6 files per climate model"""    
#     historic_dir = Path(os.path.join(dir_cmip, 'historical'))
#     future_dir = Path(os.path.join(dir_cmip, 'future'))
    
#     climate_model_files = {}
#     # Step 1: Determine climate model used in historic files and store historic file path
#     for root, _, files in os.walk(historic_dir):
#         for filename in files:
#             parts = filename.split('_')
#             if len(parts) >= 2:
#                 climate_model = parts[1]
#                 historic_file = os.path.join(root, filename)
#                 climate_model_files[climate_model] = {"historic_file": historic_file, "future_files": []}
    
#     # Step 2: Find future files for each climate model
#     for climate_model, data in climate_model_files.items():
#         for root, _, files in os.walk(future_dir):
#             for filename in files:
#                 if filename.startswith(f"tasmax_{climate_model}"):
#                     future_file = os.path.join(root, filename)
#                     data["future_files"].append(future_file)

#     return climate_model_files




def get_cmip6_per_model(dir_cmip):
    """Get dictionary with all CMIP6 files per climate model"""    
    historic_dir = Path(os.path.join(dir_cmip, 'historical'))
    future_dir = Path(os.path.join(dir_cmip, 'future'))
    
    climate_model_files = {}
    
    # Step 1: Determine climate model used in historic files and store historic file path
    for root, _, files in os.walk(historic_dir):
        for filename in files:
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

def get_monthly_sub(t_his_path, lon_bounds_t, lat_bounds, crops, variable='tasmax', nr_years=40):
    t_monthly = sm_metrics(t_his_path, lon_bounds_t, lat_bounds, crops, variable)
    t_monthly = t_monthly.drop_vars('height')
    #mean over the chosen time period
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

def get_data_model(climate_model_files, model, quantile, lon_bounds_t, lat_bounds, crops, variable):
    his_path = climate_model_files[model]['historic_file']
    var_historic = get_monthly_sub(his_path, lon_bounds_t, lat_bounds, crops, variable=variable)
    # compute the quantile
    his_quantile = var_historic.quantile(q=quantile, dim='year') 
    mean_his = var_historic.mean(dim='year')
    cooccurence_his = compute_cooccurence(var_historic, his_quantile, crops)
    
    list_data = []
    for file in climate_model_files[model]['future_files']:
        ssp = extract_ssp(file)
        var_future = get_monthly_sub(file, lon_bounds_t, lat_bounds, crops, variable=variable)
        
        mean_fut = var_future.mean(dim='year')
        diff_mean = mean_fut-mean_his
        # Create a dictionary to map old variable names to new variable names
        rename_dict = {var: f"delta_{variable}_{ssp}_{var}" for var in diff_mean.data_vars}
        # Rename the variables in the Dataset
        diff_mean_renamed = diff_mean.rename(rename_dict)
        
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
        
        
        merged = xr.merge([diff_mean_renamed, diff_freq_his_renamed, diff_freq_fut_renamed])
        
        list_data.append(merged)
    
    data_model = xr.merge(list_data)  
    
    return data_model
