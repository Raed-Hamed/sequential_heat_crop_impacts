import os
import re
import numpy as np
import xarray as xr
from pathlib import Path


SEASONS = ['spring', 'summer']
DENSITY_WATER = 1000 #kg/m3


def extract_bbox(file, lon_bounds, lat_bounds):
    
    dataset  = xr.open_dataset(file)
    
    # Adjust the longitude values from 0-360 to -180-180
    dataset['lon'] = xr.where(dataset['lon'] > 180, dataset['lon'] - 360, dataset['lon'])
    
    # Sort the longitude coordinates
    ds = dataset.sortby('lon')
    
    study_area = ds.sel(lon=slice(lon_bounds[0], lon_bounds[1]), 
                            lat=slice(lat_bounds[0], lat_bounds[1])) 
    
    return study_area


def monthly_means_from_annual_files(path, lon_bounds, lat_bounds, months, years, save=True, filenames_dict=None):
    files = [f.name for f in path.iterdir() if f.is_file() and not f.name.startswith('.')]
    files.sort()

    datasets = []
    for file in files:
        print("Processing file", file)
        file_path = Path(path, file)
        ds_file = compute_monthly_mean(file_path, lon_bounds, lat_bounds, months)
        datasets.append(ds_file)
        
    results_ds = xr.concat(datasets, dim='time')
    
    t_monthly = results_ds.rename_dims({'time': 'year'})
    t_monthly['year'] = np.arange(years[0], years[1])
    
    if save:
        # Ensure the directory exists
        Path(filenames_dict['calc_his']['dir_calc_his']).mkdir(parents=True, exist_ok=True)
        t_monthly.to_netcdf(filenames_dict['calc_his']['monthly_file'].format(variable='tasmax'))
        t_monthly.close()
    
    return t_monthly

def compute_monthly_mean(file, lon_bounds, lat_bounds, months, variable='tmax'):
                                         
    study_area = extract_bbox(file, lon_bounds, lat_bounds)
    
    # Select data for the months from March to August
    months_selection = study_area.sel(time=study_area['time.month'].isin(months))
    
    # Compute mean per month
    monthly_mean = months_selection.groupby('time.month').mean(dim='time')
    t_monthly = monthly_mean.rename_vars({variable: 'monthly_mean'})  
    
    return t_monthly


def annual_seasonal_means(t_monthly, crops_dict, save=True, filenames_dict=None):
    """Simplified version of annual_t_metrics only computing seasonal mean over two months"""
    sum_metric = []
    name_metric = []
    crops = list(crops_dict.keys())
    for crop in crops:
        crop_var = crops_dict[crop]
        for season in SEASONS:
            #monthly mean across months
            metric_values = t_monthly['monthly_mean']
            sum_metric.append(metric_values.sel(month=slice(*crop_var[f'{season}'])).groupby('year').mean('month'))
            name_metric.append(crop+'_'+season)
            
    # Create an empty dataset
    t_annual = xr.Dataset()
    # Add variables to the dataset
    for name, metric in zip(name_metric, sum_metric):
        t_annual[name] = metric
    
    t_annual.to_netcdf(filenames_dict['calc_his']['annual_file'].format(variable='tasmax'))
    t_annual.close()
    
    return t_annual


def metrics_from_timperiod(file, lon_bounds, lat_bounds, crops_dict, variable='SMroot', output_file=None):

    study_area = extract_bbox(file, lon_bounds, lat_bounds)
    
    sum_metric = []
    name_metric = []
    crops = list(crops_dict.keys())
    for crop in crops:
        crop_var = crops_dict[crop]
        for season in SEASONS:
            #monthly mean across months
            months = crop_var[f'{season}']
            months_selection = study_area.sel(time=study_area['time.month'].isin(months))
            monthly_means = months_selection.groupby('time.year').mean('time')[variable]
            
            if variable == 'mrsos':
                sm_m3 = monthly_means / (DENSITY_WATER*months_selection.depth.values)
                sum_metric.append(sm_m3)
            else:
                sum_metric.append(monthly_means)
            name_metric.append(crop+'_'+season)

    # Create an empty dataset
    sm_annual = xr.Dataset()
    # Add variables to the dataset
    for name, metric in zip(name_metric, sum_metric):
        sm_annual[name] = metric
    
    if output_file is not None:
        sm_annual.to_netcdf(output_file)
        sm_annual.close()
    
    return sm_annual


"""CMIP 6"""

def get_model_list(config_data, var_input='tasmax'):
    """Get dictionary with all CMIP6 files per climate model"""    
    dir_cmip = config_data['input']['dir_cmip_var'].format(variable=var_input)
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
                file_structure = f'{var_input}_([^_]+)_ssp\d+.*\.nc'
                match = re.match(file_structure, filename)
                if match and match.group(1) == climate_model:
                    future_file = os.path.join(root, filename)
                    climate_model_files[climate_model]["future_files"].append(future_file)

    return climate_model_files

def get_monthly_sub(t_his_path, lon_bounds, lat_bounds, crops_dict, variable='tasmax', nr_years=40):
    
    t_monthly = metrics_from_timperiod(t_his_path, lon_bounds, lat_bounds, crops_dict, variable)
    
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
        coocurrence = both.sum(dim='year')
        
        # Store the result in the dictionary
        concurrence_results[f"{crop}_coocurrence"] = coocurrence
    
    # Create a new Dataset to store these results
    concurrence_dataset = xr.Dataset(concurrence_results, coords={"lat": binary_data.lat, "lon": binary_data.lon})
    
    return concurrence_dataset

def extract_ssp(file_path):
    match = re.search(r'ssp(\d{3})', file_path)
    return match.group(0) if match else None

def get_data_model(climate_model_files, quantile, lon_bounds, lat_bounds, 
                   crops_dict, variable, filename_dict, nr_years=40, adapt_months=0):
    
    his_path = climate_model_files['historic_file']
    var_historic = get_monthly_sub(his_path, lon_bounds, lat_bounds, crops_dict, variable, nr_years)
    # compute the quantile
    his_quantile = var_historic.quantile(q=quantile, dim='year') 
    mean_his = var_historic.mean(dim='year')
    
    crops = list(crops_dict.keys())
    
    cooccurence_his = compute_cooccurence(var_historic, his_quantile, crops)
    
    list_data_delta = []
    list_data_frequency = []
    for file in climate_model_files['future_files']:
        ssp = extract_ssp(file)
        crops_dict_fut = adapt_growing_season(crops_dict, adapt_months)
        var_future = get_monthly_sub(file, lon_bounds, lat_bounds, crops_dict_fut, variable, nr_years)
        
        mean_fut = var_future.mean(dim='year')
        diff_mean = mean_fut-mean_his
        # Create a dictionary to map old variable names to new variable names
        rename_dict = {var: ('{var}_{ssp}').format(var=var, ssp=ssp) for var in diff_mean.data_vars}
        # Rename the variables in the Dataset
        diff_mean_renamed = diff_mean.rename(rename_dict)
        list_data_delta.append(diff_mean_renamed)
        
        cooccurence_fut_ref_his = compute_cooccurence(var_future, his_quantile, crops)
        diff_freq_his = (cooccurence_fut_ref_his/cooccurence_his)-1
        # Create a dictionary to map old variable names to new variable names
        rename_dict = {var: ('{var}_{ssp}_his-reference').format(var=var, ssp=ssp) for var in diff_freq_his.data_vars}
        # Rename the variables in the Dataset
        diff_freq_his_renamed = diff_freq_his.rename(rename_dict)
        
        fut_quantile = var_future.quantile(q=quantile, dim='year') 
        cooccurence_fut_ref_fut = compute_cooccurence(var_future, fut_quantile, crops)
        diff_freq_fut = (cooccurence_fut_ref_fut/cooccurence_his)-1
        # Create a dictionary to map old variable names to new variable names
        rename_dict = {var: ('{var}_{ssp}_fut-reference').format(var=var, ssp=ssp) for var in diff_freq_fut.data_vars}
        # Rename the variables in the Dataset
        diff_freq_fut_renamed = diff_freq_fut.rename(rename_dict)
        
        merged = xr.merge([diff_freq_his_renamed, diff_freq_fut_renamed])
        
        list_data_frequency.append(merged)
    
    data_delta = xr.merge(list_data_delta)  
    data_frequency = xr.merge(list_data_frequency)  
    
    return data_delta, data_frequency

def adapt_growing_season(crops_dict, nr_months):
    new_dict = {}
    for crop, seasons in crops_dict.items():
        new_dict[crop] = {}
        for season, values in seasons.items():
            # Subtract 1 from each integer in the list and store in the new dictionary
            new_dict[crop][season] = [value - nr_months for value in values]
    return new_dict

def wrap_delta_frequency(config_data, var_input, nr_months):
    
    lon = config_data['study_area']['lon']
    lat = config_data['study_area']['lat']
    nr_years = config_data['study_area']['nr_years']
    quantiles = config_data['study_area']['quantiles']
    crops_dict = config_data['study_area']['crops_dict']
    
    climate_model_files = get_model_list(config_data, var_input)
    climate_models = list(climate_model_files.keys())

    for quantile in quantiles: 
        Path(config_data['cmip6_grid']['dir_delta'].format(variable=var_input, nr_months=nr_months)).mkdir(parents=True, exist_ok=True)
        Path(config_data['cmip6_grid']['dir_frequency'].format(nr_months=nr_months, percentile=int(quantile*100)), 
             ).mkdir(parents=True, exist_ok=True)
        for model in climate_models:  
            data_delta, data_frequency = get_data_model(climate_model_files[model], quantile, lon, 
                                                            lat, crops_dict, var_input, config_data, nr_years,
                                                            nr_months)

            data_delta.to_netcdf(config_data['cmip6_grid']['delta_per_model'].format(nr_months=nr_months, 
                                                                                     model=model, 
                                                                                     variable=var_input))
            data_frequency.to_netcdf(config_data['cmip6_grid']['frequency_per_model'].format(nr_months=nr_months, 
                                                                                             percentile=int(quantile*100), 
                                                                                             model=model, variable=var_input))
            
def wrap_delta_sm(config_data, var_input='mrsos', nr_months=0):
    
    lon = config_data['study_area']['lon']
    lat = config_data['study_area']['lat']
    nr_years = config_data['study_area']['nr_years']
    crops_dict = config_data['study_area']['crops_dict']
    
    climate_model_files = get_model_list(config_data, var_input)
    climate_models = list(climate_model_files.keys())


    Path(config_data['cmip6_grid']['dir_delta'].format(variable=var_input, nr_months=nr_months)).mkdir(parents=True, exist_ok=True)
    
    quantile = 0.5
    for model in climate_models:  
        data_delta, data_frequency = get_data_model(climate_model_files[model], quantile, lon, 
                                                        lat, crops_dict, var_input, config_data, nr_years,
                                                        nr_months)

        data_delta.to_netcdf(config_data['cmip6_grid']['delta_per_model'].format(nr_months=nr_months, 
                                                                                 model=model, 
                                                                                 variable=var_input))

