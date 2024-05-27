import numpy as np
import xarray as xr
from pathlib import Path


"""Seasonal average / sum"""
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





def compute_metrics(file, lon_bounds, lat_bounds, months, thresholds):
    
    ds  = xr.open_dataset(file)
    

    
    if lon_bounds[0] > lon_bounds[1]:
        continental_us = ds.sel(lon=slice(lon_bounds[0], 360), 
                                     lat=slice(lat_bounds[0], lat_bounds[1])
                                     ).combine_first(ds.sel(
                                         lon=slice(0, lon_bounds[1]),
                                         lat=slice(lat_bounds[0], lat_bounds[1])))
    else:
        # Select data within the specified bounds
        continental_us = ds.sel(lon=slice(lon_bounds[0], lon_bounds[1]), 
                                lat=slice(lat_bounds[0], lat_bounds[1])) 
                                         
    
    # Select data for the months from March to August
    months_selection = continental_us.sel(time=continental_us['time.month'].isin(months))
    
    # Combine year and month into a single string
    # time_str = [f"{year}-{month:02d}" for year, month in zip(months_selection['time.year'].values, months_selection['time.month'].values)]
    
    # Assign the combined year-month string as the new coordinate for the time dimension
    # months_selection['time'] = pd.to_datetime(time_str)
    
    # Compute mean per month
    monthly_mean = months_selection.groupby('time.month').mean(dim='time')
    output_ds = monthly_mean.rename_vars({'tmax': 'monthly_mean'})
    
    
    """Sum above threshold"""
    for threshold in thresholds: 
        # Compute difference between values and threshold, and keep positive differences
        values_above_threshold = (months_selection - threshold).clip(min=0).fillna(0)
        
        # Sum values above threshold per month
        values_above_threshold_per_month = values_above_threshold.groupby('time.month').sum(dim='time')
        
        output_ds['values_above_thr_'+str(threshold)] = values_above_threshold_per_month.tmax
        """Days above threshold"""
        # Create a boolean mask indicating values above the threshold
        above_threshold_mask = months_selection > threshold
        # Count the number of days per month that exceed the threshold
        days_above_threshold_per_month = above_threshold_mask.groupby('time.month').sum(dim='time')
        
        output_ds['days_above_thr_'+str(threshold)] = days_above_threshold_per_month.tmax
    
    # # Convert time values to 'year-month' format
    # year = np.unique(continental_us['time.year'].values)[0]
    # year_month = [f"{year}-{month:02d}" for month in output_ds.month.values]
    # output_ds['month'] = year_month
    
    return output_ds

def monthly_metrics(path, lon_bounds, lat_bounds, months, thresholds):
    files = [f.name for f in path.iterdir() if f.is_file() and not f.name.startswith('.')]
    files.sort()

    datasets = []
    for file in files:
        print(file)
        file_path = Path(path, file)
        ds_file = compute_metrics(file_path, lon_bounds, lat_bounds, months, thresholds)
        datasets.append(ds_file)
        
    results_ds = xr.concat(datasets, dim='time')
    
    output = results_ds.rename_dims({'time': 'year'})
    output['year'] = np.arange(1980, 2022)
    
    return output

def extract_metrics(results_ds, crops, seasons, metrics, season_defs):
    sum_metric = []
    name_metric = []
    for crop in crops:
        crop_var = CROPS_DICT[crop]
        for season in seasons:
            for metric in metrics:
                metric_values = results_ds[metric+str(crop_var[f'threshold_{season}'])]
                sum_metric.append(metric_values.sel(month=slice(*crop_var[f'{season}_season'])).groupby('year').sum('month'))
                name_metric.append(crop+'_'+season+'_'+metric+'season')
                
                #per month
                sum_metric.append(metric_values.sel(month=crop_var[f'{season}_month']))
                name_metric.append(crop+'_'+season+'_'+metric+'month')
            
            #monthly mean across months
            metric_values = results_ds['monthly_mean']
            sum_metric.append(metric_values.sel(month=slice(*crop_var[f'{season}_season'])).groupby('year').mean('month'))
            name_metric.append(crop+'_'+season+'_'+'monthly_mean_'+'season')
            
            #monthly mean per month
            metric_values = results_ds['monthly_mean']
            sum_metric.append(metric_values.sel(month=crop_var[f'{season}_month']))
            name_metric.append(crop+'_'+season+'_'+'monthly_mean_'+'month')
            
    # Create an empty dataset
    final_dataset = xr.Dataset()
    # Add variables to the dataset
    for name, metric in zip(name_metric, sum_metric):
        final_dataset[name] = metric

    final = final_dataset.drop_vars('month')
    
    
    return final


def extract_soilmoisture(path_moisture, lon_bounds, lat_bounds):
    
    ds_moisture = xr.open_dataset(path_moisture)
    
    continental_us = ds_moisture.sel(lon=slice(lon_bounds[0], lon_bounds[1]), 
                            lat=slice(lat_bounds[0], lat_bounds[1])) 

    sum_metric = []
    name_metric = []
    for crop in crops:
        crop_var = CROPS_DICT[crop]
        for season in seasons:
            #monthly mean across months
            months = crop_var[f'{season}_season']
            months_selection = continental_us.sel(time=continental_us['time.month'].isin(months))
            sum_metric.append(months_selection.groupby('time.year').mean('time')['SMroot'])
            name_metric.append(crop+'_'+season+'_'+'monthly_mean_'+'season')
            
            #monthly mean per month
            months = crop_var[f'{season}_month']
            months_selection = continental_us.sel(time=continental_us['time.month'].isin(months))
            sum_metric.append(months_selection.groupby('time.year').mean('time')['SMroot'])
            name_metric.append(crop+'_'+season+'_'+'monthly_mean_'+'month')
            
    # Create an empty dataset
    final_dataset = xr.Dataset()
    # Add variables to the dataset
    for name, metric in zip(name_metric, sum_metric):
        final_dataset[name] = metric
    
    return final_dataset



"""Compute metrics for geographical area and all months needed for extracting 
final metrics of interest"""
path = Path('/Users/carmenst/Desktop/Sequential_Heat/tmax_365')
months = [3, 4, 5, 6, 7, 8]
thresholds = [15, 20.5, 28.8, 30.5]

# # Define latitude and longitude bounds for the continental US
# lon_bounds = [235, 294]
# # lon_bounds = [-125, -66]
# lat_bounds = [50, 24]
# geo_area = 'USA'

#europe
lon_bounds = [350, 40]
lat_bounds = [70, 35]
geo_area = 'EU'


"""Compute for a single file"""
file = '/Users/carmenst/Desktop/tmax.1980.nc'
dict_us = compute_metrics(file, lon_bounds+360, lat_bounds, months, thresholds)
dict_us['days_above_thr_15'].isel(month=0).plot()

"""Compute for all files / years"""
file_metrics = Path('/Users/carmenst/Desktop/Sequential_Heat', f'EDD_{geo_area}.nc')
try:
    results_ds = xr.open_dataset(file_metrics)
except:
    results_ds = monthly_metrics(path, lon_bounds, lat_bounds, months, thresholds)
    results_ds.to_netcdf(file_metrics)
    results_ds.close()


"""Extract certain values for the different crops"""
crops = ['maize', 'wheat']

seasons = ['spring', 'summer']
metrics = ['values_above_thr_', 'days_above_thr_']
season_defs = ['season']

file_extract = Path('/Users/carmenst/Desktop/Sequential_Heat', f'EDD_annual_{geo_area}.nc')
try:
    final = xr.open_dataset(file_extract)
except:
    final = extract_metrics(results_ds, crops, seasons, metrics, season_defs)
    final.to_netcdf(file_extract)
    final.close()
    
    
"""Extract soil moisture for the different seasons"""
path_moisture = '/Users/carmenst/Desktop/Sequential_Heat/SMroot_1980-2021_GLEAM_v3.6a_daily_remap05.nc'
lon_bounds = [-10, 40]
sm_ds = extract_soilmoisture(path_moisture, lon_bounds, lat_bounds)
file_moisture = Path('/Users/carmenst/Desktop/Sequential_Heat', f'SM_annual_{geo_area}.nc')
try:
    sm_ds = xr.open_dataset(file_moisture)
except:
    sm_ds.to_netcdf(file_moisture)
    sm_ds.close()



