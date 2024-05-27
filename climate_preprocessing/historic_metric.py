import numpy as np
import xarray as xr
from pathlib import Path


def monthly_t_metrics(path, lon_bounds, lat_bounds, months, thresholds):
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
    output['year'] = np.arange(1980, 2022)
    
    return output

def compute_t_metrics(file, lon_bounds, lat_bounds, months, thresholds):
    
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
    output_ds = monthly_mean.rename_vars({'tmax': 'monthly_mean'})
    
    

    for threshold in thresholds: 
        """Sum above threshold"""
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
    
    
    return output_ds

def annual_t_metrics(t_monthly, crops, seasons, metrics, crops_dict):
    sum_metric = []
    name_metric = []
    for crop in crops:
        crop_var = crops_dict[crop]
        for season in seasons:
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


def sm_metrics(path_moisture, lon_bounds, lat_bounds, crops, seasons, crops_dict):
    
    ds_moisture = xr.open_dataset(path_moisture)
    
    study_area = ds_moisture.sel(lon=slice(lon_bounds[0], lon_bounds[1]), 
                            lat=slice(lat_bounds[0], lat_bounds[1])) 

    sum_metric = []
    name_metric = []
    for crop in crops:
        crop_var = crops_dict[crop]
        for season in seasons:
            #monthly mean across months
            months = crop_var[f'{season}_season']
            months_selection = study_area.sel(time=study_area['time.month'].isin(months))
            sum_metric.append(months_selection.groupby('time.year').mean('time')['SMroot'])
            name_metric.append(crop+'_'+season+'_'+'monthly_mean_'+'season')
            
            #monthly mean per month
            months = crop_var[f'{season}_month']
            months_selection = study_area.sel(time=study_area['time.month'].isin(months))
            sum_metric.append(months_selection.groupby('time.year').mean('time')['SMroot'])
            name_metric.append(crop+'_'+season+'_'+'monthly_mean_'+'month')
            
    # Create an empty dataset
    final_dataset = xr.Dataset()
    # Add variables to the dataset
    for name, metric in zip(name_metric, sum_metric):
        final_dataset[name] = metric
    
    return final_dataset
