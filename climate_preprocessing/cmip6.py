from pathlib import Path
import os
import sys
import pathlib
import json
import pandas as pd
import xarray as xr

dir_preprocessing = os.path.join(str(pathlib.Path().resolve()))
sys.path.append(dir_preprocessing)

import functions as fct

geo_area = 'us'
var_input = 'tasmax'

with open(f'config_{geo_area}.json', 'r') as file:
    config_data = json.load(file)


lon = config_data['study_area']['lon']
lat = config_data['study_area']['lat']
months = config_data['study_area']['months']
years = config_data['study_area']['his_years']
nr_years = config_data['study_area']['nr_years']
quantiles = config_data['study_area']['quantiles']
crops_dict = config_data['study_area']['crops_dict']
crops = list(config_data['study_area']['crops_dict'].keys())

#Get all climate models
climate_model_files = fct.get_model_list(config_data['input']['dir_cmip_var'].format(variable=var_input))
climate_models = list(climate_model_files.keys())
print(climate_models)

"""Wrapper on CMIP6 grid"""
for quantile in quantiles: 
    Path(config_data['calc_cmip6_grid']['dir_delta'].format(variable=var_input)).mkdir(parents=True, exist_ok=True)
    Path(config_data['calc_cmip6_grid']['dir_frequency'].format(percentile=int(quantile*100))).mkdir(parents=True, exist_ok=True)
    for model in climate_models:  
        data_delta, data_frequency = fct.get_data_model(climate_model_files, model, quantile, lon, 
                                                        lat, crops_dict, var_input, config_data, nr_years)

        # data_delta.to_netcdf(config_data['calc_cmip6_grid']['delta_per_model'].format(model=model, variable=var_input))
        data_frequency.to_netcdf(config_data['calc_cmip6_grid']['frequency_per_model'].format(percentile=int(quantile*100), model=model, variable=var_input))


percentile = 0.9
data_delta  = xr.open_dataset(config_data['calc_cmip6_grid']['frequency_per_model_county'].format(model=model, variable=var_input, percentile=percentile))
 
"""CSV per model and crop"""
for model in climate_models:
    data_delta  = xr.open_dataset(config_data['calc_cmip6_grid']['delta_per_model'].format(model=model, variable=var_input))
   
    ds_var_list = list(data_delta.data_vars.keys())

    for crop in crops:
        list_data_crop = []
        shape_usa_file = config_data['input']['shape_crop'].format(crop=crop)
        crop_vars = [var for var in ds_var_list if crop in var]
        
        for var_name in crop_vars:
                
            crop, season, ssp = var_name.split('_')
            
            gdf_model = fct.fct_mean_per_county_future(data_delta, var_name, shape_usa_file)
            # result = fct.spatial_interpolation_future(result0, shape_usa_file, var_name, data_delta)
            
            final_ds = gdf_model.rename(columns={var_name:'T'})
            final_ds['season'] = season
            final_ds['ssp'] = ssp
            
            list_data_crop.append(final_ds)
        
        data_crop = pd.concat(list_data_crop, ignore_index=True)
        final_data_csv = data_crop.drop(columns='geometry')
        
        
        Path(config_data['calc_cmip6_counties']['dir_delta_county'].format(crop=crop, variable=var_input)).mkdir(parents=True, exist_ok=True)
        final_data_csv.to_csv(config_data['calc_cmip6_counties']['delta_per_model_county'].format(variable=var_input, model=model, crop=crop), index=False)
    

"""CSV per crop"""
for crop in crops:
    list_data_crop = []
    for model in climate_models:
        filename= config_data['calc_cmip6_counties']['delta_per_model_county'].format(variable=var_input, model=model, crop=crop)
        df = pd.read_csv(filename)
        df['model'] = model
        list_data_crop.append(df)

    data_crop = pd.concat(list_data_crop, ignore_index=True)
    data_crop.to_csv(config_data['calc_cmip6_counties']['csv_fut'].format(geo_area='us', crop=crop), index=False)
    







"""In case I need to repeat the calculations for soil moisture"""
# for idx, variable in enumerate(list_variables): 
#     dir_variable = Path(os.path.join(dir_cmip, variable))
#     climate_model_files = fct.get_model_list(dir_variable)
#     climate_models = list(climate_model_files.keys())
    
#     delta_dir = Path(os.path.join(dir_data, 'output', 'future', f'delta_{variable}'))
#     Path(delta_dir).mkdir(parents=True, exist_ok=True)

    




