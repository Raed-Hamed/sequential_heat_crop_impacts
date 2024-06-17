from pathlib import Path
import os
import sys
import pathlib
import numpy as np

dir_preprocessing = os.path.join(str(pathlib.Path().resolve()))
sys.path.append(dir_preprocessing)
import historic_metric as fct

dir_data = os.path.join(dir_preprocessing, 'data')



lon_bounds = [-125, -66]
lat_bounds = [24, 50]
geo_area = 'us'
years = (1980, 2022)
nr_years = 40 
crops = ['maize', 'wheat', 'soybean']
#iterate over each climate model 
variable = 'tasmax'
# quantiles = [0.5, 0.75, 0.9]
quantiles = [0.5]

from historic_metric import CROPS_DICT
crops_dict = CROPS_DICT[geo_area]

filenames_dict = fct.config_FILEPATHS(dir_data, geo_area)


climate_model_files = fct.get_model_list(filenames_dict['input']['dir_cmip_var'].format(variable=variable))
climate_models = list(climate_model_files.keys())


for quantile in quantiles: 
    Path(filenames_dict['calc']['dir_delta'].format(variable=variable)).mkdir(parents=True, exist_ok=True)
    Path(filenames_dict['calc']['dir_frequency'].format(percentile=int(quantile*100))).mkdir(parents=True, exist_ok=True)
    for model in climate_models:  
        data_delta, data_frequency = fct.get_data_model(climate_model_files, model, quantile, lon_bounds, 
                                                        lat_bounds, crops_dict, variable, filenames_dict, nr_years)

        data_delta.to_netcdf(filenames_dict['calc']['delta_per_model'].format(model=model, variable=variable))
        data_frequency.to_netcdf(filenames_dict['calc']['frequency_per_model'].format(percentile=int(quantile*100), model=model, variable=variable))







import pandas as pd
import xarray as xr

var_input = 'tasmax'


for model in climate_models[-5:]:
    data_delta  = xr.open_dataset(filenames_dict['calc']['delta_per_model'].format(model=model, variable=var_input))
    
    ds_var_list = list(data_delta.data_vars.keys())

    for crop in crops:
        list_data_crop = []
        
        shape_usa_file = filenames_dict['input']['shape_crop'].format(crop=crop)
        
        crop_vars = [var for var in ds_var_list if crop in var]
        
        for var_name in crop_vars:
                
            crop, season, ssp = var_name.split('_')
            
        
            result0 = fct.fct_mean_per_county_future(data_delta, var_name, shape_usa_file)
            result = fct.spatial_interpolation_future(result0, shape_usa_file, var_name, data_delta)
            
            final_ds = result.rename(columns={var_name:'T'})
            final_ds['season'] = season
            final_ds['ssp'] = ssp
            
            list_data_crop.append(final_ds)
        
    
        data_crop = pd.concat(list_data_crop, ignore_index=True)
        
        final_data_csv = data_crop.drop(columns='geometry')
        
        Path(filenames_dict['calc']['dir_delta_county'].format(crop=crop, variable=var_input)).mkdir(parents=True, exist_ok=True)
        final_data_csv.to_csv(filenames_dict['calc']['delta_per_model_county'].format(variable=var_input, model=model, crop=crop), index=False)
    


for crop in crops:
    list_data_crop = []
    for model in climate_models:
        filename= filenames_dict['calc']['delta_per_model_county'].format(variable=var_input, model=model, crop=crop)
        df = pd.read_csv(filename)
        df['model'] = model
        list_data_crop.append(df)

    data_crop = pd.concat(list_data_crop, ignore_index=True)
    data_crop.to_csv(filenames_dict['cmip6']['csv_fut'].format(geo_area='us', crop=crop), index=False)
    


# import matplotlib.pyplot as plt
# # Plot your shapefile
# fig, ax = plt.subplots(figsize=(10, 8))
# result.plot('delta_tasmax_ssp245_maize_spring', ax=ax)

# # Set the xlim and ylim to the bbox
# plt.xlim(lon_bounds[0], lon_bounds[1])
# plt.ylim(lat_bounds[0], lat_bounds[1])

# # Inspect the result
# print(result)




# def compute_counties_future(filenames_dict, ds_var, var_abbrev, crop, 
#                           seasons=['spring', 'summer']):
#     melted_list = []
#     shape_usa_file = filenames_dict['input']['shape_crop'].format(crop=crop)
    
#     crop_vars = [var for var in ds_var_list if crop in var]
#     for variable in enumerate(crop_vars):
#         print('variable: ', variable)

#         merged_data = create_gdf(ds_var, variable, shape_usa_file, years)
        
#         # Melt the DataFrame to have a single column for the years
#         id_vars = [col for col in merged_data.columns if col not in years]
#         melted_data = merged_data.melt(id_vars=id_vars, var_name='year', value_name=var_abbrev[idx_var])

#         # Convert year column to numeric
#         melted_data['year'] = pd.to_numeric(melted_data['year'])
        
#         #melted_data
        
#         melted_list.append(melted_data)
    
#     return melted_list

# import geopandas as gpd

# def create_gdf(ds_var, var_name, shape_usa_file, years):
#     # var_name = crop + '_' + season
#     # crop, season = var_name.split('_')
#     merged_data = gpd.read_file(shape_usa_file)
#     mean_per_county_t = fct.fct_mean_per_county(ds_var, var_name, shape_usa_file)
#     #years = np.arange(1980, 2022) #needs to be adapted
#     for year in years: 
#         year_data = mean_per_county_t[year]
#         merged_data = merged_data.merge(year_data, left_index=True, right_index=True)

#     # merged_data['season'] = season
    
#     return merged_data

























# from pathlib import Path
# import os
# import sys
# import pathlib
# dir_preprocessing = os.path.join(str(pathlib.Path().resolve()))
# dir_data = os.path.join(dir_preprocessing, 'data')
# dir_cmip = Path(os.path.join(dir_data, 'input', 'CMIP6'))

# sys.path.append(dir_preprocessing)
# import historic_metric as fct

# thresholds = [ ]
# lon_bounds = [-125, -66]
# lat_bounds = [24, 50]
# geo_area = 'USA'
# years = (1980, 2022)
# nr_years = 40 
# crops = ['maize', 'wheat']
# #iterate over each climate model 
# list_variables = ['tasmax','mrsos'] 
# coord_drop = ['height', 'depth'] 
# quantiles = [0.5, 0.75, 0.9]


# for idx, variable in enumerate(list_variables): 
#     dir_variable = Path(os.path.join(dir_cmip, variable))
#     climate_model_files = fct.get_model_list(dir_variable)
#     climate_models = list(climate_model_files.keys())
    
#     delta_dir = Path(os.path.join(dir_data, 'output', 'future', f'delta_{variable}'))
#     Path(delta_dir).mkdir(parents=True, exist_ok=True)
    
#     for quantile in quantiles: 
#         quantile_dir = Path(os.path.join(dir_data, 'output', 'future', f'quantile_{quantile}'))
#         Path(quantile_dir).mkdir(parents=True, exist_ok=True)
#         for model in climate_models:  
#             data_delta, data_frequency = fct.get_data_model(climate_model_files, model, quantile, lon_bounds, lat_bounds, crops, variable, coord_drop[idx],nr_years)
            
#             output_delta = Path(os.path.join(delta_dir, f'{model}_{variable}_delta.nc'))
#             data_delta.to_netcdf(output_delta)
            
#             if variable == 'tasmax':
#                 output_frequency = Path(os.path.join(quantile_dir, f'{model}_{variable}_frequency_percentile{quantile}.nc'))
#                 data_frequency.to_netcdf(output_frequency)
  
    



# shape_usa = Path(os.path.join(dir_data, 'input', 'shapefiles', 'gadm36_USA_2.shp'))

# import xarray as xr

# # Example usage
# # Define the path to the NetCDF file and shapefile
# file_path = '/Users/carmenst/Documents/CLIMADA/own_projects/sequential_heat_crop_impacts/climate_preprocessing/data/output/future/delta_tasmax/ACCESS-CM2_tasmax_delta.nc'

# # Load the dataset and shapefile
# ds = xr.open_dataset(file_path)




