import os
import sys
import pathlib
import json

dir_preprocessing = os.path.join(str(pathlib.Path().resolve()))
sys.path.append(dir_preprocessing)

import functions as fct
import functions_counties as fct_counties

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


adaptation_months = [0, 1]
for nr_months in adaptation_months:
    nr_months = 0
    fct.wrap_delta_frequency(config_data, var_input, nr_months)
    
    data_source = 'delta'
    fct_counties.wrap_cmip_counties(config_data, var_input, data_source, nr_months=nr_months)
    fct_counties.wrap_cmip_crop(config_data, var_input, data_source, nr_months=nr_months)
    
    data_source = 'percentile'
    for percentile in quantiles:
        fct_counties.wrap_cmip_counties(config_data, var_input, data_source, percentile, nr_months)
        fct_counties.wrap_cmip_crop(config_data, var_input, data_source, percentile, nr_months)

        
    

"""In case I need to repeat the calculations for soil moisture"""
# for idx, variable in enumerate(list_variables): 
#     dir_variable = Path(os.path.join(dir_cmip, variable))
#     climate_model_files = fct.get_model_list(dir_variable)
#     climate_models = list(climate_model_files.keys())
    
#     delta_dir = Path(os.path.join(dir_data, 'output', 'future', f'delta_{variable}'))
#     Path(delta_dir).mkdir(parents=True, exist_ok=True)

    




