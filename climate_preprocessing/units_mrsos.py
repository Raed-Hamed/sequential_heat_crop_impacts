import os
import sys
import pathlib
import json
import netCDF4
import xarray as xr

dir_preprocessing = os.path.join(str(pathlib.Path().resolve()))
sys.path.append(dir_preprocessing)


geo_area = 'us'
var_input = 'mrsos'


with open(f'config_{geo_area}.json', 'r') as file:
    config_data = json.load(file)


dir_mrsos = config_data['input']['dir_cmip_var'].format(variable=var_input)
historic_dir = os.path.join(dir_mrsos, 'historical')

files = os.listdir(historic_dir)
files = [file for file in files if os.path.isfile(os.path.join(historic_dir, file)) and file != '.DS_Store']


for file in files:
    _, model, _, _, _ = file.split('_')
    file_path = os.path.join(historic_dir, file)
    nc_file = netCDF4.Dataset(file_path, mode='r')
    mrsos_var = nc_file.variables['mrsos']
    
    
    dataset  = xr.open_dataset(file_path)
    
    try: 
        print(model, ',units soil moisture:',  mrsos_var.units, ', depth: ', dataset.depth.values)
    except:
        print(model, ',units soil moisture:',  mrsos_var.units, ', no depth')


