import os
import json

CROPS_DICT = {'us':
                  {'wheat': {'spring': [3,4],
                            'summer': [5,6]
                            },
                  'maize': {'spring': [5,6],
                            'summer': [7,8],
                           },
                  'soybean': {'spring': [5,6],
                            'summer': [7,8],
                           }
                  },
              'eu': 
                  {'wheat': {'spring': [3,4],
                            'summer': [5,6]
                            },
                  'maize': {'spring': [5,6],
                            'summer': [7,8],
                           },
                  'sunflower': {'spring': [5,6],
                            'summer': [7,8],
                           }},
              }


def create_config(dir_data, geo_area, months, lon_bounds, lat_bounds, his_years=(1980, 2022), 
                  nr_years_delta=40, quantiles=[0.5, 0.75, 0.9]):
    
    crops_dict = CROPS_DICT[geo_area]
    
    dir_input = os.path.join(f'{dir_data}', 'input')
    tasmax = os.path.join(dir_input, 'CPC')
    gleam = os.path.join(dir_input, 'SMroot_1980-2021_GLEAM_v3.6a_daily_remap05.nc')
    shape_crop = os.path.join(dir_input, 'shapefiles', f'{geo_area}_'+'{crop}_shapefile', f'{geo_area}_'+'{crop}_cropping_regions.shp')
    dir_cmip_var = os.path.join(dir_input,'CMIP6', '{variable}')
    
    dir_calc = os.path.join(f'{dir_data}', 'output', 'calc')
    
    dir_calc_his = os.path.join(dir_calc, 'historic')
    monthly_file = os.path.join(dir_calc_his, '{variable}_monthly_'+f'{geo_area}.nc')
    annual_file = os.path.join(dir_calc_his, '{variable}_annual_'+f'{geo_area}.nc')
    
    dir_calc_cmip6 = os.path.join(dir_calc, 'future', f'{geo_area}')
    dir_calc_cmip6grid = os.path.join(dir_calc_cmip6, 'CMIP6_grid')
    dir_delta = os.path.join(dir_calc_cmip6grid, 'delta_{variable}') 
    delta_per_model = os.path.join(dir_delta, '{model}_{variable}_delta.nc') 
    dir_frequency = os.path.join(dir_calc_cmip6grid, 'frequency_percentile_{percentile}') 
    frequency_per_model = os.path.join(dir_frequency, '{model}_{variable}_frequency_{percentile}.nc')

    dir_calc_cmip6counties = os.path.join(dir_calc_cmip6, 'Per_county')
    dir_calc_cmip6_county_crop = os.path.join(dir_calc_cmip6counties, '{crop}')
    
    dir_delta_county = os.path.join(dir_calc_cmip6_county_crop, 'delta_{variable}') 
    delta_per_model_county = os.path.join(dir_delta_county, '{model}_{variable}_delta.csv') 
    dir_frequency_county = os.path.join(dir_calc_cmip6grid, 'frequency_percentile_{percentile}') 
    frequency_per_model_county = os.path.join(dir_frequency_county, '{model}_{variable}_frequency_{percentile}.csv')
    
    dir_geo = os.path.join(f'{dir_data}', 'output', f'{geo_area}')
    dir_geo_crop = os.path.join(dir_geo, '{crop}')
    csv_his = os.path.join(dir_geo_crop, 'historic_{crop}_'+f'{geo_area}.csv')
    
    csv_fut = os.path.join(dir_geo_crop, 'delta_{crop}_'+f'{geo_area}.csv')
    
    
    config_data = {'input': {'tasmax': tasmax,
                                'gleam': gleam,
                                'shape_crop':shape_crop,
                                'dir_cmip_var': dir_cmip_var,
                                },
                      'study_area': {'geo_area': geo_area,
                                     'crops_dict': crops_dict,
                                     'months': months,
                                     'lon': lon_bounds, 
                                     'lat': lat_bounds,
                                     'his_years': his_years,
                                     'nr_years': nr_years_delta,
                                     'quantiles': quantiles,
                                },
                      
                      'calc_his': {'dir_calc_his': dir_calc_his,
                                'monthly_file': monthly_file, 
                                'annual_file': annual_file, 
                                'dir_geo': dir_geo,
                                'dir_geo_crop': dir_geo_crop,
                                'csv_his': csv_his,
                                },
                      'calc_cmip6_grid':{
                                
                                'dir_calc_cmip6':dir_calc_cmip6,
                                'dir_calc_cmip6grid': dir_calc_cmip6grid,
                                
                                'dir_delta': dir_delta,
                                'dir_frequency':dir_frequency,
                                'delta_per_model': delta_per_model,
                                'frequency_per_model': frequency_per_model,
                                },
                        'calc_cmip6_counties':{
                                'dir_delta_county': dir_delta_county,
                                'dir_frequency_county':dir_frequency_county,
                                'dir_calc_cmip6_county_crop': dir_calc_cmip6_county_crop,
                                'delta_per_model_county': delta_per_model_county,
                                'frequency_per_model_county': frequency_per_model_county,
                                'csv_fut':csv_fut,
                                },
        }
    
    
    with open(f'config_{geo_area}.json', 'w') as file:
        json.dump(config_data, file, indent=4)
    
    
    return config_data