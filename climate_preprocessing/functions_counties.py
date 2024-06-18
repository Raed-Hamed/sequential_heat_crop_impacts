import numpy as np
from pathlib import Path
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import xarray as xr

import functions as fct


"""MEAN PER COUNTRY FOR HISTORICAL DATA"""

def compute_counties_hist(filenames_dict, ds_var_list, var_abbrev, crop, 
                          seasons=['spring', 'summer'], years=np.arange(1980, 2022)):
    "Historic wrapper for extracting the variables on county level"
    
    melted_list = []
    shape_usa_file = filenames_dict['input']['shape_crop'].format(crop=crop)
    for idx_var, ds_var in enumerate(ds_var_list):
        print('variable: ', var_abbrev[idx_var])
        for season in seasons:
            print('season: ', season)
            var_name = crop + '_' + season
            merged_data = merge_historic_gdf(ds_var, var_name, shape_usa_file, years)
            merged_data['season'] = season
            
            # Melt the DataFrame to have a single column for the years
            id_vars = [col for col in merged_data.columns if col not in years]
            melted_data = merged_data.melt(id_vars=id_vars, var_name='year', value_name=var_abbrev[idx_var])
    
            # Convert year column to numeric
            melted_data['year'] = pd.to_numeric(melted_data['year'])
            
            melted_list.append(melted_data)
    
    return melted_list


def create_historic_csv(melted_list, crop=None, filenames_dict=None, save=True):

    final_data_t = pd.concat((melted_list[0], melted_list[1]), ignore_index=True)
    final_data_m = pd.concat((melted_list[2], melted_list[3]), ignore_index=True)
    
    common_columns = [col for col in final_data_t.columns if col not in ['T']]
    final_data = pd.merge(final_data_t, final_data_m, on=common_columns, how='inner', suffixes=('_T', '_M'))
    
    # Display the final merged DataFrame
    print(final_data.head())
    
    if save:
        final_data_csv = final_data.drop(columns='geometry')
        Path(filenames_dict['calc_his']['dir_geo_crop'].format(crop=crop)).mkdir(parents=True, exist_ok=True)
        final_data_csv.to_csv(filenames_dict['calc_his']['csv_his'].format(crop=crop), index=False)
    
    return final_data

 
def merge_historic_gdf(ds_var, var_name, shape_usa_file, years):
    merged_data = gpd.read_file(shape_usa_file)
    mean_per_county_t = fct_mean_per_county(ds_var, var_name, shape_usa_file)

    for year in years: 
        year_data = mean_per_county_t[year]
        merged_data = merged_data.merge(year_data, left_index=True, right_index=True)
    
    return merged_data

def fct_mean_per_county(ds, var_name, shape_region):
    crops_shapefile = gpd.read_file(shape_region)
    
    # Create a GeoDataFrame from the dataset
    gdf = create_gdf_from_ds(ds, var_name)
    
    # Spatial join between the points and the county polygons
    points_polys = gpd.sjoin(gdf, crops_shapefile, how="left", predicate='intersects')
    
    # Prepare to find the nearest point
    lon_t, lat_t = np.meshgrid(ds.lon.values, ds.lat.values)
    coord_t = np.vstack((lon_t.reshape(-1), lat_t.reshape(-1))).T
    tree = cKDTree(coord_t)
    
    data = ds[var_name].values.reshape(ds[var_name].shape[0], -1)
    
    # Initialize a list to store the results
    results = []

    # Get the year indices
    year_indices = {year: idx for idx, year in enumerate(ds.year.values)}

    # Loop through each county and compute the mean for each timestep
    for county_index in range(len(crops_shapefile)):
        county_result = []
        for year in ds.year.values:
            year_data = points_polys[points_polys['year'] == year]
            if county_index in year_data.index_right.unique():
                cell_indices = year_data[year_data['index_right'] == county_index].index
                county_mean = year_data.loc[cell_indices, var_name].mean()
            else:
                # Find the centroid of the polygon
                polygon = crops_shapefile.geometry[county_index]
                centroid = polygon.centroid
                centroid_coords = np.array([[centroid.x, centroid.y]])

                # Find the nearest point in the dataset
                nearest_point_idx = find_nearest_point(centroid_coords, tree)
                county_mean = data[year_indices[year], nearest_point_idx]

            county_result.append(county_mean)
        results.append(county_result)
    
    # Convert the results to a DataFrame
    mean_per_county = pd.DataFrame(results, columns=ds.year.values, index=crops_shapefile.index)

    return mean_per_county
    
def create_gdf_from_ds(ds, var_name):
    # Create a meshgrid of the latitude and longitude
    lon1, lat1 = np.meshgrid(ds.lon.values, ds.lat.values)

    # Flatten the meshgrid and create points
    points = [Point(lon, lat) for lon, lat in zip(lon1.flatten(), lat1.flatten())]
    
    if 'year' in ds.dims:
        # Flatten the values and repeat points for each time step
        flattened_values = ds[var_name].values.reshape(ds[var_name].shape[0], -1)
        repeated_points = points * ds[var_name].shape[0]
    
        # Create a GeoDataFrame
        gdf = gpd.GeoDataFrame({var_name: flattened_values.flatten()}, crs="EPSG:4326", 
                               geometry=repeated_points)
    
        # Add a time dimension to the GeoDataFrame
        gdf['year'] = np.repeat(ds['year'].values, ds[var_name].shape[1] * ds[var_name].shape[2])
    
    else:
        gdf = gpd.GeoDataFrame({var_name: ds[var_name].values.flatten()}, crs="EPSG:4326", geometry=points)

    return gdf

def find_nearest_point(centroid_coords, tree):
    """
    Find the index of the nearest point in the dataset to the given coordinates using cKDTree.
    """
    dists, idxs = tree.query(centroid_coords, k=1)
    return idxs[0]


def fct_mean_per_county_future(ds, var_name, shape_region):
    
    crops_shapefile = gpd.read_file(shape_region)
    
    # Create a GeoDataFrame from the dataset
    gdf = create_gdf_from_ds(ds, var_name)
    
    # Spatial join between the points and the county polygons
    points_polys = gpd.sjoin(gdf, crops_shapefile, how="left", predicate='intersects')
    
    # Compute weighted mean for each county
    mean_per_county = points_polys.groupby('index_right')[var_name].mean()
    
    # Merge the results back to the original shapefile
    result0 = pd.merge(crops_shapefile, mean_per_county, left_index=True, right_index=True, how='outer')
    
    result = spatial_interpolation_future(result0, shape_region, var_name, ds)
    
    
    return result

def spatial_interpolation_future(gdf, shape_region, var_name, ds):
    """
    Interpolate values for points within each county polygon.
    """
    
    crops_shapefile = gpd.read_file(shape_region)
    
    lon_t, lat_t = np.meshgrid(ds.lon.values, ds.lat.values)
    coord_t = np.vstack((lon_t.reshape(-1), lat_t.reshape(-1))).T  
    tree = cKDTree(coord_t)
    
    data = ds[var_name].values.reshape(-1)
    
    # Extract the values
    values = gdf[var_name].values
    
    # Prepare a list to store the results
    interpolated_values = []
    
    # Get the indices of polygons with NaN values after mean computation
    nan_indices = gdf.index[gdf[var_name].isnull()].tolist()

    # Interpolate values for polygons with NaN values
    for index, row in crops_shapefile.iterrows():
        polygon = row.geometry
        
        if index in nan_indices:
            # Find the centroid of the polygon
            centroid = polygon.centroid
            centroid_coords = np.array([[centroid.x, centroid.y]])
            
            # Find the nearest point in the temperature dataset
            nearest_point_idx = find_nearest_point(centroid_coords, tree)
            interpolated_value = data[nearest_point_idx]
        else:
            # Use the mean value for other polygons
            interpolated_value = values[index]
        
        interpolated_values.append(interpolated_value)
    
    crops_shapefile[var_name] = interpolated_values
    
    return crops_shapefile



def wrap_cmip_counties(config_data, var_input, data_source, percentile=None, nr_months=0):
    climate_model_files = fct.get_model_list(config_data['input']['dir_cmip_var'].format(variable=var_input))
    climate_models = list(climate_model_files.keys())
    
    crops = list(config_data['study_area']['crops_dict'].keys())
    
    """CSV per model and crop"""
    for model in climate_models:
        
        if data_source == 'delta':
            data_delta  = xr.open_dataset(config_data['cmip6_grid']['delta_per_model'].format(model=model, 
                                                                                              variable=var_input,
                                                                                              nr_months=nr_months))
            
        elif data_source == 'percentile':
            data_delta  = xr.open_dataset(config_data['cmip6_grid']['frequency_per_model'].format(model=model, 
                                                                                                  variable=var_input, 
                                                                                                  percentile=percentile,
                                                                                                  nr_months=nr_months))
       
        ds_var_list = list(data_delta.data_vars.keys())

        for crop in crops:
            list_data_crop = []
            shape_usa_file = config_data['input']['shape_crop'].format(crop=crop)
            crop_vars = [var for var in ds_var_list if crop in var]
            
            for var_name in crop_vars:
                    
                if data_source == 'delta':
                    crop, season, ssp = var_name.split('_')
                    gdf_model = fct_mean_per_county_future(data_delta, var_name, shape_usa_file)
                    final_ds = gdf_model.rename(columns={var_name:'T'})
                    final_ds['season'] = season
                
                elif data_source == 'percentile':
                    crop, _, ssp, reference = var_name.split('_')
                    gdf_model = fct_mean_per_county_future(data_delta, var_name, shape_usa_file)
                    final_ds = gdf_model.rename(columns={var_name:f'percentile_{percentile}'})
                    final_ds['reference'] = reference
            
                
                final_ds['ssp'] = ssp
                
                list_data_crop.append(final_ds)
            
            data_crop = pd.concat(list_data_crop, ignore_index=True)
            final_data_csv = data_crop.drop(columns='geometry')
            
            if data_source == 'delta':
                Path(config_data['cmip6_counties']['dir_delta'].format(crop=crop, variable=var_input, 
                                                                              nr_months=nr_months)).mkdir(parents=True, exist_ok=True)
                final_data_csv.to_csv(config_data['cmip6_counties']['delta_per_model'].format(variable=var_input, 
                                                                                                     model=model, crop=crop,
                                                                                                     nr_months=nr_months), index=False)
            elif data_source == 'percentile':
                Path(config_data['cmip6_counties']['dir_frequency'].format(crop=crop, variable=var_input, 
                                                                                  percentile=percentile, nr_months=nr_months)).mkdir(parents=True, exist_ok=True)
                final_data_csv.to_csv(config_data['cmip6_counties']['frequency_per_model'].format(variable=var_input, model=model, 
                                                                                                         crop=crop, percentile=percentile,
                                                                                                         nr_months=nr_months), index=False)


def wrap_cmip_crop(config_data, var_input, data_source, percentile=None, nr_months=0):
    climate_model_files = fct.get_model_list(config_data['input']['dir_cmip_var'].format(variable=var_input))
    climate_models = list(climate_model_files.keys())
    
    # geo_area = config_data['study_area']['geo_area']
    
    crops = list(config_data['study_area']['crops_dict'].keys())
    
    for crop in crops:
        list_data_crop = []
        for model in climate_models:
            if data_source == 'delta':
                filename= config_data['cmip6_counties']['delta_per_model'].format(nr_months=nr_months, 
                                                                                         variable=var_input, 
                                                                                         model=model, crop=crop)
                data_source_name = data_source
            elif data_source == 'percentile':
                filename = config_data['cmip6_counties']['frequency_per_model'].format(nr_months=nr_months, 
                                                                                              variable=var_input, 
                                                                                              model=model, 
                                                                                              crop=crop, 
                                                                                              percentile=percentile)
                data_source_name = data_source+'_'+str(percentile)
                
                
            df = pd.read_csv(filename)
            df['model'] = model
            list_data_crop.append(df)

        data_crop = pd.concat(list_data_crop, ignore_index=True)

        data_crop.to_csv(config_data['cmip6_counties']['csv_fut'].format(data_source=data_source_name, 
                                                                         crop=crop, nr_months=nr_months), index=False)



def plot_counties(df, year, var_name, season, crop, lon, lat):
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    df_year = df[df['year'] == year][df['season'] == season]
    df_year.plot(column=var_name, cmap='viridis', ax=ax, legend=True)
    ax.set_title(f'{var_name} in {season} {year} in counties growing {crop}')
    plt.xlim(lon[0], lon[1])
    plt.ylim(lat[0], lat[1])
    plt.show()