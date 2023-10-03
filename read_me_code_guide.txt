> Load_functions:
Contains general helper functions to support code processing
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
> Load_crop_data:
Process USDA sourced crop yield datasets for soy and maize.
Linear data detrend per county
Add spatial dimension to csv files
Save output as sf aware files
--
Main outputs:
*us_corn_dtr_sf.rds
*us_soy_dtr_sf.rds
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
> Load_climate_data:
Process GLEAM and CPC soil moisture and maximum temperature early (AM) and late (JA) growing season periods
Cut data based on crop yield shapefile
Save output as dataframe for easy processing into county scale in next step
--Main outputs 
*climate_sf_corn.rds
*climate_sf_soy.rds
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
> Join_crop_climate_data:
Tranform raster climate data into shapefiles based on crop data resolution
Save output as one model data frame
--Main outputs:
*model_data_soy.rds
*model_data_corn.rds
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
> Model_fitting:
Filter based on at least 30 crop yield data points per county
(optional) Standardize unit variables
Fit linear models at county scale level of the form:
detrended_yield ~ sm_spring*sm_summer + tmax_spring*tmax_summer