#======================================================================
#clean work space
#======================================================================
rm(list = ls())
graphics.off()
gc()
#======================================================================
#load libraries
#======================================================================
library(tidyverse)
library(sf)
library(raster)
library(colorspace)
#======================================================================
#load directories
#======================================================================
root_dir <- getwd()
dir_model <- file.path(root_dir, "data/model_data")
#======================================================================
#load crop and climate data
#======================================================================
#load yield dataset
#[us_corn_dtr_sf.rds|us_soy_dtr_sf.rds|global_wheat_dtr_sf.rds]
crop_yield_sf <-
  list.files(dir_model, "global_wheat_dtr_sf.rds", full.names = TRUE) %>%
  readRDS() %>%
  filter(Year > 1979) %>%
  rename(crop_yield = Value) %>%
  st_as_sf()

#load climate df
#[climate_sf_corn.rds|climate_sf_soy.rds|climate_sf_wheat.rds]
climate_sf <-
  list.files(dir_model, "climate_sf_wheat.rds", full.names = TRUE)  %>%
  readRDS() %>%
  st_as_sf(x = .,
           coords = c("x", "y"),
           crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
#======================================================================
#Join climate and crop sf data sets (i.e. average raster into counties)
#======================================================================
#split climate data per year
df_point_list <- split(dplyr::select(climate_sf, -Year), 
                       climate_sf$Year)

#Split sf crop data per year
df_poly_list <- split(crop_yield_sf, crop_yield_sf$Year)

#join yield and climate sf, NAN are reproduced as climate data is filtered
full_Sf_yield_climate <- map2_dfr(
  df_poly_list,
  df_point_list,
  ~ .x %>%
    st_join(.y, left = FALSE, join = st_nearest_feature) %>%
    group_by(County, Year, State) %>%
    mutate_at(
      vars(-geometry,-zone,-n,-crop_yield,-County,-Year,-State),
      ~ (mean(., na.rm = TRUE))
    )) %>%
  group_by(County, State) %>%
  mutate(group_id = cur_group_id()) %>%
  ungroup() %>% 
  distinct()

#save model data
saveRDS(full_Sf_yield_climate, file.path(dir_model,"model_data_wheat.rds"))

