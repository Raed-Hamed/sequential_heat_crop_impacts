#======================================================================
#clean workspace
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
library(USAboundaries)
library(colorspace)
#======================================================================
#load directories
#======================================================================
root_dir <- getwd()
dir_data <- file.path(root_dir, "data/cpc")
dir_crop <- file.path(root_dir, "data/usda-nass")
dir_model <-file.path(root_dir, "data/model_df")
source("load_functions.R")
#======================================================================
#prepare shapefiles for quick raster processing
#======================================================================
US_shp <-us_states() %>%  us_spatial_sf %>% 
  filter(!State %in% c("ALASKA","HAWAII","PUERTO RICO")) %>% 
  dplyr::select(geometry)
#======================================================================
#load climate data
#======================================================================
#load cpc kdd spring dataset
kdd_spring <-
  load_cpc_t(filename = "kdd_gdd",
             varname = "kdd_spring",
             shapefile = US_shp)

#load cpc kdd summer dataset
kdd_summer <-   load_cpc_t(filename = "kdd_gdd",
                           varname = "kdd_summer",
                           shapefile = US_shp)

#load cpc kdd spring dataset
gdd_spring <-
  load_cpc_t(filename = "kdd_gdd",
             varname = "gdd_spring",
             shapefile = US_shp)

#load cpc kdd summer dataset
gdd_summer <-   load_cpc_t(filename = "kdd_gdd",
                           varname = "gdd_summer",
                           shapefile = US_shp)

#join all in one df
climate_data_full <- kdd_spring %>%
  rename(kdd_spring = value) %>%
  inner_join(kdd_summer) %>%
  rename(kdd_summer = value) %>% 
  inner_join(gdd_summer) %>%
  rename(gdd_summer = value) %>% 
  inner_join(gdd_spring) %>%
  rename(gdd_spring = value) 
#======================================================================
#transform df into sf object (i.e. include spatial dimension)
#======================================================================
climate_data_sf <- climate_data_full %>%
  st_as_sf(x = .,
           coords = c("x", "y"),
           crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
#======================================================================
#load yield dataset
crop_yield_sf <-list.files(dir_crop,"usda_wheat_sf", full.names = TRUE) %>%
  readRDS()%>%
  filter(n > 33) %>%
  dplyr::select(-n) %>%
  rename(crop_yield = Value)
#======================================================================
#Join climate and crop sf datasets (i.e. average raster into counties)
#======================================================================
#split climate df per year
df_point_list <- split(dplyr::select(climate_data_sf, -Year), 
                       climate_data_sf$Year)

#Split sf crop data per year
df_poly_list <- split(crop_yield_sf, crop_yield_sf$Year)

#join yield and climate sf, na are reproduced as climate_df is filtered
full_Sf_yield_climate <- map2_dfr(
  df_poly_list,
  df_point_list,
  ~ .x %>%
    st_join(.y, left = FALSE) %>%
    group_by(County, Year, State) %>%
    mutate_at(vars(-geometry,-State,-County,-Year),  ~
                (mean(., na.rm = TRUE)))
)%>%
  group_by(County, State) %>%
  mutate(group_id = cur_group_id())%>% 
  ungroup()

#save model data
saveRDS(full_Sf_yield_climate, file.path(dir_model,"model_data.rds"))

