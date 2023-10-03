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
library(colorspace)
library(furrr)
#======================================================================
#load directories
#======================================================================
root_dir <- getwd()
dir_data <- file.path(root_dir, "data/cpc")
dir_usda <- file.path(root_dir, "data/usda-nass")
dir_model <-file.path(root_dir, "data/model_data")
#======================================================================
#load helper functions and parallel setup
#======================================================================
source("load_functions.R")
#set up parallel plan
no_cores <- availableCores() - 1
plan(multisession, workers = no_cores)
#======================================================================
#prepare shape files for quick raster processing
#======================================================================
shp_corn <-  list.files(dir_model,"us_corn_dtr_sf.rds", full.names = TRUE) %>% 
  read_rds() %>% 
  dplyr::select(geometry,zone) %>% 
  tibble() %>% 
  distinct() %>% 
  st_as_sf()

shp_soy<-  list.files(dir_model,"us_soy_dtr_sf.rds", full.names = TRUE) %>% 
  read_rds() %>% 
  dplyr::select(geometry,zone) %>% 
  tibble() %>% 
  distinct() %>% 
  st_as_sf()

shp_wheat<-  list.files(dir_model,"global_wheat_dtr_sf", full.names = TRUE) %>% 
  read_rds() %>% 
  dplyr::select(geometry,zone) %>% 
  tibble() %>% 
  distinct() %>% 
  st_as_sf()
#======================================================================
#load climate data
#======================================================================
select_vars <- c("mean_sm_spring",
                 "mean_sm_summer",
                 "mean_tmax_spring",
                 "mean_tmax_summer")
#======================================================================
#get clean climate arranged in a data frame (approx . 35 sec)
climate_data_corn <-  select_vars %>%
  future_map_dfr(
    .,
    ~ load_cpc_t(
      filename = "maize_AM_JA",
      varname = .x,
      shapefile = shp_corn
    ) %>%
      mutate(var_name = .x)
  ) %>%
  pivot_wider(names_from = var_name, values_from = value)

#save processed climate data frame
saveRDS(climate_data_corn, file.path(dir_model,"climate_sf_corn.rds"))
#======================================================================
#get clean climate arranged in a data frame (approx . 35 sec)
climate_data_soy <-  select_vars %>%
  future_map_dfr(
    .,
    ~ load_cpc_t(
      filename = "soy_AM_JA",
      varname = .x,
      shapefile = shp_soy
    ) %>%
      mutate(var_name = .x)
  ) %>%
  pivot_wider(names_from = var_name, values_from = value)

#save processed climate data frame
saveRDS(climate_data_soy, file.path(dir_model,"climate_sf_soy.rds"))
#======================================================================
#get clean climate arranged in a data frame (approx . 35 sec)
climate_data_wheat <-  select_vars %>%
  future_map_dfr(
    .,
    ~ load_cpc_t(
      filename = "wheat_MA_JJ",
      varname = .x,
      shapefile = shp_wheat
    ) %>%
      mutate(var_name = .x)
  ) %>%
  pivot_wider(names_from = var_name, values_from = value)

#save processed climate data frame
saveRDS(climate_data_wheat, file.path(dir_model,"climate_sf_wheat.rds"))
