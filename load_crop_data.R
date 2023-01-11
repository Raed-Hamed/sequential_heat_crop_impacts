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
library(USAboundaries)
library(colorspace)
#======================================================================
#load directories
#======================================================================
root_dir <- getwd()
dir_usda <- file.path(root_dir, "data/usda-nass")
dir_eurostat <- file.path(root_dir, "data/eurostat")
dir_model <- file.path(root_dir, "data/model_data")
#======================================================================
#load helper functions
#======================================================================
source("load_functions.R")
#======================================================================
#load shapefiles (US for USAboundaries)
#======================================================================
#load EU shapefile
shp_eu <- list.files(dir_eurostat, "shp", full.names = TRUE) %>%
  read_sf() %>%
  filter(LEVL_CODE %in% c(0, 1, 2, 3)) %>%
  dplyr::select(NUTS_ID, CNTR_CODE) %>%
  tibble() 
#======================================================================
#load crop yields 
#======================================================================
#load USDA wheat yield data
usda_wheat_df <- list.files(dir_usda, "usda-nass-wheat", full.names = TRUE) %>%
  map_dfr(read.csv) %>%
  dplyr::select(Year, State, County, Data.Item, Value) %>%
  filter(Data.Item == "WHEAT, WINTER - YIELD, MEASURED IN BU / ACRE") %>%
  filter(County != "OTHER (COMBINED) COUNTIES") %>% 
  mutate(Value = Value / 14.87) %>%
  dplyr::select(-Data.Item) %>% 
  arrange(Year) %>% 
  filter(Value != 0) %>% 
  group_by(State,County) %>% 
  nest() %>% 
  mutate(model_trend = map(data,~lm(Value ~ Year, data = . ))) %>% 
  mutate(augment_trend  = map(model_trend,  ~ broom::augment(.))) %>%
  unnest(cols = c(augment_trend)) %>%
  ungroup() %>%
  dplyr::select(State,County,Year,.resid) %>%
  rename(Value = .resid)

#load USDA corn yield data
usda_corn_df <- list.files(dir_usda, "usda-nass-corn", full.names = TRUE) %>%
  map_dfr(read.csv) %>%
  dplyr::select(Year, State, County, Data.Item, Value) %>%
  filter(Data.Item == "CORN, GRAIN - YIELD, MEASURED IN BU / ACRE") %>%
  filter(County != "OTHER (COMBINED) COUNTIES") %>% 
  filter(County != "OTHER COUNTIES") %>% 
  mutate(Value = Value / 14.87) %>%
  arrange(Year) %>% 
  dplyr::select(-Data.Item) %>% 
  filter(Value != 0) %>% 
  group_by(State,County) %>% 
  nest() %>% 
  mutate(model_trend = map(data,~lm(Value ~ Year, data = . ))) %>% 
  mutate(augment_trend  = map(model_trend,  ~ broom::augment(.))) %>%
  unnest(cols = c(augment_trend)) %>%
  ungroup() %>%
  dplyr::select(State,County,Year,.resid) %>%
  rename(Value = .resid)
  
#load eurostat wheat yield data
eurostat_wheat_df <- list.files(dir_eurostat, "cropstats", full.names = TRUE) %>%
  read.csv() %>%
  filter(TYPE == "Yield") %>% 
  filter(CROP_NAME == "C1100") %>% 
  dplyr::select(IDREGION,YEAR,VALUE) %>% 
  distinct() %>% 
  rename(NUTS_ID = IDREGION) %>% 
  arrange(YEAR) %>% 
  filter(VALUE != 0) %>% 
  group_by(NUTS_ID) %>% 
  nest() %>% 
  mutate(model_trend = map(data,~lm(VALUE ~ YEAR, data = . ))) %>% 
  mutate(augment_trend  = map(model_trend,  ~ broom::augment(.))) %>%
  unnest(cols = c(augment_trend)) %>%
  ungroup() %>%
  dplyr::select(NUTS_ID,YEAR,.resid) %>%
  rename(Value = .resid)
#======================================================================
#transform dataframe into SF object (i.e. add spatial dimension)
#======================================================================
#add spatial dimension to us wheat dataframe
usda_wheat_sf <- us_counties() %>%  us_spatial_sf %>%
  left_join(usda_wheat_df, by = c("County", "State")) %>%
  drop_na() %>%
  arrange(Year) %>%
  tibble() %>%
  group_by(State, County) %>%
  mutate(n = n()) %>%
  filter(n < 44) %>%
  ungroup() %>% 
  mutate(zone = "US")

#add spatial dimension to us corn dataframe
usda_corn_sf <- us_counties() %>%  us_spatial_sf %>%
  left_join(usda_corn_df, by = c("County", "State")) %>%
  drop_na() %>%
  arrange(Year) %>%
  tibble() %>%
  group_by(State, County) %>%
  mutate(n = n()) %>%
  filter(n < 44) %>%
  ungroup() %>% 
  mutate(zone = "US")

#add spatial dimension to eu dataframe
eurostat_wheat_sf <- eurostat_wheat_df %>%
  inner_join(shp_eu) %>%
  group_by(NUTS_ID) %>%
  mutate(n = n()) %>%
  rename(
    County = 1,
    Year = 2,
    State = 4,
    Value = 3
  ) %>%
  ungroup() %>% 
  mutate(zone = "EU")

#join both into one spatial informed dataframe
global_wheat_sf <- eurostat_wheat_sf %>%
  bind_rows(usda_wheat_sf) %>%
  st_as_sf()
#======================================================================
#Diagnostic plot: How many years do we have per county
#======================================================================
usda_corn_sf %>%
  dplyr::select(n, geometry, zone) %>%
  tibble() %>%
  distinct() %>%
  st_as_sf() %>%
  ggplot() +
  geom_sf(aes(fill = n)) +
  facet_grid(~zone)+
  theme_bw(base_size = 20) +
  scale_fill_continuous_sequential("Viridis")

global_wheat_sf %>%
  dplyr::select(n, geometry, zone) %>%
  tibble() %>%
  distinct() %>%
  st_as_sf() %>%
  ggplot() +
  geom_sf(aes(fill = n)) +
  facet_grid(~zone)+
  theme_bw(base_size = 20) +
  scale_fill_continuous_sequential("Viridis")
#======================================================================
#save crop yield data
#======================================================================
saveRDS(global_wheat_sf, file.path(dir_model,"global_wheat_dtr_sf.rds"))
saveRDS(usda_corn_sf, file.path(dir_model,"us_corn_dtr_sf.rds"))
